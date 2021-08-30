//
// Created by 刘凯鑫 on 2019-06-05.
//

#include "../include/graph.h"

//读图，然后计算kcore
Graph::Graph(const string graph_path) {
    Timer timer(READ_GRAPH);
    INFO("Reading graph ...");
    this->data_folder = graph_path;
    init_nm();
    adj_list = vector<vector<int>>(n, vector<int>());
    string graph_file = data_folder + FILESEP + "undirect_graph.txt";
    assert_file_exist("undirect_graph file", graph_file);
    FILE *fin = fopen(graph_file.c_str(), "r");
    int t1, t2;
    while (fscanf(fin, "%d%d", &t1, &t2) != EOF) {
        if (t1 == t2) continue;
        adj_list[t1].push_back(t2);
        //adj_list[t2].push_back(t1);
    }
    adj_list_copy=adj_list;
    result.n = this->n;
    result.m = this->m;
    cout << "init graph n: " << this->n << " m: " << this->m << endl;
}

void Graph::init_nm() {
    string attribute_file = data_folder + FILESEP + "attribute.txt";
    assert_file_exist("attribute file", attribute_file);
    ifstream attr(attribute_file);
    string line1, line2;
    char c;
    while (true) {
        attr >> c;
        if (c == '=') break;
    }
    attr >> n;
    while (true) {
        attr >> c;
        if (c == '=') break;
    }
    attr >> m;
}

//临时检测引入锚点带来的kcore
void Graph::detect_kcore() {
    Timer timer(COMPUTE_KCORE);
    result.deg_more_k = 0;
    kcore_size = this->n;
    k_tag = vector<int>(n, 1);//初始时都认为是k core中的，然后逐个检测删除
    del_nei = vector<int>(n, 0);//删除的邻居数量
    //将度小于k的初始化一下
    nodes_less_k.clear();
    for (int i = 0; i < this->n; i++) {
        if (adj_list[i].size() >= config.inputK) {
            result.deg_more_k++;
        }
        if (adj_list[i].size() - del_nei[i] < config.inputK) {
            nodes_less_k.emplace_back(i);
            k_tag[i] = -1;
        }
    }
    for (int i = 0; i < nodes_less_k.size(); i++) {
        int node = nodes_less_k[i];
        k_tag[node] = 0;
        kcore_size--;
        for (int nei: adj_list[node]) {
            if (k_tag[nei] > 0) {
                del_nei[nei]++;
                if (adj_list[nei].size() - del_nei[nei] < config.inputK) {
                    nodes_less_k.emplace_back(nei);
                    k_tag[nei] = -1;
                }
            }
        }
    }
    kc_nei = vector<int>(n, 0);//保存每个节点中拥有的k core中的邻居数量
    for (int j = 0; j < n; ++j) {
        if (k_tag[j] == 1) {
            for (int nei:adj_list[j]) {
                kc_nei[nei]++;
            }
        }
    }
    result.kcore_num = kcore_size;
    result.deg_more_k -= kcore_size;
}

void Graph::convert_to_undirected_graph(string graph_path) {
    data_folder = config.graph_location + graph_path;
    init_nm();
    vector<set<int>> undirect_adj_list(n);
    string graph_file = data_folder + FILESEP + "graph.txt";
    assert_file_exist("graph file", graph_file);
    FILE *fin = fopen(graph_file.c_str(), "r");
    char line[1000];
    long long t1, t2;
    int min_id = INT_MAX, max_id = -1, final_n = 0, final_m = 0, tmp = 0;
    bool need_index = false;
    unordered_map<long, int> str_index;
    //先读第一遍，确定最大值和最小值
    while (fgets(line, 1000, fin) != NULL) {
        if (line[0] > '9' || line[0] < '0') { continue; }
        sscanf(line, "%lld\t%lld\n", &t1, &t2);
        if (t1 >= n || t2 >= n) {
            need_index = true;
            break;
        }
    }
    fin = fopen(graph_file.c_str(), "r");
    while (fgets(line, 1000, fin) != NULL) {
        if (line[0] > '9' || line[0] < '0') { continue; }
        sscanf(line, "%lld\t%lld\n", &t1, &t2);
        min_id = min(min_id, (int) t1);
        min_id = min(min_id, (int) t2);
        max_id = max(max_id, (int) t1);
        max_id = max(max_id, (int) t2);
        if (need_index) {
            if (str_index.find(t1) == str_index.end()) {
                int index = (int) str_index.size();
                str_index[t1] = index;
            }
            if (str_index.find(t2) == str_index.end()) {
                int index = (int) str_index.size();
                str_index[t2] = index;
            }
            t1 = str_index[t1];
            t2 = str_index[t2];
        }
        undirect_adj_list[t1].insert(t2);
        undirect_adj_list[t2].insert(t1);
    }
    string undirect_graph_file = data_folder + FILESEP + "undirect_graph.txt";
    FILE *fout = fopen(undirect_graph_file.c_str(), "w");
    for (int j = 0; j < undirect_adj_list.size(); ++j) {
        if (!undirect_adj_list[j].empty()) {
            final_n++;
            for (int des:undirect_adj_list[j]) {
                final_m++;
                fprintf(fout, "%d\t%d\n", j, des);
            }
        }
    }
    fclose(fout);
    string undirect_graph_attr_file = data_folder + FILESEP + "undirect_graph_attribute.txt";
    fout = fopen(undirect_graph_attr_file.c_str(), "w");
    fprintf(fout, "n=%d\nm=%lld\n", n, m);
    fprintf(fout, "final_n=\t%d\nfinal_m=\t%d\n", final_n, final_m);
    fprintf(fout, "min_node=\t%d\nmax_node=\t%d\n", min_id, max_id);
    for (auto item:str_index) {
        fprintf(fout, "%ld\t%d\n", item.first, item.second);
    }
    fclose(fin);
    fclose(fout);
}

void Graph::core_decompsition() {
    core_num = vector<int>(n, 0);
    vector<int> deg = vector<int>(n, 0);
    vector<bool> del_flag = vector<bool>(n, 0);
    unordered_map<int, vector<int>> active_nodes;
    for (int i = 0; i < n; ++i) {
        deg[i] = adj_list[i].size();
        active_nodes[deg[i]].emplace_back(i);
    }
    for (int k = 0; k < n; ++k) {
        while (!active_nodes[k].empty()) {
            int seed = active_nodes[k].back();
            active_nodes[k].pop_back();
            if (k == deg[seed]) {
                core_num[seed] = k;
                for (int nei:adj_list[seed]) {
                    if (!del_flag[nei] && deg[nei] > k) {
                        deg[nei]--;
                        active_nodes[deg[nei]].emplace_back(nei);
                    }
                }
                del_flag[seed] = true;
            }
        }
    }
    //check
    for (int i = 0; i < n; ++i) {
        if (k_tag[i] == true) {
            assert(core_num[i] >= config.inputK);
        }
    }
}

int Graph::detect_kcore_with_anc(vector<int> &b) {
    Timer timer(COMPUTE_KCORE_WITH_ANC);
    int kcore_size_tmp = int (kmb_core_size);
    vector<int> k_tag_tmp1=k_tag_tmp;
    del_nei = vector<int>(n, 0);//删除的邻居数量
    for (int anc:b){
        del_nei[anc] = -1*n;
        k_tag_tmp1[anc]=2;
    }
    vector<int> nodes_less_k_tmp=nodes_less_k;
    for (int i = 0; i < nodes_less_k_tmp.size(); ++i) {
        //nodes_less_k_tmp中有锚点
        int node=nodes_less_k_tmp[i];
        if (k_tag_tmp1[node]!=2){
            k_tag_tmp1[node]=0;
            kcore_size_tmp--;
            for (int nei: adj_list[node]){
                if(k_tag_tmp1[nei]>0){
                    del_nei[nei]++;
                    if (adj_list[nei].size()-del_nei[nei]<config.inputK){
                        assert(core_num[nei]<config.inputK);
                        nodes_less_k_tmp.emplace_back(nei);
                        k_tag_tmp1[nei]=-1;
                    }
                }
            }
        }
    }
    return kcore_size_tmp;
}

void Graph::generate_subgraph_rand() {
    //
    int chosen_nodes_num=0;
    vector<bool> chosen_nodes(n,0);
    while (chosen_nodes_num < config.new_graph_size) {
        unsigned int x = n;
        while (x >= n) {
            x = rand() / ((RAND_MAX + 1u) / n);
        }
        if (!chosen_nodes[x]){
            chosen_nodes[x]=true;
            chosen_nodes_num++;
        }
    }
    m=0;
    vector<int> result_nodes;
    // update adj_list
    for (int i = 0; i < n; ++i) {
        if (chosen_nodes[i]){
            result_nodes.emplace_back(i);
            for (int j = 0; j < adj_list[i].size(); ++j) {
                int nei = adj_list[i][j];
                if (!chosen_nodes[nei]){
                    adj_list[i][j]=adj_list[i].back();
                    adj_list[i].pop_back();
                }else{
                    m++;
                }
            }
        } else {
            adj_list[i].clear();
        }
    }
    m/=2;
    INFO(config.new_graph_size,m);
}

void Graph::generate_subgraph_bfs(){
    adj_list=adj_list_copy;
    // 第一步选点后BFS，如果选到孤立点，则重新再选
    vector<int> chosen_nodes_queue;
    vector<bool> chosen_nodes;
    while (chosen_nodes_queue.size()<config.new_graph_size){
        unsigned int seed = n;
        while (seed>=n){
            seed = rand()/((RAND_MAX + 1u)/n);
        }
        vector<bool> tmp(n,0);
        swap(chosen_nodes,tmp);
        chosen_nodes[seed]=true;
        chosen_nodes_queue.emplace_back((int)seed);
        int active_pointer=0;
        while (chosen_nodes_queue.size()<config.new_graph_size && active_pointer<chosen_nodes_queue.size()){
            int active_node=chosen_nodes_queue[active_pointer++];
            unsigned time_seed = std::chrono::system_clock::now ().time_since_epoch ().count ();
            shuffle(adj_list[active_node].begin(),adj_list[active_node].end(),std::default_random_engine (time_seed));
            //遍历邻居
            for (int nei: adj_list[active_node]){
                if (chosen_nodes_queue.size()<config.new_graph_size && !chosen_nodes[nei]){
                    chosen_nodes_queue.emplace_back(nei);
                    chosen_nodes[nei]=true;
                }
            }
        }
    }

    m=0;
    vector<int> result_nodes;
    // update adj_list
    for (int i = 0; i < n; ++i) {
        if (chosen_nodes[i]){
            result_nodes.emplace_back(i);
            for (int j = 0; j < adj_list[i].size(); ++j) {
                int nei = adj_list[i][j];
                if (!chosen_nodes[nei]){
                    adj_list[i][j--]=adj_list[i].back();
                    adj_list[i].pop_back();
                }else{
                    m++;
                }
            }
        } else {
            adj_list[i].clear();
        }
    }
    m/=2;
    INFO(config.new_graph_size,m);
}

void Graph::save_subgraph(string outputFile){
    //重新编号
    int current=0;
    vector<int> new_id(n,-1);
    for (int i = 0; i < n; ++i) {
        if (!adj_list[i].empty()){
            new_id[i]=current++;
        }
    }
    //
    vector<vector<int>> sub_adj_list=vector<vector<int>> (current,vector<int>());
    for (int i = 0; i < n; ++i) {
        if (!adj_list[i].empty()){
            int nid=new_id[i];
            sub_adj_list[nid]=adj_list[i];
            for (int j = 0; j < sub_adj_list[nid].size(); ++j) {
                sub_adj_list[nid][j]=new_id[sub_adj_list[nid][j]];
            }
            sort(sub_adj_list[nid].begin(),sub_adj_list[nid].end());
        }
    }
    FILE *fout;
    string graph_path=outputFile+"/undirect_graph.txt";
    if ((fout = fopen(graph_path.c_str(), "w+")) != NULL) {
        INFO("Save graph to file ", outputFile);
        for (int i = 0; i < sub_adj_list.size(); ++i) {
            for (int j = 0; j < sub_adj_list[i].size(); ++j) {
                fprintf(fout, "%d\t%d\n", i, sub_adj_list[i][j]);
                assert(sub_adj_list[i][j]<config.new_graph_size);
            }
        }
        fclose(fout);
        INFO("save graph finished");
    } else {
        INFO("Can not find output file, please test path!");
    }
    string attribute_path=outputFile+"/attribute.txt";
    if ((fout = fopen(attribute_path.c_str(), "w+")) != NULL) {
        INFO("Save graph to file ", outputFile);
        fprintf(fout,"n=%d\n",config.new_graph_size);
        fprintf(fout,"m=%d\n",m);
        INFO("save graph finished");
    } else {
        INFO("Can not find output file, please test path!");
    }
}
/*


    unsigned int seed = n;
    while (seed >= n) {
        seed = rand() / ((RAND_MAX + 1u) / n);
    }
    // BFS 直到得到子图大小符合
    vector<bool> chosen_nodes(n,0);
    vector<int> chosen_nodes_queue{int (seed)};
    chosen_nodes[seed]= true;
    int active_pointer=0;
    while (chosen_nodes_queue.size()<config.new_graph_size){
        int active_node=chosen_nodes_queue[active_pointer++];
        unsigned time_seed = std::chrono::system_clock::now ().time_since_epoch ().count ();
        shuffle(adj_list[active_node].begin(),adj_list[active_node].end(),std::default_random_engine (time_seed));
        //遍历邻居
        for (int nei: adj_list[active_node]){
            if (chosen_nodes_queue.size()<config.new_graph_size && !chosen_nodes[nei]){
                chosen_nodes_queue.emplace_back(nei);
                chosen_nodes[nei]=true;
            }
        }
    }
 */