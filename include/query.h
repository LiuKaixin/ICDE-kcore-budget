//
// Created by 刘凯鑫 on 2020/1/20.
//

#ifndef GETUNDIRECTGRAPH_QUERY_H
#define GETUNDIRECTGRAPH_QUERY_H

#include "myalgo.h"
#include "anchoredkcore.h"

//exact
void combination(Graph &graph, const vector<int> &a, vector<int> b, int n, int m, const int &M, int &max_fols_and_ancs,
                 vector<int> &final_b) {
    if (m > 0) {
        for (int i = n; i >= m; i--) {
            b[m - 1] = a[i - 1];
            //由于确定了当前位置的值
            //则下次递归c(n-1,m-1)
            combination(graph, a, b, i - 1, m - 1, M, max_fols_and_ancs, final_b);
        }
    } else {
        //对b中的点锚定，然后计算kcore，得到
        int followers_and_ancs = graph.detect_kcore_with_anc(b) - graph.kcore_size;
        if (followers_and_ancs > max_fols_and_ancs) {
            max_fols_and_ancs = followers_and_ancs;
            final_b = b;
        }
    }
}

pair<double, double> exact_query(Graph &graph) {
    graph.core_decompsition();
    double start_time = (double) clock() / CLOCKS_PER_SEC;
    //kmq shell及其邻居都有可能是！
    Timer timer1(EXACT_QUERY);
    pair<double, double> ans;
    vector<bool> candidate_nodes(graph.n, 0);
    graph.k_tag_tmp = vector<int>(graph.n, 1);
    graph.nodes_less_k.clear();
    graph.kmb_core_size = graph.n;
    for (int i = 0; i < graph.n; ++i) {
        //如果属于kmb shell，及其邻居可以作为候选点。
        //if (graph.core_num[i] >= config.inputK - config.inputB && graph.core_num[i] < config.inputK) {
        if (graph.core_num[i] >= config.inputK - config.inputB && graph.core_num[i] < config.inputK) {
            candidate_nodes[i] = true;
            for (int nei:graph.adj_list[i]) {
                if (graph.core_num[nei] < config.inputK - config.inputB) {
                    candidate_nodes[nei] = true;
                }
            }
        }
    }
    //删除点，对于所有core number < k，且不在候选中的
    for (int i = 0; i < graph.n; ++i) {
        if (!candidate_nodes[i] && graph.core_num[i] < config.inputK) {
            graph.k_tag_tmp[i] = -1;
            graph.kmb_core_size--;
            for (int nei:graph.adj_list[i]) {
                for (int j = 0; j < graph.adj_list[nei].size(); ++j) {
                    if (graph.adj_list[nei][j] == i) {
                        graph.adj_list[nei].erase(graph.adj_list[nei].begin() + j);
                    }
                }
            }
            graph.adj_list[i].clear();
        }
    }
    vector<int> a, b(config.inputB);
    for (int i = 0; i < graph.n; ++i) {
        if (candidate_nodes[i]) {
            if (graph.adj_list[i].size() < config.inputK){
                graph.nodes_less_k.emplace_back(i);
                graph.k_tag_tmp[i] = -1;
            }
            a.emplace_back(i);
        }
    }
    int max_fols_and_ancs = -1;
    vector<int> final_b;
    INFO(a.size());
    /*
    for (int i = 0; i < graph.n; ++i) {
        if(!graph.adj_list[i].empty()){
            INFO(i,graph.k_tag[i],candidate_nodes[i]);
            assert(graph.k_tag[i]!=candidate_nodes[i]);
        }
    }
    for (int i = 0; i < graph.n; ++i) {
        if (candidate_nodes[i]) {
            INFO(i,graph.k_tag[i],graph.core_num[i]);
        }
    }*/
    if (a.size()<config.inputB){
        max_fols_and_ancs=a.size();
    }else{
        combination(graph, a, b, a.size(), config.inputB, config.inputB, max_fols_and_ancs, final_b);
        INFO(final_b);
    }
    INFO(max_fols_and_ancs);
    double run_time = (double) clock() / CLOCKS_PER_SEC - start_time;
    INFO(run_time);
    result.follower_num = max_fols_and_ancs;
    return make_pair(double(max_fols_and_ancs), run_time);
}
/*
pair<double, double> multi_thread_exact(Graph &graph){

    unsigned NUM_CORES = std::thread::hardware_concurrency()-4;
    assert(NUM_CORES >= 2);

    int num_thread = min(10, int(NUM_CORES));


    int avg_queries_per_thread = 20/num_thread;

    vector<vector<int>> source_for_all_core(num_thread);
    vector<unordered_map<int, vector<pair<int ,double>>>> ppv_for_all_core(num_thread);

    for(int tid=0; tid<num_thread; tid++){
        int s = tid*avg_queries_per_thread;
        int t = s+avg_queries_per_thread;

        if(tid==num_thread-1)
            t+=query_size%num_thread;

        for(;s<t;s++){
            // cout << s+1 <<". source node:" << queries[s] << endl;
            source_for_all_core[tid].push_back(queries[s]);
        }
    }
    {
        INFO("Computing Exact results...");
        std::vector< std::future<void> > futures(num_thread);
        for(int tid=0; tid<num_thread; tid++){
            futures[tid] = std::async( std::launch::async, exact_query, std::ref(graph)) ;
        }
        std::for_each( futures.begin(), futures.end(), std::mem_fn(&std::future<void>::wait));
    }
}
*/
//LKX
void lkx(Graph &graph) {
    //首先取得非kcore中度大于等于k的点及其邻居
    Timer timer1(LKX_QUERY);
    SubGraph subGraph(graph);
    return subGraph.main_algorithm();
}

void set_result(Graph &graph, int used_counter) {
    //给定锚点后，重新计算一遍followers
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - result.startTime).count();
    result.whole_time_usage = duration / TIMES_PER_SEC;;
    result.total_mem_usage = get_proc_memory() / 1000.0;
    result.total_time_usage = Timer::used(used_counter);
    result.kcore_num = graph.kcore_size;
    vector<int> k_tag_old = graph.k_tag;
    if (config.algo != RAND && config.algo != EXACT) {
        //result.follower_num = graph.detect_kcore(anchors, true) - result.kcore_num - config.inputB;
        result.followers.clear();
        for (int j = 0; j < graph.n; ++j) {
            if (graph.k_tag[j] == 1 && k_tag_old[j] == 0) {
                result.followers.emplace_back(j);
            }
        }
        INFO(result.follower_num);
    }
    if (used_counter == LKX_QUERY) {
        result.merge_time = Timer::used(MERGE_STEP);
        result.merge_time_ratio = Timer::used(MERGE_STEP) * 100 / Timer::used(used_counter);
        result.supporter_time = Timer::used(SUPPORTER_COMPUTE);
        result.supporter_time_ratio = Timer::used(SUPPORTER_COMPUTE) * 100 / Timer::used(used_counter);
    }
}

void dataOutput() {
    //write
    char record[100];
    FILE *fs;
    if (config.algo == LKX) {
        fs = fopen(config.outfile.c_str(), "a");
    } else if (config.algo == EXACT) {
        fs = fopen(config.exactfile.c_str(), "a");
    }

    char fsch;

    if (fs == NULL) {
        printf("ERROR!");
        exit(1);
    } else {
        sprintf(record, "%s", config.graph_alias.c_str());
        fwrite(record, sizeof(*record), strlen(record), fs);
        fsch = putc('\t', fs);
        sprintf(record, "k%d", config.inputK);
        fwrite(record, sizeof(*record), strlen(record), fs);
        fsch = putc('\t', fs);
        sprintf(record, "inputB%d", config.inputB);
        fwrite(record, sizeof(*record), strlen(record), fs);
        fsch = putc('\t', fs);
        sprintf(record, "quota%.0f", result.follower_num);
        fwrite(record, sizeof(*record), strlen(record), fs);
        fsch = putc('\t', fs);
        sprintf(record, "t%.1f", config.threshold);
        fwrite(record, sizeof(*record), strlen(record), fs);
        fsch = putc('\t', fs);
        sprintf(record, "%.3lfs", result.whole_time_usage);
        fwrite(record, sizeof(*record), strlen(record), fs);
        fsch = putc('\t', fs);
        sprintf(record, "deg_more_k%d", result.deg_more_k);
        fwrite(record, sizeof(*record), strlen(record), fs);
        fsch = putc('\t', fs);
        sprintf(record, "kc%d", result.kcore_num);
        fwrite(record, sizeof(*record), strlen(record), fs);
        fsch = putc('\t', fs);
        sprintf(record, "mode%d", config.mode);
        fwrite(record, sizeof(*record), strlen(record), fs);
        fsch = putc('\t', fs);
        sprintf(record, "supporters%d", config.supporter_num);
        fwrite(record, sizeof(*record), strlen(record), fs);
        fsch = putc('\t', fs);
        sprintf(record, "neighbor_num%d", config.neighbor_num);
        fwrite(record, sizeof(*record), strlen(record), fs);
        fsch = putc('\t', fs);

        for (double i = 0.1; i < 1; i += 0.1) {
            int key = int(i * 10 + 0.1);
            if (result.fol_anc_time.find(key) != result.fol_anc_time.end()) {
                sprintf(record, "fol:%.1f\tanc:%.2f\ttime:%f\t", i, result.fol_anc_time[key].first,
                        result.fol_anc_time[key].second);
                fwrite(record, sizeof(*record), strlen(record), fs);
                fsch = putc('\t', fs);
            }
        }
        fsch = putc('\n', fs);
        fsch = putc('\n', fs);
    }
}

void query(Graph &graph) {
    INFO(graph.kcore_size);
    INFO(result.deg_more_k);
    result.startTime = std::chrono::steady_clock::now();
    INFO(config.inputK, config.inputB, config.threshold);
    int used_counter = 0;
    double anc_num = -1;
    if (config.algo == LKX) {
        used_counter = LKX_QUERY;
        lkx(graph);
    } else if (config.algo == EXACT) {
        used_counter = EXACT_QUERY;
        pair<double, double> ans = exact_query(graph);
        result.fol_anc_time[config.inputB] = ans;
    }

    set_result(graph, used_counter);
    dataOutput();
}

#endif //GETUNDIRECTGRAPH_QUERY_H s_v(node_group[anc_node])=-316459 -315867 -315247 -314695 -313225 -312691 -310284 -307167 -306966 -259525 -242191 -236627 255679 273754 312709 316661
