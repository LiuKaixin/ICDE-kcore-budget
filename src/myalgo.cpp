//
// Created by 刘凯鑫 on 2020/1/24.
//

#include "../include/myalgo.h"
//template<typename T>
bool cmp_pair(const pair<int, double> &p1, const pair<int, double> &p2) {
    if (p1.second < p2.second) return true;
    if (p1.second > p2.second) return false;
    return p1.first < p2.first;
    //return p1.second < p2.second; //升序排列
}

//先按照second降序排列，相等时按照first升序排列
bool cmp_pair_dec(const pair<int, int> &p1, const pair<int, int> &p2) {
    if (p1.second > p2.second) return true;
    if (p1.second < p2.second) return false;
    return p1.first < p2.first;
}

bool compare_state(vector<pair<int, int>> &state_vec, int node, int function_counter, int ag_id) {
    if (state_vec[node].first == function_counter && state_vec[node].second == ag_id) {
        return true;
    } else {
        state_vec[node].first = function_counter;
        state_vec[node].second = ag_id;
        return false;
    }
}
bool cmp_fol_deg_merge(const fol_deg_merge &a, const fol_deg_merge &b) {
    if (a.fol_deg > b.fol_deg) {
        return true;
    }
    if (a.fol_deg == b.fol_deg && a.score < b.score) {
        return true;
    }
    return false;
    //return (a.fol_deg>=b.fol_deg)&&(a.score<b.score);
}

vector<int> s_v(const set<int> & tmp){
    vector<int> tmp_v;
    tmp_v.assign(tmp.begin(),tmp.end());
    return tmp_v;
}
unordered_map<int, int> get_different_nodes(const vector<int> &vec1, const vector<int> &vec2){
    unordered_map<int, int> tmp;
    vector<int> res;
    for (int node1:vec1) {
        tmp[node1] = 1;
    }
    for (int node2:vec2) {
        tmp[node2] -= 1;
    }
    return tmp;
}
vector<int> vec1_vec2(const vector<int> &vec1, const vector<int> &vec2) {
    unordered_map<int, int> tmp;
    vector<int> res;
    if (vec1.empty()){
        return res;
    }
    for (int node1:vec1) {
        tmp[node1] = 1;
    }
    for (int node2:vec2) {
        tmp[node2] -= 1;
    }
    for(int fol:vec1){
        if (tmp[fol]==1){
            res.emplace_back(fol);
        }
    }
    return res;
}


vector<int> merge_vec(const vector<int> &vec1, const vector<int> &vec2, vector<int> &vec_erase) {
    vector<int> res = vec1;
    res.insert(res.end(), vec2.begin(), vec2.end());
    sort(res.begin(), res.end());
    auto pos = unique(res.begin(), res.end());
    vec_erase.assign(pos, res.end());
    res.erase(pos, res.end());
    return res;
}

vector<int> add_vec(vector<int> &vec1, const vector<int> &vec2) {
    vec1.insert(vec1.end(), vec2.begin(), vec2.end());
    sort(vec1.begin(), vec1.end());
    auto pos = unique(vec1.begin(), vec1.end());
    vector<int> vec_erase(pos, vec1.end());
    vec1.erase(pos, vec1.end());
    return vec_erase;
}

void add_score_to_nodes(vector<int> anchores, unordered_map<int, double> &node_score, double score) {
    for (auto node:anchores) {
        if (node_score.find(node) == node_score.end()) {
            node_score[node] = score;
        } else {
            node_score[node] += score;
        }
    }
}

void add_score_to_nodes(set<int> anchores, unordered_map<int, double> &node_score, double score) {
    for (auto node:anchores) {
        if (node_score.find(node) == node_score.end()) {
            node_score[node] = score;
        } else {
            node_score[node] += score;
        }
    }
}

queue<int> push_neighbor_queue(const vector<int> &neighbors) {
    queue<int> nei_q;
    for (int i = 0; i < neighbors.size(); i++) {
        nei_q.push(neighbors[i]);
    }
    return nei_q;
}

queue<int> push_neighbor_queue(set<int> neighbors) {
    queue<int> nei_q;
    for (int node:neighbors) {
        nei_q.push(node);
    }
    return nei_q;
}

int compare_double(const double &score1, const double &score2,double thre) {
    if ((score1 - score2) > thre) {
        return 1;
    } else if ((score2 - score1) > thre) {
        return -1;
    } else {
        return 0;
    }
}
/*
set<int> del_nodes_safely(set<int> deleted_vec, unordered_map<int, vector<int>> &g) {
    queue<int> active_node = push_neighbor_queue(deleted_vec);
    set<int> del_node_set(deleted_vec.begin(), deleted_vec.end());
    while (!active_node.empty()) {
        int node = active_node.front();
        active_node.pop();
        vector<int> neighbors = g.at(node);
        del_node(node, g);
        for (int nei:neighbors) {
            //判断node_kcoreNei是否大于0，来避免删除锚点
            int nei_num = g.at(nei).size();
            if (get_value(NODE_KCORENEI, nei) >= 0 && nei_num + get_value(NODE_KCORENEI, nei) == config.inputK - 1 &&
                del_node_set.find(nei) == del_node_set.end()) {
                active_node.push(nei);
                del_node_set.emplace(nei);
            }
        }
    }
    return del_node_set;
}
*/
vector<int> del_node(int deleted_node, unordered_map<int, vector<int>> &g) {
    vector<int> neighbors;
    if (g.find(deleted_node) != g.end()) {
        neighbors = g[deleted_node];
        for (int i = 0; i < g[deleted_node].size(); i++) {
            int neighbor = g[deleted_node][i];
            for (int j = 0; j < g[neighbor].size(); j++) {
                if (g[neighbor][j] == deleted_node) {
                    g[neighbor].erase(g[neighbor].begin() + j);
                    break;
                }
            }
        }
        g.erase(deleted_node);
    }
    return neighbors;
}
/*
bool insert_node(int ins_node, const vector<int> &neighbors, unordered_map<int, vector<int>> &g) {
    if (g.find(ins_node) != g.end())
        return false;
    vector<int> legal_neighbors;
    for (int neighbor : neighbors) {
        if (g.find(neighbor) != g.end()) {
            g[neighbor].push_back(ins_node);
            legal_neighbors.push_back(neighbor);
        }
    }
    g[ins_node] = legal_neighbors;
    return true;
}
*/
//得到去掉k-core的子图
unordered_map<int, vector<int>> get_margin_subgraph(const vector<vector<int>> &graph, const vector<int> &kTag) {
    unordered_map<int, vector<int>> g;
    for (int i = graph.size() - 1; i >= 0; i--) {
        if (kTag[i] == 0) {
            vector<int> neighbors;
            for (int j = 0; j < graph[i].size(); j++) {
                if (kTag[graph[i][j]] == 0) {
                    neighbors.push_back(graph[i][j]);
                }
            }
            g[i] = neighbors;
        }
    }
    return g;
}
//持续删除度在[lbound,ubound)中的点

vector<int>
continuous_del_branch_node(unordered_map<int, vector<int>> &g, queue<int> active_nodes,
                           const unordered_map<int, int> &node_kcoreNei, int ubound) {
    vector<int> deleted_nodes;
    while (!active_nodes.empty()) {
        int node = active_nodes.front();
        active_nodes.pop();
        vector<int> neighbors = g[node];
        del_node(node, g);
        for (int nei:neighbors) {
            //判断node_kcoreNei是否大于0，来避免删除锚点
            if (get_value(node_kcoreNei, nei) >= 0 && g[nei].size() + get_value(node_kcoreNei, nei) == ubound - 1) {
                active_nodes.push(nei);
                deleted_nodes.push_back(nei);
            }
        }
    }
    return deleted_nodes;
}

vector<int> continuous_del_branch_node(unordered_map<int, vector<int>> &g, queue<int> active_nodes,
                                       const unordered_map<int, int> &node_kcoreNei,
                                       unordered_map<int, vector<int>> &del_adj_list, int ubound) {
    vector<int> deleted_nodes;
    while (!active_nodes.empty()) {
        int node = active_nodes.front();
        active_nodes.pop();
        vector<int> neighbors = g[node];
        del_adj_list[node] = del_node(node, g);
        for (int nei:neighbors) {
            //判断node_kcoreNei是否大于0，来避免删除锚点
            if (get_value(node_kcoreNei, nei) >= 0 && g[nei].size() + get_value(node_kcoreNei, nei) == ubound - 1) {
                active_nodes.push(nei);
                deleted_nodes.push_back(nei);
            }
        }
    }
    return deleted_nodes;
}

///去掉度小于K的点
unordered_map<int, vector<int>>
refine_margin_graph(const unordered_map<int, vector<int>> &graph, const unordered_map<int, int> &node_kcoreNei) {
    unordered_map<int, vector<int>> g = graph;
    queue<int> active_nodes;
    for (auto v: graph) {
        //如果度不满足条件，则删除该点，并删除和邻居的边
        if (get_value(node_kcoreNei, v.first) >= 0 &&
            v.second.size() + get_value(node_kcoreNei, v.first) < config.inputK) {
            del_node(v.first, g);
        }
    }
    return g;
}

///////将图转为多个子图存储
vector<unordered_map<int, vector<int>>> get_disjoint_subgraphs(const unordered_map<int, vector<int>> &graph) {
    vector<unordered_map<int, vector<int>>> subgraphs;
    unordered_map<int, bool> detected_node;
    for (auto v: graph) {
        unordered_map<int, vector<int>> subg;
        if (detected_node.find(v.first) != detected_node.end()) {
            continue;
        }
        queue<int> nei_q;
        nei_q.emplace(v.first);
        while (!nei_q.empty()) {
            int nei = nei_q.front();
            detected_node[nei] = true;
            nei_q.pop();
            subg[nei] = graph.at(nei);
            for (int node:graph.at(nei)) {
                if (detected_node.find(node) == detected_node.end()) {
                    nei_q.push(node);
                }
            }
        }
        if (subg.size() > 0)
            subgraphs.push_back(subg);
    }
    return subgraphs;
}

///为每个子图计算相应的node_kcoreNei
unordered_map<int, int>
get_subgraphs_node_kcoreNei(unordered_map<int, vector<int>> subgraph, unordered_map<int, int> node_kcoreNei) {
    unordered_map<int, int> subgraph_node_kcoreNei;
    for (auto v:subgraph) {
        subgraph_node_kcoreNei[v.first] = node_kcoreNei[v.first];
    }
    return subgraph_node_kcoreNei;
}

int vec_find(const vector<int> &vec, int node) {
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i] == node)
            return i;
    }
    return -1;
}

bool vec_erase(vector<int> &vec, int index) {
    if (index >= 0 & index < vec.size()) {
        vec[index] = vec.back();
        vec.pop_back();
        return true;
    }
    return false;
}

/*
 * 1. 对refine函数的优化
    1. 贪心选择
        1. 每次找度最大的
        2. 每次找节点最容易变成kcore的
    2. 设计剪枝
        1. 如果剩下的锚点数/邻居最大的度数量都
        2.
2. 并行处理
3. profile 进一步优化

void compute_degree_candidate_anchores_sets(vector<vector<int>> &result, int &min_size, vector<int> &alive_raw_anchors,
                                            vector<int> new_anchores,
                                            unordered_map<int, vector<int>> &node_nei,
                                            unordered_map<int, vector<int>> &nei_node,
                                            unordered_map<int, vector<int>> subg,
                                            const unordered_map<int, int> &node_kcoreNei,
                                            const unordered_map<int, vector<int>> &g_margin) {
    //找度最大的，注意点可能有多个
    int max_degree = 0, max_score = 0;
    vector<int> nei_id_set;
    set<int> max_node_set;
    for (auto test_node:alive_raw_anchors) {
        int kcoreNei = node_kcoreNei.find(test_node) == node_kcoreNei.end() ? 0
                                                                            : node_kcoreNei.find(test_node)->second;
        //除了分数还要求，有共享的邻居
        for (auto nei: node_nei[test_node]) {
            if (nei_node.find(nei) != nei_node.end()) {
                if (int((subg[test_node].size() + kcoreNei)) > max_score) {
                    max_score = subg[test_node].size() + kcoreNei;
                    max_node_set.clear();
                }
                if (subg[test_node].size() + kcoreNei == max_score) {
                    max_node_set.emplace(test_node);
                }
                break;
            }
        }

    }
    if (config.inputK - max_score > threshold_compute_candidate_sets) {
        return;
    }
    //对所有nei找到拥有最多最大点邻居的，这里可以再优化，不对所有的nei而是对有可能的那些
    for (auto iter:nei_node) {
        int nei_max_node = 0;
        for (int node:iter.second) {
            if (max_node_set.find(node) != max_node_set.end()) {
                ++nei_max_node;
            }
        }
        if (nei_max_node > max_degree) {
            nei_id_set.clear();
            max_degree = nei_max_node;
        }
        if (nei_max_node == max_degree) {
            nei_id_set.emplace_back(iter.first);
        }
    }
    unordered_map<int, vector<int>> nei_node_copy = nei_node, node_nei_copy = node_nei;
    for (int nei_id:nei_id_set) {
        // alive_raw_anchors
        vector<int> alive_raw_anchors_copy = alive_raw_anchors;
        insert_node(nei_id, g_margin.find(nei_id)->second, subg);
        for (int test_node: nei_node.find(nei_id)->second) {
            int kcoreNei = node_kcoreNei.find(test_node) == node_kcoreNei.end() ? 0
                                                                                : node_kcoreNei.find(test_node)->second;
            if (subg[test_node].size() + kcoreNei == config.inputK) {
                vec_erase(alive_raw_anchors, vec_find(alive_raw_anchors, test_node));
            }
        }
        // new_anchores
        new_anchores.emplace_back(nei_id);
        if (alive_raw_anchors.size() + new_anchores.size() < min_size) {
            min_size = alive_raw_anchors.size() + new_anchores.size();
            result.clear();
        }
        if (alive_raw_anchors.size() + new_anchores.size() == min_size) {
            vector<int> candidate_result(new_anchores);
            candidate_result.insert(candidate_result.end(), alive_raw_anchors.begin(), alive_raw_anchors.end());
            result.emplace_back(candidate_result);
        }
        if (!alive_raw_anchors.empty()) {
            for (int affected_node:nei_node.at(nei_id)) {
                vec_erase(node_nei.at(affected_node), vec_find(node_nei.at(affected_node), nei_id));
            }
            nei_node.erase(nei_id);

            compute_degree_candidate_anchores_sets(result, min_size, alive_raw_anchors, new_anchores, node_nei,
                                                   nei_node, subg,
                                                   node_kcoreNei, g_margin);
            alive_raw_anchors = alive_raw_anchors_copy;
            new_anchores.pop_back();
        }
        del_node(nei_id, subg);
    }
    node_nei = node_nei_copy;
    nei_node = nei_node_copy;
}
 */
void saveGraph(unordered_map<int, vector<int>> g, string outputFile) {
    FILE *fout;
    if ((fout = fopen(outputFile.c_str(), "w+")) != NULL) {
        INFO("Save graph to file ", outputFile);
        long long m = 0;
        int t1, t2;
        for (auto v:g) {
            for (auto nei:v.second) {
                fprintf(fout, "%d\t%d\n", v.first, nei);
            }
        }
        fclose(fout);
        INFO("save graph finished");
    } else {
        INFO("Can not find output file, please test path!");
    }
}

