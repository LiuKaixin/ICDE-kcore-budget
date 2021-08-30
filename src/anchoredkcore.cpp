//
// Created by 刘凯鑫 on 2020/1/24.
//

/// 下一步要做的是只有当前点集，进行融合！！！

#include "../include/anchoredkcore.h"

SubGraph::SubGraph(Graph &graph, bool merge, bool divide) {
    //初始化
    fol_thre = result.deg_more_k * 0.9;
    thre_detect = config.threshold;
    config.supporter_num = ceil((max_supporters + 1) * thre_detect);
    config.neighbor_num = ceil((max_supporters + 1) * thre_detect);

    kc_nei = graph.kc_nei; // 每个点的邻居中属于kcore的数量
    k_tag_anchor = graph.k_tag; // 是否属于kcore
    whole_anchor = 0;
    node_group = vector<set<int>>(graph.n, set<int>{}); // 节点所在的锚组，如果是锚点就只有一个锚组，如果是支持者则可能有多个锚组
    nei_score_state = vector<long long>(graph.n); //两个标志位，似乎用于对锚点进行探索
    sub_nei_score_state = vector<long long>(graph.n);
    del_nodes_state = vector<long long>(graph.n); //判断节点是否被删除

    neighbor_score = vector<int>(graph.n);  //似乎是为了计算supporter时，给的得分
    sub_neighbor_score = vector<int>(graph.n);
    delete_nodes = vector<bool>(graph.n);
    need_update_nodes = vector<int>(graph.m, -1);
    //计算度不小于K的非kcore及其邻居
    remove_core_adj = vector<vector<int>>(graph.n, vector<int>());
    remove_core_adj_copy = vector<vector<int >>(graph.n, vector<int>());

    supporter_record = vector<int>(graph.n);
    new_supporter_record = vector<int>(graph.n);
    supporter_order = vector<int>(graph.n);
    new_supporter_order = vector<int>(graph.n);

    int counter1 = 0, counter2 = 0;
    //去除掉kcore中的点。remove_core_adj中只包含度大于k的点，remove_core_adj_copy包含所有点
    for (int j = 0; j < graph.n; ++j) {
        if (graph.k_tag[j] == 0) {
            if (graph.adj_list[j].size() < config.inputK) {
                k_tag_anchor[j] = -1;
            }
            const vector<int> &neis = graph.adj_list[j];
            for (int nei:neis) {
                if (graph.k_tag[nei] == 0) {
                    counter1++;
                    remove_core_adj_copy[j].emplace_back(nei);
                    if (graph.adj_list[j].size() >= config.inputK && graph.adj_list[nei].size() >= config.inputK) {
                        remove_core_adj[j].emplace_back(nei);
                        counter2++;
                    }
                }
            }
        }
    }
    vector<bool> tmp(graph.n, false);
    for (int k = 0; k < graph.n; ++k) {
        if (k_tag_anchor[k] == 0 && remove_core_adj[k].size() + kc_nei[k] < config.inputK) {
            tmp[k] = true;
        }
    }
    for (int l = 0; l < remove_core_adj_copy.size(); ++l) {
        if (k_tag_anchor[l] == -1 && int(graph.adj_list[l].size()) >= 0) {//config.inputK - 5) {
            //先判断邻居中有无锚点，单纯的度并无意义，主要是和锚点之间的边的多少，
            bool flag = false;
            for (int nei_id:remove_core_adj_copy[l]) {
                if (tmp[nei_id]) {//重要
                    flag = true;
                    break;;
                }
            }
            //有锚点的话，添加到该点
            if (flag) {
                for (int nei_id:remove_core_adj_copy[l]) {
                    if (remove_core_adj_copy[nei_id].size() + kc_nei[nei_id] >= config.inputK) {
                        remove_core_adj[nei_id].emplace_back(l);
                        remove_core_adj[l].emplace_back(nei_id);
                    }
                }
                k_tag_anchor[l] = 0;
            }
        }
    }
    //saveGraph("raw.csv");
    //expand_input_graph();
    //到这里应该对比一下两种情况的子图大小。
    unordered_map<int, int> init_anchor_deg, add_anchor_deg;
    //anchor_group_vec.emplace_back(AnchorGroup());
    for (int k = 0; k < graph.n; ++k) {
        if (k_tag_anchor[k] == 0) {
            if (remove_core_adj[k].size() + kc_nei[k] < config.inputK) {
                whole_anchor++;
                k_tag_anchor[k] = 2;
                node_group[k].emplace(anchor_group_vec.size());
                anchor_group_vec.emplace_back(AnchorGroup(k));
                if (graph.adj_list[k].size() < config.inputK) {
                    add_anchor_deg[graph.adj_list[k].size()]++;
                }
                init_anchor_deg[graph.adj_list[k].size()]++;
            } else {
                whole_follow++;
            }
        }
    }
}

vector<int> SubGraph::delete_node(int deln) {
    const vector<int> &neighbors = remove_core_adj[deln];
    for (int nei:neighbors) {
        const vector<int> &neineighbors = remove_core_adj[nei];
        for (int j = 0; j < neineighbors.size(); ++j) {
            if (remove_core_adj[nei][j] == deln) {
                remove_core_adj[nei].erase(remove_core_adj[nei].begin() + j);
                break;
            }
        }
    }
    remove_core_adj[deln].clear();
    k_tag_anchor[deln] = -1;
    return neighbors;
}

void SubGraph::insert_node(int ins_node) {
    //assert(k_tag_anchor[ins_node] == -1);
    for (int nei:remove_core_adj_copy[ins_node]) {
        if (k_tag_anchor[nei] == -1 || k_tag_anchor[nei] == 1)continue;
        remove_core_adj[nei].emplace_back(ins_node);
        remove_core_adj[ins_node].emplace_back(nei);
    }
    k_tag_anchor[ins_node] = 0;
}


vector<int> SubGraph::continuous_del_nodes(int first_node) {
    vector<int> del_nodes{first_node};
    for (int i = 0; i < del_nodes.size(); ++i) {
        int node = del_nodes[i];
        vector<int> neighbors = remove_core_adj[node];
        delete_node(node);
        for (int nei:neighbors) {
            if (k_tag_anchor[nei] == 0 && kc_nei[nei] + remove_core_adj[nei].size() == config.inputK - 1) {
                del_nodes.emplace_back(nei);
            }
        }
    }
    vector<int> final_del_nodes = del_nodes;
    return del_nodes;
}


//逻辑是新添加的点的周围需要更新支持者，而本锚组的支持者关系到的只需要更新得分
void SubGraph::add_anchor_group(int anc, unordered_map<int, double> &ag_sup_update,
                                unordered_map<int, double> &ag_score_update) {
    //如果是支持者，那么其锚点组肯定是有值的！
    anchor_group_vec.emplace_back(AnchorGroup(anc));
    int ag_id = anchor_group_vec.size() - 1;
    detect_ag_supporters(ag_id);
    anchor_group_vec[ag_id].supporters = anchor_group_vec[ag_id].new_supporters;
    if (k_tag_anchor[anc] == 0) {
        k_tag_anchor[anc] = 2;
        whole_anchor++;
    }
    ///锚点的node_group,只有因锚点引起的度变化才会导致重新计算支持者
    if (node_group[anc].empty()) {// merge with neighbors如果是外部点，那么其邻居所在的组都需要更新
        for (int nei:remove_core_adj[anc]) {
            for (int ag_id_tmp:node_group[nei]) {
                ag_sup_update[abs(ag_id_tmp)]++;
            }
        }
    }
    for (int ag:node_group[anc]) {
        ag_sup_update[abs(ag)]++;
    }
    node_group[anc].clear();
    node_group[anc].emplace(ag_id);
    ag_sup_update[ag_id]++;
    ///支持者的node_group，节点的node_group发生变化，那么分数就该变化，deg也该变化
    for (int sup: anchor_group_vec[ag_id].supporters) {
        node_group[sup].emplace(ag_id);
        for (int ag_id_tmp:node_group[sup]) {
            ag_score_update[abs(ag_id_tmp)]++;
        }
    }
}


/*
* 1. 可能删除多个锚点组
* 2. 对其他锚点组支持者和得分的更新
* 3. 解除锚定和更新node_group的state
 * 4. 考虑锚点和支持者两个。
*/

int SubGraph::del_anchor_groups(const vector<int> &ags_need_del, unordered_map<int, double> &ag_sup_update,
                                unordered_map<int, double> &ag_score_update, bool partial) {
    unordered_map<int, int> del_nodes_whole;
    int return_value = 0;
    for (int ag_id:ags_need_del) {
        if (anchor_group_vec[ag_id].state == -1)continue;
        //先修改现在支持者和锚点的node_group
        //assert(anchor_group_vec[ag_id].state != -1);
        anchor_group_vec[ag_id].state = -1;
        vector<int> test_nodes = anchor_group_vec[ag_id].supporters;
        //一个条件吧，主要是删除的点过于多，那么就把k-1中的点设置为锚点。
        whole_anchor--;
        int anc = anchor_group_vec[ag_id].anc;
        test_nodes.emplace_back(anc);
        node_group[anc].erase(ag_id);
        //assert(node_group_deg[anc] == anchor_group_vec[ag_id].anchors.size());
        k_tag_anchor[anc] = 0;
        vector<int> del_nodes = continuous_del_nodes(anc);

        for (int deln:del_nodes) {
            del_nodes_whole[deln]++;
        }
        if (whole_anchor % 100 == 0) {
            INFO(whole_anchor, whole_follow);
        }
        whole_follow -= del_nodes.size() - 1;
        while (whole_follow + whole_anchor < fol_thre) {
            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - result.startTime).count();
            INFO(duration / TIMES_PER_SEC);
            result.fol_anc_time[int(fol_thre_number * 10 + 0.1)] = make_pair(whole_anchor + 1,
                                                                             duration / TIMES_PER_SEC);
            fol_thre_number -= 0.1;
            fol_thre = result.deg_more_k * fol_thre_number;
        }
        //anchor_group_vec[ag_id].print();
        return_value += del_nodes.size() - 1;
        for (int sup: anchor_group_vec[ag_id].supporters) {
            node_group[sup].erase(ag_id);
        }
    }
    //锚点的邻居以及锚点的邻居的邻居都需要更新支持者。
    unordered_map<int, int> all_nodes;
    for (const auto &item:del_nodes_whole) {
        all_nodes[item.first] = 1;
        const vector<int> &neis = remove_core_adj_copy[item.first];
        for (int nei:neis) {
            all_nodes[nei] = 1;
            const vector<int> &neineis = remove_core_adj[nei];
            for (int neinei:neineis) {
                all_nodes[neinei] = 1;
            }
        }
    }
    for (const auto &item: all_nodes) {
        if (k_tag_anchor[item.first] == -1) {
            node_group[item.first].clear();
        } else {
            for (int ag_tmp:node_group[item.first]) {
                if (anchor_group_vec[abs(ag_tmp)].state != -1)
                    ag_sup_update[abs(ag_tmp)] = 1;
            }
        }
    }
    return return_value;
}
// update anchor group 分为计算支持者和得分
/*
 * 1. 对一组锚点集可能删除可能添加更新支持者和得分
 * 2. 更新node_group 和索引。
 * 3. 如果是单个锚点，且是添加的，那么需要检测是否可以再加一个邻居的。！！！
 * 4. 既有删除的也有增加的
 * 5. 需要更新支持者，删除的锚点组，添加的锚点在图中则需要其本身，否则包含其邻居
 * anchors	supporters	new_supporters	score	state	whole_anchor	k_tag_anchor	set_index2	node_group（锚点和支持者）	group_deg	adj_list	anchor_group_vec
 */

void SubGraph::campare_new_old_sup(int ag_id) {
    for (int sup: anchor_group_vec[ag_id].supporters) {
        supporter_record[sup] = function_counter;
        supporter_order[sup] = ag_id;
    }
    for (int sup: anchor_group_vec[ag_id].new_supporters) {
        new_supporter_record[sup] = function_counter;
        new_supporter_order[sup] = ag_id;
    }
}

unordered_map<int, int>
SubGraph::update_supporters(unordered_map<int, double> &ag_sup_update, unordered_map<int, double> &ag_score_update) {
    vector<int> need_divide_groups, anchors_for_new_group;
    unordered_map<int, int> ag_real_change_sup, nodes_ag_changed;
    //先计算锚组的支持者，然后先删除node_groups更新supporters
    for (const auto &item:ag_sup_update) {
        int ag_id = item.first;
        AnchorGroup &anchorGroup = anchor_group_vec[ag_id];
        //assert(anchorGroup.state != -1);
        detect_ag_supporters(ag_id);
        campare_new_old_sup(ag_id);
        int counter = 0;
        for (int sup : anchor_group_vec[ag_id].new_supporters) {
            if (supporter_record[sup] == new_supporter_record[sup]) {
                counter++;
            } else {
                ag_real_change_sup[ag_id] = 1;
                nodes_ag_changed[sup] = 1;
                node_group[sup].emplace(new_supporter_order[sup]);
            }
        }
        //旧的supporters可能成为了锚点
        for (int sup:anchor_group_vec[ag_id].supporters) {
            if (abs(supporter_record[sup]) != abs(new_supporter_record[sup])) {
                ag_real_change_sup[ag_id] = 1;
                nodes_ag_changed[sup] = 1;
                if (k_tag_anchor[sup] != 2) {
                    node_group[sup].erase(supporter_order[sup]);
                }
            }
        }
        anchorGroup.supporters = anchorGroup.new_supporters;
    }
    sort(need_divide_groups.begin(), need_divide_groups.end());
    for (auto item:nodes_ag_changed) {
        for (int ag:node_group[item.first]) {
            ag_score_update[abs(ag)] = 1;
        }
    }
    return ag_real_change_sup;
}

vector<int>
SubGraph::update_anchor_group(unordered_map<int, double> &ag_sup_update, unordered_map<int, double> &ag_score_update) {
    unordered_map<int, int> ag_real_change_sup = update_supporters(ag_sup_update, ag_score_update);
    vector<int> ag_sup_changes;
    for (const auto &item: ag_real_change_sup) {
        ag_sup_changes.emplace_back(item.first);
    }
    //对于得分需要重新计算的，重新计算得分放入索引中。
    vector<int> ags_index_update;
    int old_ag_score_size = ag_score_update.size();
    ag_score_update.insert(ag_sup_update.begin(), ag_sup_update.end());
    unordered_map<int, bool> all_sup;
    for (auto &item:ag_score_update) {
        int ag = item.first;
        if (anchor_group_vec[ag].state >= 0) {

            item.second = anchor_group_vec[ag].score;

            double old_score = anchor_group_vec[ag].score;
            compute_score(ag);
            for (int sup:anchor_group_vec[ag].supporters) {
                all_sup[sup] = true;
            }
        }
    }
    for (const auto &item:ag_score_update) {
        int ag = item.first;
        if (anchor_group_vec[ag].state < 0) continue;
        if (anchor_group_vec[ag].state == 0 || compare_double(
                anchor_group_vec[ag].score,
                item.second) != 0) {
            index_push(ag);
            merge_index1_push(ag);//////////这个是新加的
            ags_index_update.emplace_back(ag);
        }
    }
    //return ags_index_update;
    return ag_sup_changes;
}


void SubGraph::detect_ag_supporters(int ag_id) {
    Timer timer(SUPPORTER_COMPUTE);
    //unordered_map<int ,bool> del_nodes;
    anchor_group_vec[ag_id].stability = 0;
    function_counter++;
    queue<int> active_nodes;
    int anc = anchor_group_vec[ag_id].anc;
    //先修改状态，然后更改值，然后查询值前，先比较状态
    del_nodes_state[anc] = function_counter;
    delete_nodes[anc] = true;
    active_nodes.push(anc);
    anchor_group_vec[ag_id].new_supporters.clear();
    while (!active_nodes.empty()) {
        int node = active_nodes.front();
        active_nodes.pop();
        const vector<int> &neis_node = remove_core_adj[node];
        for (int nei:neis_node) {
            //如果在当前的锚组中已经被删除了，那么就算了，否则的话，
            if (del_nodes_state[nei] == function_counter) {
                if (delete_nodes[nei]) { continue; }
            } else {
                del_nodes_state[nei] = function_counter;
                delete_nodes[nei] = false;
            }
            anchor_group_vec[ag_id].stability += 1;
            if (nei_score_state[nei] == function_counter) {
                --neighbor_score[nei];
            } else {
                nei_score_state[nei] = function_counter;
                neighbor_score[nei] = -1;
            }
            int real_deg = remove_core_adj[nei].size() + kc_nei[nei] + neighbor_score[nei];
            //判断node_kcoreNei是否大于0，来避免删除锚点，上面判断是否是已删除的点，然后结合现在的得分计算
            if (k_tag_anchor[nei] == 0 && real_deg < config.inputK) {
                delete_nodes[nei] = true;
                active_nodes.push(nei);
                anchor_group_vec[ag_id].new_supporters.emplace_back(nei);
                assert(k_tag_anchor[nei] == 0);
                if (anchor_group_vec[ag_id].new_supporters.size() > max_supporters) {
                    return;
                }
            } else if (k_tag_anchor[nei] == 0 && real_deg == config.inputK &&
                       sub_nei_score_state[nei] != function_counter) {
                sub_nei_score_state[nei] = function_counter;
            }
        }
    }
}

bool SubGraph::detect_ag_supporters_without_threold(int ag_id) {
    Timer timer(SUPPORTER_COMPUTE);
    //unordered_map<int ,bool> del_nodes;
    anchor_group_vec[ag_id].stability = 0;
    function_counter++;
    queue<int> active_nodes;
    int anc = anchor_group_vec[ag_id].anc;
    //先修改状态，然后更改值，然后查询值前，先比较状态
    del_nodes_state[anc] = function_counter;
    delete_nodes[anc] = true;
    active_nodes.push(anc);
    anchor_group_vec[ag_id].new_supporters.clear();
    while (!active_nodes.empty()) {
        int node = active_nodes.front();
        active_nodes.pop();
        const vector<int> &neis_node = remove_core_adj[node];
        for (int nei:neis_node) {
            //如果在当前的锚组中已经被删除了，那么就算了，否则的话，
            if (del_nodes_state[nei] == function_counter) {
                if (delete_nodes[nei]) { continue; }
            } else {
                del_nodes_state[nei] = function_counter;
                delete_nodes[nei] = false;
            }
            anchor_group_vec[ag_id].stability += 1;
            if (nei_score_state[nei] == function_counter) {
                --neighbor_score[nei];
            } else {
                nei_score_state[nei] = function_counter;
                neighbor_score[nei] = -1;
            }
            int real_deg = remove_core_adj[nei].size() + kc_nei[nei] + neighbor_score[nei];
            //判断node_kcoreNei是否大于0，来避免删除锚点，上面判断是否是已删除的点，然后结合现在的得分计算
            if (k_tag_anchor[nei] == 0 && real_deg < config.inputK) {
                delete_nodes[nei] = true;
                active_nodes.push(nei);
                anchor_group_vec[ag_id].new_supporters.emplace_back(nei);
                if (anchor_group_vec[ag_id].new_supporters.size() > factor_max_supporters * max_supporters) {
                    more1000_group_count[ag_id]++;
                    return true;
                }

            } else if (k_tag_anchor[nei] == 0 && real_deg == config.inputK &&
                       sub_nei_score_state[nei] != function_counter) {
                sub_nei_score_state[nei] = function_counter;
            }
        }
    }
    return false;
}

void SubGraph::compute_score(int ag_id) {
    if (more1000_group_count.find(ag_id) != more1000_group_count.end()) {
        return;
    }
    AnchorGroup &ag = anchor_group_vec[ag_id];
    ag.score = 0;
    //ag.score=anchor_group_vec[ag_id].supporters.size()*1.0/anchor_group_vec[ag_id].anchors.size();

    for (int sup:ag.supporters) {
        assert(node_group[sup].size() != 0);
        if (node_group[sup].size() == 1) {
            ag.score += 1.0;
        } else {
            ag.score += 1.0 / node_group[sup].size();
        }
    }
    ag.old_score = ag.score;
    //ag.score = anchor_group_vec[ag_id].supporters.size();
    if (ag.score < -0.1) {
        cout << ag_id << endl;
    } else if (ag.score > 10000) {
        cout << ag_id << endl;
    }
    return;
}


void SubGraph::index_push(int ag_id) {
    anchor_group_vec[ag_id].state = 1;
    set_index2.push(ans_index(ag_id, 1, anchor_group_vec[ag_id].score,
                              anchor_group_vec[ag_id].supporters.size(), anchor_group_vec[ag_id].stability));
}

void SubGraph::merge_index1_push(int ag_id) {
    anchor_group_vec[ag_id].tag = tag;
    merge_index1.push(
            ans_index_big_first(ag_id, 1, anchor_group_vec[ag_id].score,
                                anchor_group_vec[ag_id].supporters.size(), anchor_group_vec[ag_id].stability));
}


candidate_anc
SubGraph::detect_merged_groups2(vector<int> &need_merged_ags, int direct) {
    function_counter++;
    vector<int> new_common_followers_all, merged_set_true;
    int merged_followers = 0, base_ag = need_merged_ags[0];
    long long neighbor_func_counter = function_counter;
    double del_anc_count = 0, del_anc_count2 = 0, del_anc_count3 = 0;
    for (int m = 0; m < need_merged_ags.size(); ++m) {
        bool flag = true;
        int anc = anchor_group_vec[need_merged_ags[m]].anc;
        function_counter++;
        //AnchorGroup &anchorGroup = anchor_group_vec[ag_id];

        double thre = min(ceil(anchor_group_vec[need_merged_ags[m]].supporters.size() * thre_detect),
                          (whole_follow * (need_merged_ags.size() * 1.0 - m)) / whole_anchor);//这里是因为最大的限制。
        //queue<int> active_nodes;
        vector<int> all_nodes{anc};
        del_nodes_state[anc] = neighbor_func_counter;
        delete_nodes[anc] = true;

        int need_update_nodes_num = 0, del_fol = 0;
        for (int i = 0; i < all_nodes.size() && flag; ++i) {
            int node = all_nodes[i];
            const vector<int> &neis_node = remove_core_adj[node];
            for (int nei:neis_node) {
                if (del_nodes_state[nei] == neighbor_func_counter) {
                    if (delete_nodes[nei]) { continue; }
                } else {
                    delete_nodes[nei] = false;
                    del_nodes_state[nei] = neighbor_func_counter;
                }
                if (nei_score_state[nei] != neighbor_func_counter) {
                    neighbor_score[nei] = 0;
                    nei_score_state[nei] = neighbor_func_counter;
                }
                if (sub_nei_score_state[nei] == function_counter) {
                    --sub_neighbor_score[nei];
                } else {
                    need_update_nodes[need_update_nodes_num++] = nei;
                    sub_nei_score_state[nei] = function_counter;
                    sub_neighbor_score[nei] = -1;
                }
                //判断node_kcoreNei是否大于0，来避免删除锚点，上面判断是否是已删除的点，然后结合现在的得分计算
                if (k_tag_anchor[nei] == 0 &&
                    remove_core_adj[nei].size() + kc_nei[nei] + neighbor_score[nei] + sub_neighbor_score[nei] <
                    config.inputK) {
                    delete_nodes[nei] = true;
                    all_nodes.emplace_back(nei);
                    del_fol++;
                    if (all_nodes.size() > thre + 1) {
                        flag = false;
                        break;
                    }
                }
            }
        }
        if (del_fol == 1) {
            del_anc_count2++;
        }
        del_anc_count3 += pow(2, del_fol * -1);
        if (merged_set_true.empty() && base_ag != need_merged_ags[m]) {//如果基础锚组也被删除，则不合并
            return candidate_anc(-1, 0, 0, 0);
        }
        if (flag) {
            merged_set_true.emplace_back(need_merged_ags[m]);
            for (int i = 1; i < all_nodes.size(); ++i) {
                new_common_followers_all.emplace_back(all_nodes[i]);
            }
            for (int j = 0; j < need_update_nodes_num; ++j) {
                int node = need_update_nodes[j];
                neighbor_score[node] += sub_neighbor_score[node];
            }
            // 尝试计算所有点的边？
            if (del_fol > 0)
                del_anc_count += pow(2, del_fol * -1);
        } else {
            for (int i = 1; i < all_nodes.size(); ++i) {
                delete_nodes[all_nodes[i]] = false;
            }
            delete_nodes[anc] = false;
        }
    }
    need_merged_ags = merged_set_true;
    if (need_merged_ags.size() < 2) {
        return candidate_anc(-1, 0, 0, 0);
    }
    if (direct == 0) {
        return candidate_anc(-1, need_merged_ags.size() - 1, new_common_followers_all.size() + 1, del_anc_count,
                             del_anc_count2,
                             del_anc_count3);
    } else if (direct == 2) {
        return candidate_anc(-1, need_merged_ags.size() - 1, new_common_followers_all.size(), del_anc_count,
                             del_anc_count2,
                             del_anc_count3);
    }
    return candidate_anc(0, 0, 0, 0);
}


candidate_anc
SubGraph::merge_with_supporters2(const unordered_map<int, double> &node_score, const unordered_map<int, int> &node_deg,
                                 vector<int> &need_merge_ags) {
    //对支持者排序并按照支持者的组数分类。
    vector<pair<int, double>> node_score_vec;
    node_score_vec.reserve(node_score.size());
    for (auto v:node_score) {
        node_score_vec.emplace_back(v);
    }
    sort(node_score_vec.begin(), node_score_vec.end(), cmp_pair);
    unordered_map<int, vector<int>> fol_pri_que, remove_out_nei;
    for (auto v:node_score_vec) {
        if (node_deg.at(v.first) == 1) continue;
        fol_pri_que[node_deg.at(v.first)].emplace_back(v.first);
    }
    //从前3个里面找出最合适的点？此处先找最合适的?
    int max_node = -1;
    candidate_anc return_value(0, 0, 0, 0);
    vector<int> final_merged_set;
    for (auto &item:fol_pri_que) {
        int test_nodes_num = config.supporter_num;
        for (int k = 0; k < test_nodes_num && k < item.second.size(); ++k) {
            //只考虑base group的supporter
            int node_id = item.second[k];
            vector<int> merged_set;
            for (int ag_id:need_merge_ags) {
                if (node_group[node_id].find(ag_id) != node_group[node_id].end()) {
                    merged_set.emplace_back(ag_id);
                }
            }
            //锚定该点后检测锚组删除情况
            k_tag_anchor[node_id] = 2;
            candidate_anc ca = detect_merged_groups2(merged_set);
            ca.id = node_id;
            k_tag_anchor[node_id] = 0;
            if (merged_set.size() < 2)continue;
            if (cmp_cand_anc(return_value, ca)) {
                return_value = ca;
                final_merged_set = merged_set;
                //break;
            }
        }
    }
    need_merge_ags = final_merged_set;
    return return_value;
}


candidate_anc
SubGraph::merge_with_neighbors2(const unordered_map<int, double> &node_score, const unordered_map<int, int> &node_deg,
                                vector<int> &need_merge_ags) {
    //对邻居按照分数排序并按照邻居相关的组数分类。
    vector<pair<int, double>> node_score_vec;
    node_score_vec.reserve(node_score.size());
    for (auto v:node_score) {
        node_score_vec.emplace_back(v);
    }
    sort(node_score_vec.begin(), node_score_vec.end(), cmp_pair);
    unordered_map<int, vector<int>> fol_pri_que, remove_out_nei;
    for (auto v:node_score_vec) {
        if (node_deg.at(v.first) == 1) continue;
        fol_pri_que[node_deg.at(v.first)].emplace_back(v.first);
    }

    int max_node = -1;
    vector<int> final_merged_set;
    candidate_anc return_value(0, 0, 0, 0);
    for (auto &item:fol_pri_que) {
        int test_nodes_num = config.neighbor_num;
        for (int k = 0; k < test_nodes_num && k < item.second.size(); ++k) {
            int node_id = item.second[k];
            vector<int> merged_set;
            //检测是否和之前的分数相同，不是要找排名前三的节点，而是要找排名前三的分数的所有节点。
            //锚定并插入该点，然后检测合并的锚点集
            insert_node(node_id);
            for (int ag_id:need_merge_ags) {
                for (int nei:remove_core_adj[node_id]) {
                    if (node_group[nei].find(ag_id) != node_group[nei].end()) {
                        merged_set.emplace_back(ag_id);
                        break;
                    }
                }
            }
            //锚定该点后检测锚组删除情况
            k_tag_anchor[node_id] = 2;
            candidate_anc ca = detect_merged_groups2(merged_set, 0);
            ca.id = node_id;
            k_tag_anchor[node_id] = 0;
            delete_node(node_id);
            if (merged_set.size() < 2)continue;
            if (cmp_cand_anc(return_value, ca)) {
                return_value = ca;
                final_merged_set = merged_set;
                break;
            }
        }
    }
    need_merge_ags = final_merged_set;
    return return_value;
}

vector<int> SubGraph::merge_anchor_groups2(vector<int> ags, bool next_phase) {
    Timer timer1(MERGE_STEP);
    int base_ag = ags[0];
    priority_queue<ans_index> tmp;//这里是为了调整锚组的顺序
    for (int k = 1; k < ags.size(); ++k) {
        tmp.push(ans_index(ags[k], 1,
                           anchor_group_vec[ags[k]].score, anchor_group_vec[ags[k]].supporters.size(),
                           anchor_group_vec[ags[k]].stability));
    }
    for (int l = 1; l < ags.size(); ++l) {
        ags[l] = tmp.top().id;
        tmp.pop();
    }
    vector<int> ags_directly = ags, ags_supporters = ags, ags_outnode = ags;
    unordered_map<int, double> sup_seq, nei_seq;
    unordered_map<int, int> sup_ag_count, nei_count, node_ag, node_count;
    for (int ag:ags) {
        int sup_size = anchor_group_vec[ag].supporters.size();
        for (int j = 0; j < sup_size; ++j) {
            int sup = anchor_group_vec[ag].supporters[j];
            node_count[sup] = 1;
            if (ag == base_ag || sup_seq.find(sup) != sup_seq.end()) {//如果是基础锚组或者是基础锚组的支持者
                sup_seq[sup] += (double(j + 1)) / sup_size;//这个应该是对于比较长的更友好吧
                ++sup_ag_count[sup];
            }
            for (int nei:remove_core_adj_copy[sup]) {
                if (k_tag_anchor[nei] != -1) continue;
                if (ag != base_ag && nei_seq.find(nei) == nei_seq.end()) continue; //如果是基础锚组或者是基础锚组的邻居
                if (node_ag.find(nei) != node_ag.end() && node_ag[nei] == ag) continue; //如果是多个支持者的邻居，则只考虑最高的
                node_ag[nei] = ag; //只计算该点连接的最高的ag
                nei_seq[nei] += double(j + 1) / sup_size;
                ++nei_count[nei];
            }
        }
    }
    //merge with supporters
    candidate_anc max_in(0, 0, 0, 0);
    if (config.mode > 1) {
        max_in = merge_with_supporters2(sup_seq, sup_ag_count, ags_supporters);
    }
    // merge with out neighbors
    candidate_anc max_out(0, 0, 0, 0);
    if (config.mode > 2) {
        max_out = merge_with_neighbors2(nei_seq, nei_count, ags_outnode);
    }
    bool cmp_result = cmp_cand_anc(max_out, max_in);
    unordered_map<int, double> ag_sup_update, ag_score_update;//有的被删除的，还是需要调整。
    if (cmp_result == 0 && ags_outnode.size() > 1 && max_out.del_ancs > 0 && config.mode > 2) {
        double single_ratio = max_out.del_fol * 1.0 / max_out.del_ancs;
        //锚定目标点并插入目标点！
        if (single_ratio >= thre_whole) {
            //把合并的这些锚组设置标志位
            //anchor_group_vec[base_ag].tag = !tag;
            for (int ag_id:ags_outnode) {
                anchor_group_vec[ag_id].tag = !tag;
            }
        } else {
            insert_node(max_out.id);
            k_tag_anchor[max_out.id] = 2;
            whole_anchor++;
            int del_num = del_anchor_groups(ags_outnode, ag_sup_update, ag_score_update);
            add_anchor_group(max_out.id, ag_score_update, ag_score_update);
            //int del_num = delete_add_ags(max_out.id, ags_outnode, ag_sup_update, ag_score_update);
            //assert(del_num==max_out.del_fol);
            merge_neighbor_num++;
        }
    } else if (cmp_result == 1 && max_in.del_ancs > 0 && ags_supporters.size() > 1 && config.mode > 1) {
        double whole_ratio = whole_follow * 1.0 / whole_anchor, single_ratio =
                max_in.del_fol * 1.0 / max_in.del_ancs, anc = max_in.del_ancs, fol = max_in.del_fol;
        if (single_ratio >= thre_whole) {
            //把合并的锚点的设置标志位
            //anchor_group_vec[base_ag].tag = !tag;
            for (int ag_id:ags_supporters) {
                anchor_group_vec[ag_id].tag = !tag;
            }
        } else {
            k_tag_anchor[max_in.id] = 2;
            whole_follow--;
            whole_anchor++;
            int del_num = del_anchor_groups(ags_supporters, ag_sup_update, ag_score_update);
            add_anchor_group(max_in.id, ag_score_update, ag_score_update);
            //int del_num = delete_add_ags(new_ancs, ags_supporters, ag_sup_update, ag_score_update);
            //assert(del_num==max_in.del_fol-1);
            merge_support_num++;
        }
    }
    merge_whole++;
    new_support_count[ag_sup_update.size()]++;
    return update_anchor_group(ag_sup_update, ag_score_update);
}

void SubGraph::adjust_parameters() {
    unordered_map<int, double> ag_score_update;
    //先调整参数
    thre_whole *= 2;
    //thre_whole += 1;
    if (max_supporters < 80) {
        max_supporters *= 2;
        unordered_map<int, int> nodes_ag_changed;
        for (int i = 0; i < anchor_group_vec.size(); ++i) {
            if (anchor_group_vec[i].state != -1) {
                merge_index1_push(i);
                index_push(i);
            }
        }
    } else {
        factor_max_supporters *= 2;
    }
    config.supporter_num = ceil((max_supporters + 1) * thre_detect);
    config.neighbor_num = ceil((max_supporters + 1) * thre_detect);
    tag = !tag;
    //把所有的锚组从头计算
}

/*
void SubGraph::merge_big_groups() {
    int anc_num = 0, group_num = 0, whole_group = 0, whole_anc_num = 0;
    unordered_map<int, int> big_fenbu, whole_fenbu;
    for (auto item:more1000_group_count) {
        int ag_id = item.first;
        if (anchor_group_vec[ag_id].state != -1) {
            anc_num += anchor_group_vec[ag_id].anchors.size();
            group_num++;
            big_fenbu[anchor_group_vec[ag_id].anchors.size()]++;
        }
    }
    for (int i = 0; i < anchor_group_vec.size(); ++i) {
        if (anchor_group_vec[i].state != -1) {
            bool flag = more1000_group_count.find(i) == more1000_group_count.end();
            if (flag) {
                //如果他没有的话
                whole_fenbu[anchor_group_vec[i].anchors.size()]++;
                whole_anc_num += anchor_group_vec[i].anchors.size();
                whole_group++;
            } else {
                cout << more1000_group_count.at(i);
            }
        }
    }
    INFO(whole_anchor, whole_follow);
    INFO(anc_num, group_num, whole_anc_num, whole_group);
    cout << "大锚组分布" << endl << "anc_num:\tcount" << endl;
    for (int i = 0; i < whole_anchor; ++i) {
        if (big_fenbu.find(i) != big_fenbu.end()) {
            cout << i << "\t" << big_fenbu[i] << endl;
        }
    }
    cout << "全部分布" << endl << "anc_num:\tcount" << endl;
    for (int i = 0; i < whole_anchor; ++i) {
        if (whole_fenbu.find(i) != whole_fenbu.end()) {
            cout << i << "\t" << whole_fenbu[i] << endl;
        }
    }
    for (int i = 0; i < 10; ++i) {
        if (result.fol_anc_time.find(i) != result.fol_anc_time.end()) {
            INFO(i, result.fol_anc_time[i].first, result.fol_anc_time[i].second);
        }
    }
}
*/
void SubGraph::main_algorithm() {
    //检测支持者，并为每个supporter更新node_group，node_group_deg
    for (int j = 0; j < anchor_group_vec.size(); ++j) {
        detect_ag_supporters(j);
        anchor_group_vec[j].supporters = anchor_group_vec[j].new_supporters;
        for (int sup:anchor_group_vec[j].supporters) {
            node_group[sup].emplace(j);
        }
    }
    for (int k = 0; k < anchor_group_vec.size(); ++k) {
        compute_score(k);
        index_push(k);
        if (anchor_group_vec[k].supporters.empty()) continue;
        merge_index1_push(k);
    }
    //这个是为了保存已经删除过的锚组的顺序，为最后一步的添加点做准备。
    bool print = true, next_phase = false;
    //先删除一遍fol为0的
    while (whole_anchor > config.inputB) {
        //合并 得分由低到高
        //INFO(whole_anchor);
        while (!merge_index1.empty() && config.mode > 0 && whole_anchor > config.inputB) {
            ans_index_big_first ansIndex = merge_index1.top();
            merge_index1.pop();
            //索引中可能有一个锚组的多个记录，找出和现在相同的，再进行计算
            if (anchor_group_vec[ansIndex.id].state == 1 && anchor_group_vec[ansIndex.id].tag == tag &&
                compare_double(anchor_group_vec[ansIndex.id].old_score, ansIndex.score) == 0) {
                //构造AGS，然后场尝试合并
                int ag_id = ansIndex.id;
                int tmp = 0;
                vector<int> ags{ag_id};
                unordered_map<int, int> group_count;
                for (int sup:anchor_group_vec[ag_id].supporters) {
                    for (int group_id:node_group[sup]) {
                        if (abs(group_id) != ag_id) {
                            ++group_count[abs(group_id)];
                        }
                    }
                }
                double thre = thre_merge_detect;
                for (auto item:group_count) {
                    if (item.second > anchor_group_vec[ag_id].supporters.size() * thre_merge_detect) {
                        //if (anchor_group_vec[item.first].anchors.size() > 1)continue;
                        ags.emplace_back(item.first);
                    }
                }
                if (ags.size() < 2)continue;
                int old_size = anchor_group_vec.size();
                vector<int> new_add_groups = merge_anchor_groups2(ags, next_phase);
                for (int id:new_add_groups) {
                    merge_index1_push(id);
                }
            }
        }
        //删除 删除并更新，如果有新加的，那么重新检测合并。
        while (whole_anchor > config.inputB) {
            if (set_index2.empty()) {
                //merge_big_groups();
                adjust_parameters();
                break;
            }
            unordered_map<int, double> ag_sup_update, ag_score_update;
            vector<int> new_add_groups;
            ans_index index_tmp = set_index2.top();
            set_index2.pop();
            //INFO(index_tmp.id);
            if (anchor_group_vec[index_tmp.id].state != 1 ||
                compare_double(anchor_group_vec[index_tmp.id].score,
                               index_tmp.score) != 0)
                continue;
/*
            if (more1000_group_count.find(index_tmp.id) != more1000_group_count.end()) {
                continue;
            } else if (detect_ag_supporters_without_threold(index_tmp.id)) {
                continue;
            }
*/
            int del_fols = del_anchor_groups(vector<int>{index_tmp.id}, ag_sup_update, ag_score_update, false);
            //INFO(whole_anchor, index_tmp.anc_num, index_tmp.fol_num);
            del_group_count[del_fols]++;
            new_support_count[ag_sup_update.size()]++;
            new_add_groups = update_anchor_group(ag_sup_update, ag_score_update);

            ag_score_update.clear();
            direct_del_fol.push(index_tmp.fol_num * 1.0 / index_tmp.anc_num);
            direct_del_fol_sum += index_tmp.fol_num * 1.0 / index_tmp.anc_num;
            if (direct_del_fol.size() > del_fol_queue_size) {
                direct_del_fol_sum -= direct_del_fol.front();
                direct_del_fol.pop();
            }
            if (direct_del_fol_sum / del_fol_queue_size > thre_whole) {
                //merge_big_groups();
                adjust_parameters();
                break;
            } else if (!new_add_groups.empty()) {
                for (int id: new_add_groups) {
                    unordered_map<int, int> group_count;
                    for (int sup:anchor_group_vec[id].supporters) {
                        for (int group_id:node_group[sup]) {
                            if (abs(group_id) != id) {
                                ++group_count[abs(group_id)];
                            }
                        }
                    }
                    for (auto item:group_count) {
                        if (item.second > anchor_group_vec[item.first].supporters.size() * 0.5) {
                            merge_index1_push(item.first);
                        }
                    }
                    merge_index1_push(id);
                }
                break;
            }

        }
    }
    result.follower_num=whole_anchor+whole_follow;

    int anc_num=0,follow_num=0;
    for (int l = 0; l < k_tag_anchor.size(); ++l) {
        if (k_tag_anchor[l] == 2) {
            anc_num++;
        }else if (k_tag_anchor[l]== 0){
            follow_num++;
        }
    }
    INFO(anc_num,whole_anchor,follow_num,whole_follow);
}


