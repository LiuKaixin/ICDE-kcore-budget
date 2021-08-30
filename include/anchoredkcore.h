//
// Created by 刘凯鑫 on 2020/1/24.
//

#ifndef GETUNDIRECTGRAPH_ANCHOREDKCORE_H
#define GETUNDIRECTGRAPH_ANCHOREDKCORE_H

#include "graph.h"
#include "anchorgroup.h"

class SubGraph {
public:
    double fol_thre = 0;
    double fol_thre_number = 0.9;

    int thre_whole = 1, del_fol_queue_size = 50;
    bool tag = true;
    double direct_del_fol_sum = 0;
    queue<double> direct_del_fol;
    int max_supporters = 5, factor_max_supporters = 10;

    double thre_detect = 0.2, thre_merge_detect = 0.5;

    int whole_anchor = 0, whole_follow = 0;
    long long function_counter = 0;
    vector<int> kc_nei, k_tag_anchor; // 0 不是kcore ，1是kcore，2是锚点

    priority_queue<ans_index> set_index2;
    priority_queue<ans_index_big_first> merge_index1;//对于比较删除的比较多的，先候补着
    vector<set<int>> node_group;

    //unordered_map<int, set<int>> node_group;
    vector<vector<int>> remove_core_adj, remove_core_adj_copy;
    //unordered_map<int, vector<int>> adj_list,adj_list_copy; //adj_list_copy作为副本提供从外部引入锚点的可能。
    vector<AnchorGroup> anchor_group_vec;
    unordered_map<int, int> node_extra_deg, node_extra_deg_backup;//在有多组节点时，计算出每个节点组实际贡献的度

    unordered_map<int, int> del_group_count, new_support_count, more1000_group_count;
    int merge_support_num = 0, merge_neighbor_num = 0, merge_whole = 0;

    vector<long long> nei_score_state, sub_nei_score_state, del_nodes_state;
    vector<int> neighbor_score, sub_neighbor_score, need_update_nodes;
    vector<bool> delete_nodes;

    vector<int> supporter_record;//用来比较新旧的支持者
    vector<int> new_supporter_record;//用来比较新旧的支持者
    vector<int> supporter_order;//用来比较新旧的支持者
    vector<int> new_supporter_order;//用来比较新旧的支持者

    void campare_new_old_sup(int ag_id);

    //先初始化adj_list和锚点集
    SubGraph(Graph &graph, bool merge = true, bool divide = true);

    vector<int> delete_node(int deln);

    void insert_node(int ins_node);

    bool static cmpIndex(const pair<int, double> &index1, const pair<int, double> &index2) {
        return index1.second < index2.second;
    }

    void index_push(int ag_id);

    void merge_index1_push(int ag_id);

    void compute_score(int ag_id);

    void detect_ag_supporters(int ag_id);

    bool detect_ag_supporters_without_threold(int ag_id);


    candidate_anc detect_merged_groups2(vector<int> &need_merged_ags, int direct = 0);

    candidate_anc
    merge_with_supporters2(const unordered_map<int, double> &node_score, const unordered_map<int, int> &node_deg,
                           vector<int> &need_merge_ags);


    candidate_anc
    merge_with_neighbors2(const unordered_map<int, double> &node_score, const unordered_map<int, int> &node_deg,
                          vector<int> &need_merge_ags);



    void adjust_parameters();

    void merge_big_groups();


    vector<int> merge_anchor_groups2(vector<int> ags, bool next_phase);

    void main_algorithm();


    vector<int> continuous_del_nodes(int first_node);

    // add anchor group
    /*
     * 1. 需要事先对新锚点组的锚点进行锚定！！！否则添加不成功的。
     * 2. 如果是单个锚点那么需要检测是否可以再加一个邻居的。
     * 3. 返回被影响的锚点集，然后对其进行更新。
     * 4. 如果引入新点 handle insert node
     */
    void add_anchor_group(int anc, unordered_map<int, double> &ag_sup_update,
                          unordered_map<int, double> &ag_score_update);

    // delete anchor group
    /*
     * 1. 可能删除多个锚点组
     * 2. 对其他锚点组支持者和得分的更新
     * 3. 解除锚定和更新node_group的state
     * 4. 为保证其他的索引不便，删除主要是指从索引中删除，即set_index2 node_group adj_list。
     */
    int del_anchor_groups(const vector<int> &ags_need_del, unordered_map<int, double> &ag_sup_update,
                          unordered_map<int, double> &ag_score_update, bool partial = true);


    unordered_map<int, int>
    update_supporters(unordered_map<int, double> &ag_sup_update, unordered_map<int, double> &ag_score_update);

    vector<int>
    update_anchor_group(unordered_map<int, double> &ag_sup_update, unordered_map<int, double> &ag_score_update);
    // 特殊的操作 divide merge
    /*
     * 1. 删除锚点组时，如果小于索引当前最小的则直接删除，否则停止划分和合并。先删除
     * 2. 需要事先对新锚点组的锚点进行锚定！！！否则添加不成功的。
     */

    struct node_score {
        int id;
        vector<int> neis;

        node_score(const int &node_id, const vector<int> &node_neis) {
            id = node_id;
            neis = node_neis;
        }

        bool operator<(const node_score &node_score2) const {
            for (int j = config.inputK - 1; j >= 0; --j) {
                if (neis[j] > node_score2.neis[j]) {
                    return false;
                } else if (neis[j] < node_score2.neis[j]) {
                    return true;
                }
            }
            return false;
        }
    };

    //输出锚点，以及不是kcore中的点有多少kcore中的邻居，输出remove_kcore。
    void saveGraph(string outputFile) {
        FILE *fout;
        outputFile = to_string(config.inputK) + "_" + to_string(config.inputB) + "_" + outputFile;
        if ((fout = fopen(outputFile.c_str(), "w+")) != NULL) {
            INFO("Save graph to file ", outputFile);
            fprintf(fout, "%s,%s\n", "source", "target");
            for (int j = 0; j < remove_core_adj.size(); ++j) {
                for (int k = 0; k < remove_core_adj[j].size(); ++k) {
                    fprintf(fout, "%d,%d\n", j, remove_core_adj[j][k]);
                }
            }
            fclose(fout);
            outputFile = "node_" + outputFile;
            if ((fout = fopen(outputFile.c_str(), "w+")) != NULL) {
                fprintf(fout, "%s,%s,%s\n", "Id", "deg", "k_tag");
                INFO("Save node to file ", outputFile);
                for (int j = 0; j < remove_core_adj.size(); ++j) {
                    if (remove_core_adj[j].size() > 0) {
                        fprintf(fout, "%d,%d,%d\n", j, kc_nei[j], k_tag_anchor[j]);
                    }
                }
                fclose(fout);
            }
            //构造新的图
            outputFile = "expand_" + outputFile;
            vector<vector<int>> remove_core_adj_candi = remove_core_adj;
            for (int l = 0; l < remove_core_adj_copy.size(); ++l) {
                if (k_tag_anchor[l] == -1) {
                    //先判断邻居中有无锚点
                    bool flag = false;
                    for (int nei_id:remove_core_adj_copy[l]) {
                        if (k_tag_anchor[nei_id] == 2) {
                            flag = true;
                            break;;
                        }
                    }
                    //有锚点的话，添加到该点
                    if (flag) {
                        for (int nei_id:remove_core_adj_copy[l]) {
                            if (remove_core_adj_copy[nei_id].size() + kc_nei[nei_id] >= config.inputK) {
                                remove_core_adj_candi[nei_id].emplace_back(l);
                                remove_core_adj_candi[l].emplace_back(nei_id);
                            }
                        }
                    }
                }
            }
            if ((fout = fopen(outputFile.c_str(), "w+")) != NULL) {
                fprintf(fout, "%s,%s\n", "source", "target");
                for (int j = 0; j < remove_core_adj_candi.size(); ++j) {
                    for (int k = 0; k < remove_core_adj_candi[j].size(); ++k) {
                        fprintf(fout, "%d,%d\n", j, remove_core_adj_candi[j][k]);
                    }
                }
                fclose(fout);
            }
            outputFile = "node_" + outputFile;
            if ((fout = fopen(outputFile.c_str(), "w+")) != NULL) {
                fprintf(fout, "%s,%s,%s\n", "Id", "deg", "k_tag");
                INFO("Save node to file ", outputFile);
                for (int j = 0; j < remove_core_adj.size(); ++j) {
                    if (remove_core_adj_candi[j].size() > 0) {
                        fprintf(fout, "%d,%d,%d\n", j, kc_nei[j], k_tag_anchor[j]);
                    }
                }
                fclose(fout);
            }
            INFO("save graph finished");
        } else {
            INFO("Can not find output file, please test path!");
        }
    }

    bool cmp_cand_anc(const candidate_anc &a, const candidate_anc &b) {
        if (b.del_ancs == 0) return false;
        if (a.del_ancs == 0) return true;

        if (a.del_fol * 1.0 / a.del_ancs < whole_follow * 1.0 / whole_anchor &&
            b.del_fol * 1.0 / b.del_ancs < whole_follow * 1.0 / whole_anchor) {
            if (a.del_ancs < b.del_ancs)return true;//想要anc更多的，但是可能anc差的不多，fol差得多啊。
            if (a.del_ancs > b.del_ancs)return false;
            if (a.del_fol > b.del_fol) return true;//想要fol更少的
            if (a.del_fol < b.del_fol)return false;
        }
        //cout<<"equal"<<endl;
        if (a.del_fol * 1.0 / a.del_ancs > b.del_fol * 1.0 / b.del_ancs)return true;//想要fol/anc更小的
        if (a.del_fol * 1.0 / a.del_ancs < b.del_fol * 1.0 / b.del_ancs)return false;
        return a.id > b.id;
    }

    int cmp_cand_anc(const candidate_anc &a, const candidate_anc &b, const candidate_anc &c) {
        //0：a最大，1：b最大，2：c最大，-1没有合适的
        //对于锚点大于0且删除的followers不多的，只比较锚点。
        if (a.del_ancs == 0 && b.del_ancs == 0 && c.del_ancs == 0)return -1;
        if (cmp_cand_anc(a, b)) {
            if (cmp_cand_anc(b, c)) {
                return 2;
            } else {
                return 1;
            }
        } else {
            if (cmp_cand_anc(a, c)) {
                return 2;
            } else {
                return 0;
            }
        }
    }


};


#endif //GETUNDIRECTGRAPH_ANCHOREDKCORE_H
