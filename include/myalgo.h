//
// Created by 刘凯鑫 on 2020/1/24.
//

#ifndef GETUNDIRECTGRAPH_MYALGO_H
#define GETUNDIRECTGRAPH_MYALGO_H

#include "mylib.h"
#include "config.h"

bool cmp_pair(const pair<int, double> &p1, const pair<int, double> &p2);

//降序排列
bool cmp_pair_dec(const pair<int, int> &p1, const pair<int, int> &p2);

bool compare_state(vector<pair<int, int>> &state_vec, int node, int function_counter, int ag_id);

struct cmpPair {
    bool operator()(const pair<int, int> &p1, const pair<int, int> &p2) {
        return p1.second > p2.second; //second的小值优先
    }
};
struct cmpPair2 {
    bool operator()(const pair<int, int> &p1, const pair<int, int> &p2) {
        return p1.second < p2.second; //second的大值优先
    }
};

vector<int> s_v(const set<int> &tmp);

vector<int> vec1_vec2(const vector<int> &vec1, const vector<int> &vec2);

unordered_map<int, int> get_different_nodes(const vector<int> &vec1, const vector<int> &vec2);

struct fol_deg_merge {
    int fol_id, score, fol_deg;
};

int compare_double(const double &score1, const double &score2,double thre=epsilon);

struct ans_index {
    int id = 0, anc_num = 0, fol_num = 0;
    double score = 0, stability = 0;

    ans_index(int x, int y, double z, int t, double st) {
        id = x;
        anc_num = y;
        score = z;
        fol_num = t;
        stability = st;
    }


    friend bool operator<(const ans_index &a, const ans_index &b) {
        if (compare_double(a.score, b.score) == 1)return true;//分数大的真
        if (compare_double(a.score, b.score) == -1)return false;
        if (a.anc_num > b.anc_num) return true;
        if (a.anc_num < b.anc_num) return false;
        //如果anc_num相等
        //如果anc_num和分数都相等
        if (a.fol_num > b.fol_num) return false;
        if (a.fol_num < b.fol_num)return true;

        if (compare_double(a.stability, b.stability) == 1)return true;
        if (compare_double(a.stability, b.stability) == -1)return false;

        return a.id < b.id;
    }
};

struct ans_index_big_first {
    int id = 0, anc_num = 0, fol_num = 0;
    double score = 0, stability;

    ans_index_big_first(int i, int a, double sc, int f, double st) {
        id = i;
        anc_num = a;
        score = sc;
        fol_num = f;
        stability = st;
    }

    friend bool operator<(const ans_index_big_first &a, const ans_index_big_first &b) {
        //如果anc_num和分数都相等
        if (compare_double(a.score, b.score) == 1)return true;
        if (compare_double(a.score, b.score) == -1)return false;
        if (a.anc_num > b.anc_num) return true;
        if (a.anc_num < b.anc_num) return false;
        //如果anc_num相等
        //如果anc_num和分数都相等
        if (a.fol_num > b.fol_num) return false;
        if (a.fol_num < b.fol_num)return true;

        if (compare_double(a.stability, b.stability) == 1)return true;
        if (compare_double(a.stability, b.stability) == -1)return false;

        return a.id > b.id;
    }
};

struct candidate_anc {
    int id = -1, del_ancs = 0, del_fol = 0;
    double affected_num = 0, score2 = 0, score3 = 0;

    candidate_anc(int id, int del_ancs, int del_fol, double affected_num, double score2 = 0, double score3 = 0) {
        this->id = id;
        this->del_ancs = del_ancs;
        this->del_fol = del_fol;
        this->affected_num = affected_num;
        this->score2 = score2;
        this->score3 = score3;
    };
};
struct deleted_node {
    int id = 0, deg = 0, group = 0;
    //double score = 0, stability = 0;

    deleted_node(int id, int deg, int group) {
        this->id = id;
        this->deg = deg;
        this->group = group;
    }


    friend bool operator<(const deleted_node &a, const deleted_node &b) {
        if (a.deg<b.deg) return false;
        if (a.deg>b.deg) return true;
        return a.id < b.id;
    }
};
bool cmp_fol_deg_merge(const fol_deg_merge &a, const fol_deg_merge &b);


inline void split_line() {
    INFO("-----------------------------");
}

static void display_time_usage() {
    double whole_time = Timer::used(COMPUTE_KCORE) + Timer::used(GET_NUM_KCORE_NEIGHBOR) + Timer::used(PROCESS);
    cout << "Total cost  (s): " << whole_time + Timer::used(READ_GRAPH) << endl;
    cout << "Total cost exclude read graph (s): " << whole_time << endl;
    cout << Timer::used(COMPUTE_KCORE) * 100.0 / whole_time << "%" << " for  k-core nodes computation cost" << endl;
    cout << Timer::used(GET_NUM_KCORE_NEIGHBOR) * 100.0 / whole_time << "%"
         << " for numbers of k-core neighbor computation cost" << endl;
    cout << (Timer::used(PROCESS) - Timer::used(COMPUTE_ANCHORS)) * 100.0 / whole_time << "%"
         << " for index construction cost" << endl;
    cout << Timer::used(COMPUTE_ANCHORS) * 100.0 / whole_time << "%" << " for queries cost" << endl;

    split_line();
    //cout << "Average query time (s):" << Timer::used(COMPUTE_ANCHORS) / config.inputBratio << endl;
    //cout << "Memory usage (MB):" << get_proc_memory()/1000.0 << endl << endl;
}

inline int get_value(const unordered_map<int, int> &m, int key) {
    if (m.find(key) != m.end()) {
        return m.at(key);
    }
    return -1;
}


vector<int> merge_vec(const vector<int> &vec1, const vector<int> &vec2, vector<int> &vec_erase);

vector<int> add_vec(vector<int> &vec1, const vector<int> &vec2);

queue<int> push_neighbor_queue(vector<int> neighbors);

queue<int> push_neighbor_queue(set<int> neighbors);


int vec_find(const vector<int> &vec, int node);

bool vec_erase(vector<int> &vec, int index);

//对每个点进行累加或减得分
void add_score_to_nodes(vector<int> anchores, unordered_map<int, double> &node_score, double score);

void add_score_to_nodes(set<int> anchores, unordered_map<int, double> &node_score, double score);

set<int> del_nodes_safely(set<int> deleted_vec, unordered_map<int, vector<int>> &g);

vector<int> del_node(int deleted_node, unordered_map<int, vector<int>> &g);


//得到去掉k-core的子图
unordered_map<int, vector<int>> get_margin_subgraph(const vector<vector<int>> &graph, const vector<int> &kTag);

//持续删除度在[lbound,ubound)中的点
vector<int>
continuous_del_branch_node(unordered_map<int, vector<int>> &g, queue<int> active_nodes,
                           const unordered_map<int, int> &node_kcoreNei, int ubound = config.inputK);

//持续删除度在[lbound,ubound)中的点，并返回删除的边
vector<int>
continuous_del_branch_node(unordered_map<int, vector<int>> &g, queue<int> active_nodes,
                           const unordered_map<int, int> &node_kcoreNei, unordered_map<int, vector<int>> &del_adj_list,
                           int ubound = config.inputK);


///去掉度小于K的点
unordered_map<int, vector<int>>
refine_margin_graph(const unordered_map<int, vector<int>> &graph, const unordered_map<int, int> &node_kcoreNei);

///////将图转为多个子图存储，然后去掉枝丫节点
vector<unordered_map<int, vector<int>>> get_disjoint_subgraphs(const unordered_map<int, vector<int>> &graph);

///为每个子图计算相应的node_kcoreNei
unordered_map<int, int>
get_subgraphs_node_kcoreNei(unordered_map<int, vector<int>> subgraph, unordered_map<int, int> node_kcoreNei);

/*
void compute_all_candidate_anchores_sets(vector<vector<int>> &result, int &min_size, vector<int> alive_raw_anchors,
                                         vector<int> new_anchores,
                                         unordered_map<int, vector<int>> &node_nei,
                                         unordered_map<int, vector<int>> &nei_node,
                                         unordered_map<int, vector<int>> subg,
                                         const unordered_map<int, int> &node_kcoreNei,
                                         const unordered_map<int, vector<int>> &g_margin);
                                         */
void compute_degree_candidate_anchores_sets(vector<vector<int>> &result, int &min_size, vector<int> &alive_raw_anchors,
                                            vector<int> new_anchores,
                                            unordered_map<int, vector<int>> &node_nei,
                                            unordered_map<int, vector<int>> &nei_node,
                                            unordered_map<int, vector<int>> subg,
                                            const unordered_map<int, int> &node_kcoreNei,
                                            const unordered_map<int, vector<int>> &g_margin);

#endif //GETUNDIRECTGRAPH_MYALGO_H
