//
// Created by 刘凯鑫 on 2020/1/24.
//

#ifndef GETUNDIRECTGRAPH_ANCHORGROUP_H
#define GETUNDIRECTGRAPH_ANCHORGROUP_H

#include "mylib.h"
#include "config.h"
#include "myalgo.h"

class AnchorGroup {
public:
    int anc;
    vector<int> supporters, new_supporters;
    double score = 0, old_score;//根据包含的第一个删除的点做调整。
    double stability = 0; //好像是节点的平均度
    int state = 0;//0代表还没有排序过，-1代表被删除，1代表已经在排序里了。
    int init_supp_num = -1;//初始时支持者的数量，对于多锚点的锚组，支持者变更太多则拆分
    bool tag = true;
    //unordered_map<int, int> node_anc_deg,node_anc_deg_backup;


    //AnchorGroup(const set<int> &anchored_nodes) { possible_anchored_sets.emplace_back(anchored_nodes); }

    AnchorGroup(int anchor) : anc(anchor) {}

    AnchorGroup() {}

    bool static cmpAnchor(const AnchorGroup &set1, const AnchorGroup &set2) {
        return set1.score < set2.score;
    }

    void print() {
        //printf("\n");
        printf("Anchore:\t");
        printf("%d\t", anc);
        //printf("\tScore:\t%f\tAdjust_score:\t%f\n", score,adjust_score);
        printf("\tScore:\t%f\n", score);
        //INFO(new_supporters);
        printf("supporters:\t");
        for (int node: supporters) {
            printf("%d\t", node);
        }
        printf("\n\n");
    }
};
#endif //GETUNDIRECTGRAPH_ANCHORGROUP_H
