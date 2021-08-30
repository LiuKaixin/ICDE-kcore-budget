#include <iostream>
#include <string>
#include <vector>
#include "./include/graph.h"
#include "./include/query.h"

using namespace std;
/*
 * 使用节点数大小的集合向量存储！
 * 可以对每次输入的两个数字排序，小的在前面。
 * 如果遇到初始不是0或者1，当做字符串。
 *
struct Graph{
    string name;
    bool index_id= false;//是否进行对节点id进行转化
    int n=0,m=0;
    int final_n=0, final_m=0;
    int min_node=-1,max_node=-1;
    Graph(string name){
        this->name=name;
    };
};
 */
int main(int argc, char *argv[]) {
    ios::sync_with_stdio(false);
    program_start(argc, argv);//输出命令内容
    Saver::init();
    srand(time(NULL));
    config.init_config(argc, argv);

    INFO(config.operation);
    if (config.operation == QUERY){
        vector<string> possibleAlgo = {EXACT, RAND, RAND1, RAND2, DEGREE, DEGREE1, LKX};
        if (find(possibleAlgo.begin(), possibleAlgo.end(), config.algo) == possibleAlgo.end()) {
            INFO("Wrong algo param: ", config.algo);
            exit(1);
        }
        INFO(config.algo);

        Graph graph(config.graph_location);
        graph.detect_kcore();

        INFO(graph.kcore_size);
        INFO(result.deg_more_k);
        query(graph);
    } else if (config.operation == GET_UNDIRECT_GRAPH) {
        //自己用直接写死好了。
        //vector<string> graphs = {"ego-facebook","loc-Brightkite","loc-gowalla","web-flickr","com-youtube","com-dblp","com-lj","pokec","orkut"};
        vector<string> graphs = {"dataset2"};
        for (string graph_file:graphs) {
            Graph graph;
            graph.convert_to_undirected_graph(graph_file);
            INFO(graph_file);
        }
    } else if (config.operation == "test") {
        Graph graph(config.graph_location);
        for (int k = 5; k < 35; k=k+5) {
            for (int i = 0; i < 10; ++i) {
                do {
                    graph.generate_subgraph_bfs();
                    config.inputK=k;
                    graph.detect_kcore();
                } while (graph.kcore_size==0);
                INFO(graph.kcore_size);
                string  path="./datasets/"+config.graph_alias+".k."+to_str(k)+".No."+to_str(i);
                mkdir(path.c_str(), 0775);
                graph.save_subgraph(path);
            }
        }

    }
    //facebook 有向--operation query --prefix ../kcore_dataset/ --dataset ego-facebook --algo lkx -k 30 -b 2

    //输出属性min_node largset_node n m

    Timer::show();
    if(config.operation == QUERY){
        Counter::show();
        auto args = combine_args(argc, argv);
        Saver::save_json(config, result, args);
    }
    return 0;
}