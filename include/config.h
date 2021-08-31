#ifndef __CONFIG_H__
#define __CONFIG_H__


#ifdef WIN32
#define FILESEP "\\"
#else
#define FILESEP "/"
#endif

//#include <boost/progress.hpp>
#include <getopt.h>
#include <boost/timer/timer.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem.hpp>
#include <unordered_map>
#include <list>
using namespace boost;
using namespace boost::property_tree;

const double epsilon = 1.0e-7;
const int kTag_anchor_value = -10;
const int NUM_ITERATION = 500;//随机抽取的轮数

const string QUERY = "query";
const string GET_UNDIRECT_GRAPH = "get_undirect_graph";

const string EXACT = "exact";
const string RAND = "random";
const string RAND1 = "random1";
const string RAND2 = "random2";
const string DEGREE = "degree";
const string DEGREE1 = "degree1";
const string LKX = "clock";


//用作计时
const int WHOLE = 0;
const int READ_GRAPH = 1;
const int COMPUTE_KCORE = 2;
const int EXACT_QUERY = 3;
const int RAND_QUERY=4;
const int DEGREE_QUERY=5;
const int LKX_QUERY=6;
const int GET_NUM_KCORE_NEIGHBOR = 7;
const int PROCESS = 8;
const int COMPUTE_ANCHORS = 9;
const int MERGE_STEP = 10;
const int DIVIDE_STEP = 11;
const int SUPPORTER_COMPUTE = 12;
const int COMPUTE_KCORE_WITH_ANC = 13;

#ifdef WIN32
const string parent_folder = "../../";
#else
const string parent_folder = string("./") + FILESEP;
#endif

class Config {
public:
    string operation = QUERY;
    string graph_alias;
    string graph_location;

    string prefix = "../kcore_dataset/";
    string outfile = "result.txt";
    string randomfile = "random_result.txt";
    string randomfile1 = "random_result1.txt";
    string randomfile2 = "random_result2.txt";
    string degreefile = "degree_result.txt";
    string degreefile1 = "degree_result1.txt";
    string exactfile = "exact_result.txt";

    string exe_result_dir = parent_folder;//结果保存路径

    string get_graph_folder() {
        return prefix + graph_alias + FILESEP;
    }

    int inputK = 5;
    int inputB = 10;
    //double inputBratio = 0.2;
    double threshold = 0.2;

    int mode = 3;//0: 无合并, 1: 直接合并, 2: 直接合并+支持者合并, 3: 直接合并+支持者合并+邻居合并
    int supporter_num = 4;// 支持者合并的候选项
    int neighbor_num = 2;// 和邻居合并

    int new_graph_size = 500;

    string algo;

    void init_config(int argc, char *argv[]) {
        int opt; // getopt_long() 的返回值
        int digit_optind = 0; // 设置短参数类型及是否需要参数

        // 如果option_index非空，它指向的变量将记录当前找到参数符合long_opts里的
        // 第几个元素的描述，即是long_opts的下标值
        int option_index = 0;
        // 设置短参数类型及是否需要参数
        const char *optstring = "k:b:t:";

        // 设置长参数类型及其简写，比如 --reqarg <==>-r
        static struct option long_options[] = {
                {"dataset", required_argument, NULL, 'd'},
                {"prefix",  required_argument, NULL, 'p'},
                {"operation", required_argument, NULL,'o'},
                {"algo", required_argument, NULL,'a'},
                {"mode", required_argument, NULL,'m'},
                {"supporter_num", required_argument, NULL,'s'},
                {"neighbor_num", required_argument, NULL,'n'},
                {"new_graph_size", required_argument, NULL,'g'},
                {0, 0, 0, 0}  // 添加 {0, 0, 0, 0} 是为了防止输入空值
        };

        while ( (opt = getopt_long(argc,argv,optstring,long_options,&option_index)) != -1) {
            switch (opt) {
                case 'k':
                    inputK = atoi(optarg);
                    break;
                case 'b':
                    inputB = atoi(optarg);
                    break;
                case 't':
                    threshold = atof(optarg);
                    break;
                case 'd':
                    graph_alias = string(optarg);
                    break;
                case 'p':
                    prefix = string(optarg);
                    break;
                case 'o':
                    operation = string(optarg);
                    break;
                case 'a':
                    algo = string(optarg);
                    break;
                case 'm':
                    mode = atoi(optarg);
                    break;
                case 's':
                    supporter_num = atoi(optarg);
                    break;
                case 'n':
                    neighbor_num = atoi(optarg);
                    break;
                case 'g':
                    new_graph_size = atoi(optarg);
                    break;
            }
        }
        graph_location=get_graph_folder();
    }
    ptree get_data() {
        ptree data;
        data.put("graph_alias", graph_alias);
        data.put("operation", operation);
        data.put("inputK", inputK);
        data.put("inputB", inputB);
        data.put("threshold", threshold);
        data.put("mode",mode);
        data.put("supporter_num",supporter_num);
        data.put("neighbor_num",neighbor_num);
        //data.put("adjust_score",adjust_score);
        data.put("algo", algo);
        data.put("result-dir", exe_result_dir);
        return data;
    }
};

class Result {
public:
    int n;
    long long m;
    int kcore_num, deg_more_k;
    double follower_num;

    unordered_map<int, pair<double, double>> fol_anc_time;

    vector<int> anchors;
    vector<int> followers;
    std::chrono::steady_clock::time_point startTime;
    double whole_time_usage;

    double total_mem_usage;
    double total_time_usage;

    double merge_time;
    double merge_time_ratio;

    double divide_time;
    double divide_time_ratio;

    double supporter_time;
    double supporter_time_ratio;

    ptree get_data() {
        ptree data;
        data.put("n", n);
        data.put("m", m);
        data.put("kcore_num", kcore_num);
        data.put("follower_num", follower_num);

        data.put("whole time usage",whole_time_usage);
        data.put("total memory usage(MB)", total_mem_usage);

        data.put("total time usage(s)", total_time_usage);
        data.put("total time on merging anchors", merge_time);
        //data.put("total time on dividing anchors", divide_time);
        data.put("total time on supporter computation", supporter_time);


        data.put("total time ratio on merging anchors", merge_time_ratio);
        //data.put("total time ratio on dividing anchors", divide_time_ratio);
        data.put("total time ratio on supporter computation", supporter_time_ratio);

        return data;
    }

};

extern Config config;
extern Result result;

bool exists_test(const std::string &name);

void assert_file_exist(string desc, string name);

namespace Saver {
    static string get_current_time_str() {
        time_t rawtime;
        struct tm *timeinfo;
        char buffer[80];

        time(&rawtime);
        timeinfo = localtime(&rawtime);

        strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", timeinfo);
        std::string str(buffer);

        return str;

    }
    //设置保存的文件名
    static string get_time_path() {
        // using namespace boost::posix_time;
        // auto tm = second_clock::local_time();
        if(!boost::algorithm::ends_with(config.exe_result_dir, FILESEP))
            config.exe_result_dir += FILESEP;
        config.exe_result_dir += "execution/";
        if(!boost::filesystem::exists(config.exe_result_dir)){
            boost::filesystem::path dir(config.exe_result_dir);
            boost::filesystem::create_directories(dir);
        }

        string filename = config.graph_alias+"."+config.algo;

        filename += "k-" + to_string(config.inputK) + ".";
        filename += "t-" + to_string(config.inputB) + ".";
        filename += "mode-" + to_string(config.mode) + ".";
        filename += "s-"+to_string(config.supporter_num)+".";
        filename += "n-"+to_string(config.neighbor_num);

        return config.exe_result_dir + filename;
        // return config.exe_result_dir + to_iso_string(tm);
    }

    static ptree combine;

    static void init() {
        combine.put("start_time", get_current_time_str());
    }


    static void save_json(Config &config, Result &result, vector<string> args) {
        ofstream fout(get_time_path() + ".json");
        string command_line = "";
        for (int i = 1; i < args.size(); i++) {
            command_line += " " + args[i];
        }
        combine.put("end_time", get_current_time_str());
        combine.put("command_line", command_line);
        combine.put_child("config", config.get_data());
        combine.put_child("result", result.get_data());
        ptree timer;
        for (int i = 0; i < (int) Timer::timeUsed.size(); i++) {
            if (Timer::timeUsed[i] > 0) {
                timer.put(to_str(i), Timer::timeUsed[i] / TIMES_PER_SEC);
            }
        }
        combine.put_child("timer", timer);

        write_json(fout, combine, true);
    }
};


#endif
