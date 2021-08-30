//
// Created by 刘凯鑫 on 2020/1/19.
//
#include "../include/mylib.h"


template<class T>
string toStr(T t) {
    stringstream ss;
    ss << t;
    return ss.str();
}

string __n_variable(string t, int n) {
    t = t + ',';
    int i = 0;
    if (n) for (; i < SIZE(t) && n; i++) if (t[i] == ',') n--;
    n = i;
    for (; t[i] != ','; i++);
    t = t.substr((unsigned long) n, (unsigned long) (i - n));
    trim(t);
    if (t[0] == '"') return "";
    return t + "=";
}
//静态对象的初始化
int Counter::cnt[1000] = {0};

vector<double> Timer::timeUsed;
vector<string> Timer::timeUsedDesc;
