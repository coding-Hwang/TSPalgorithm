#pragma once
#include <cmath>

using namespace std;

const int INF = 1e9;

// point structure
struct Point {
    double x, y;
};

// 두 점 사이의 거리
inline double dist(Point a, Point b) {
    return hypot(a.x - b.x, a.y - b.y);
}

struct TSPResult {
    vector<int> path;
    double total_cost;
};