#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
using namespace std;

#include "insertion_method.h"

// 경로의 총 거리 계산
double total_cost(const vector<int>& path, const vector<Point>& points) {
    double cost = 0;
    for (int i = 0; i < path.size() - 1; ++i) {
        cost += dist(points[path[i]], points[path[i + 1]]);
    }
    return cost;
}

// 삽입 위치 계산
int bestInsertionPosition(const vector<int>& path, int city, const vector<Point>& points) {
    int bestPos = -1;
    double minIncrease = numeric_limits<double>::max();

    for (int i = 0; i < path.size() - 1; ++i) {
        int u = path[i];
        int v = path[i + 1];
        double increase = dist(points[u], points[city]) + dist(points[city], points[v]) - dist(points[u], points[v]);
        if (increase < minIncrease) {
            minIncrease = increase;
            bestPos = i + 1;
        }
    }

    return bestPos;
}

// Nearest Insertion 알고리즘
TSPResult nearestInsertion(const vector<Point>& points) {
    int n = points.size();
    vector<bool> visited(n, false);

    // 시작점
    int start = 0;
    visited[start] = true;

    // 시작점에서 가장 가까운 점 찾기
    int nearest = -1;
    double minDist = numeric_limits<double>::max();
    for (int i = 0; i < n; ++i) {
        if (i != start) {
            double d = dist(points[start], points[i]);
            if (d < minDist) {
                minDist = d;
                nearest = i;
            }
        }
    }

    visited[nearest] = true;
    vector<int> path = { start, nearest, start };

    // 나머지 도시 삽입
    for (int count = 2; count < n; ++count) {
        int nextCity = -1;
        double nearestDist = numeric_limits<double>::max();
        for (int city = 0; city < n; ++city) {
            if (visited[city]) continue;
            for (int inPath : path) {
                double d = dist(points[city], points[inPath]);
                if (d < nearestDist) {
                    nearestDist = d;
                    nextCity = city;
                }
            }
        }

        int pos = bestInsertionPosition(path, nextCity, points);
        path.insert(path.begin() + pos, nextCity);
        visited[nextCity] = true;
    }

    double cost = total_cost(path, points);
    return {path, cost};
}
