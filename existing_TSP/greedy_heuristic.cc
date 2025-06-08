#include <iostream>
#include <vector>
#include <cmath>
#include <limits>

#include "greedy_heuristic.h"

using namespace std;


// 가장 가까운 미방문 도시 찾기
int findNearest(int current, const vector<Point>& points, const vector<bool>& visited) {
    int nearest = -1;
    double minDist = INF;

    for (int i = 0; i < points.size(); ++i) {
        if (!visited[i]) {
            double d = dist(points[current], points[i]);
            if (d < minDist) {
                minDist = d;
                nearest = i;
            }
        }
    }
    return nearest;
}

// Greedy TSP 알고리즘
TSPResult greedy(const vector<Point>& points) {
    int n = points.size();
    vector<bool> visited(n, false);
    int current = 0;
    visited[current] = true;

    vector<int> path;
    path.push_back(current);

    double totalCost = 0;

    for (int i = 1; i < n; ++i) {
        int next = findNearest(current, points, visited);
        if (next == -1) break;

        totalCost += dist(points[current], points[next]);
        visited[next] = true;
        path.push_back(next);
        current = next;
    }

    // 마지막 -> 시작
    totalCost += dist(points[current], points[0]);
    path.push_back(0);

    return {path, totalCost};
}
