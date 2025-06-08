#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

#include "Held-Karp.h"

using namespace std;

// Held-Karp DP
TSPResult held_karp(const vector<Point>& points) {
    int n = points.size();
    int full = 1 << n; // 2^n -> 모든 부분 집합 (비트 마스크 이용)
    vector<vector<double>> dp(full, vector<double>(n, INF));
    // dp[mask][u] : mask 부분 집합에 방문 + 마지막 도시가 u

    vector<vector<int>> parent(full, vector<int>(n, -1)); // 역추적용

    dp[1][0] = 0; // start at node 0

    for (int mask = 1; mask < full; ++mask) { // 부분 집합
        for (int u = 0; u < n; ++u) { // 마지막 방문 도시 후보
            if (!(mask & (1 << u))) continue; // u 포함 시에만 실행
            for (int v = 0; v < n; ++v) { // 다음 방문 도시 후보
                if (mask & (1 << v)) continue; // 포함 안되어있을 때만
                int nextMask = mask | (1 << v); // v 포함한 새로운 집합

                double newCost = dp[mask][u] + dist(points[u], points[v]);
                if (newCost < dp[nextMask][v]) {
                    dp[nextMask][v] = newCost; // 최소 비용 갱신
                    parent[nextMask][v] = u; // v 로 가기 전에 u 였음
                }
                // 마지막 도시 u 별 v 로 방문하는 경우의 수를 갱신
            }
        }
    }

     // 최종 비용과 마지막 도시 u
    double minCost = INF;
    int lastCity = -1;
    for (int u = 1; u < n; ++u) {
        double cost = dp[full - 1][u] + dist(points[u], points[0]);
        if (cost < minCost) {
            minCost = cost;
            lastCity = u;
        }
    }

    // 경로 재구성 (역추적)
    vector<int> path;
    int mask = full - 1;
    int curr = lastCity;
    while (curr != -1) {
        path.push_back(curr);
        int temp = parent[mask][curr];
        mask ^= (1 << curr);
        curr = temp;
    }
    path.push_back(0); // 시작 도시로 돌아감
    reverse(path.begin(), path.end());

    return { path, minCost };
}