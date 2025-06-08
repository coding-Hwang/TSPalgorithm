// MST-based 2-approximation algorithm for TSP
#include <iostream>
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

#include "MST-based.h"


using namespace std;

// 1. build MST
//  - prim's algorithm
vector<vector<int>> buildMST(const vector<Point>& points) {
    int n = points.size();
    vector<double> minDist(n, INF); // 각 노드까지의 최소 거리 (여기다가 다 저장)
    vector<int> parent(n, -1); // 각 노드의 부모 노드드
    vector<bool> inMST(n, false); // MST 포함 여부 (여기다가 다 저장)
    minDist[0] = 0; // 시작 노드의 거리는 0으로 초기화

    for (int i = 0; i < n; ++i) {
        int u = -1; // input 에서 선택된 정점의 인덱스스
        double minVal = INF; // 일단 큰걸로 잡고 시작
        for (int j = 0; j < n; ++j) { 
            // MST에 포함되지 않은 노드 중에서 가장 작은 노드 찾기
            if (!inMST[j] && minDist[j] < minVal) {
                minVal = minDist[j];
                u = j;
            }
        }
        inMST[u] = true; // 위에 찾은 노드 MST에 포함함
        for (int v = 0; v < n; ++v) { // points[u]와 연결된 모든 노드의 거리 저장
            double d = dist(points[u], points[v]); 
            if (!inMST[v] && d < minDist[v]) {
                minDist[v] = d;
                parent[v] = u;
            } // 다음 루프에서 최소 엣지 찾기 위해 minDist 업데이트 -> 계속 덮어 씌어짐짐
        }
    }

    // Convert to adjacency list
    vector<vector<int>> adj(n);
    for (int i = 1; i < n; ++i) {
        adj[i].push_back(parent[i]); // 자식 : 부모 추가 (인접 노드)
        adj[parent[i]].push_back(i); // 부모 : 자식 추가 (인접 노드드)
    }
    return adj;
}


// 2. DFS
// - MST에서 DFS (eulerian tour) 로 경로 찾기

void dfs_eulerian(int u,
                  const vector<vector<int>>& adj,
                  vector<vector<bool>>& used,
                  vector<int>& path) {
    for (int v : adj[u]) {
        if (used[u][v]) continue;           // 이미 방문한 간선이면 skip
        used[u][v] = used[v][u] = true;     // 간선 방문 처리 (무방향이므로 양방향 체크)
        dfs_eulerian(v, adj, used, path);  // 다음 정점으로 이동
    }
    path.push_back(u); // 후위 순회로 경로 저장
}


vector<int> shortcut (int n, vector<int> path) {
    vector<bool> seen(n, false);
    vector<int> shortcut_path;
    for (int v : path) {
        if (!seen[v]) {
            shortcut_path.push_back(v);
            seen[v] = true;
        }
    }
    return shortcut_path;
}

double calculate_total_distance(const vector<int>& path, const vector<Point>& points) {
    double total = 0.0;
    for (int i = 0; i < path.size() - 1; ++i) {
        total += dist(points[path[i]], points[path[i + 1]]);
    }
    return total;
}


// 2-approximation TSP
TSPResult mst_2approx(const vector<Point>& points) {
    vector<vector<int>> mst = buildMST(points); // 1. MST 생성
    // 2. DFS로 경로 찾기 (MST에 대한 오일러 경로)
    // 오일러 패스 경로 저장
    int n = points.size();
    vector<vector<bool>> used_edges(n, vector<bool>(n, false));
    vector<int> path;

    dfs_eulerian(0, mst, used_edges, path);
    reverse(path.begin(), path.end());

    vector<int> final_path = shortcut(n, path); // 중복 제거

    final_path.push_back(path[0]); // 시작점으로 돌아가기 위해 시작점 추가
    double cost = calculate_total_distance(final_path, points); 

    return {final_path, cost}; 
}


// ref - https://gazelle-and-cs.tistory.com/18

