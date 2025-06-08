#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>

#include "branch_bound.h"
#include "__point.h"
using namespace std;


struct Node {
    vector<int> path;
    vector<bool> visited;
    int level;
    double cost;
    double bound;

    bool operator>(const Node& other) const {
        return bound > other.bound;
    }
};

// 유클리드 거리 행렬 만들기
vector<vector<double>> makeDistanceMatrix(const vector<Point>& points) {
    int n = points.size();
    vector<vector<double>> mat(n, vector<double>(n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            mat[i][j] = dist(points[i], points[j]);
    return mat;
}

// Bound 계산: 현재까지의 경로 비용 + 남은 노드에 대한 최소 거리 추정
double calculateBound(const Node& node, const vector<vector<double>>& costMatrix) {
    double bound = node.cost;
    int n = costMatrix.size();

    for (int i = 0; i < n; ++i) {
        if (!node.visited[i]) {
            double minEdge = numeric_limits<double>::max();
            for (int j = 0; j < n; ++j) {
                if (i != j && !node.visited[j]) {
                    minEdge = min(minEdge, costMatrix[i][j]);
                }
            }
            bound += (minEdge == numeric_limits<double>::max() ? 0 : minEdge);
        }
    }
    return bound;
}

double branchAndBoundTSP(const vector<Point>& points) {
    int n = points.size();
    vector<vector<double>> costMatrix = makeDistanceMatrix(points);

    priority_queue<Node, vector<Node>, greater<Node>> pq;

    Node root;
    root.path.push_back(0);
    root.visited = vector<bool>(n, false);
    root.visited[0] = true;
    root.level = 1;
    root.cost = 0;
    root.bound = calculateBound(root, costMatrix);

    pq.push(root);

    double minCost = numeric_limits<double>::max();
    vector<int> bestPath;

    while (!pq.empty()) {
        Node current = pq.top();
        pq.pop();

        if (current.bound >= minCost) continue;

        if (current.level == n) {
            double finalCost = current.cost + costMatrix[current.path.back()][0];
            if (finalCost < minCost) {
                minCost = finalCost;
                bestPath = current.path;
                bestPath.push_back(0);
            }
            continue;
        }

        for (int i = 0; i < n; ++i) {
            if (!current.visited[i]) {
                Node next = current;
                next.path.push_back(i);
                next.visited[i] = true;
                next.level = current.level + 1;
                next.cost += costMatrix[current.path.back()][i];
                next.bound = calculateBound(next, costMatrix);

                if (next.bound < minCost)
                    pq.push(next);
            }
        }
    }

    cout << "Best TSP path (Branch and Bound): ";
    for (int city : bestPath)
        cout << city << " ";
    cout << "\nTotal cost: " << minCost << endl;

    return minCost;
}