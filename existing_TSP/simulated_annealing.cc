// 완전 랜덤

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <random>

#include "simulated_annealing.h"


using namespace std;



// 경로의 총 길이 계산
double totalCost(const vector<int>& path, const vector<Point>& points) {
    double cost = 0;
    for (int i = 0; i < path.size() - 1; ++i) {
        cost += dist(points[path[i]], points[path[i + 1]]);
    }
    cost += dist(points[path.back()], points[path[0]]); // 원점 복귀
    return cost;
}



// 랜덤으로 경로 초기화
vector<int> createInitialPath(int n) {
    vector<int> path(n - 1);  // 마지막 노드는 중복되므로 n-1개만 셔플
    for (int i = 0; i < n - 1; ++i) path[i] = i + 1; // 1부터 n-1까지
    random_shuffle(path.begin(), path.end());

    path.insert(path.begin(), 0);  // 시작 노드 0 삽입
    path.push_back(0);    // 끝 노드도 0 삽입 (순환)
    return path;
}

// 경로 일부 swap (2-opt)
vector<int> getNeighbor(const vector<int>& path) {
    vector<int> newPath = path;
    int n = path.size();

    // 인덱스 1 ~ n-2 (시작과 끝은 0으로 고정)
    int i = 1 + rand() % (n - 2);
    int j = 1 + rand() % (n - 2);
    if (i > j) swap(i, j);

    reverse(newPath.begin() + i, newPath.begin() + j + 1);
    return newPath;
}

// Simulated Annealing 본체
TSPResult simulatedAnnealing(const vector<Point>& points) {
    int n = points.size();
    vector<int> currPath = createInitialPath(n); // 초기 경로 설정정
    double currCost = totalCost(currPath, points); // 초기 경로의 비용용
    vector<int> bestPath = currPath; // 최적 경로 초기화
    double bestCost = currCost;

    double T = 1000.0;           // 초기 온도 (높을수록 나쁜 해도 쉽게 수락락)
    double T_min = 1e-6;         // 최소 온도
    double alpha = 0.995;        // 냉각률

    while (T > T_min) {
        for (int iter = 0; iter < 100; ++iter) {
            vector<int> newPath = getNeighbor(currPath);
            double newCost = totalCost(newPath, points);
            double delta = newCost - currCost; // 새로운해가 얼마나 더 나빠졌는지지

            if (delta < 0 || exp(-delta / T) > (double)rand() / RAND_MAX) {
                currPath = newPath;
                currCost = newCost;
                if (currCost < bestCost) {
                    bestPath = currPath;
                    bestCost = currCost;
                }
            }
        }
        T *= alpha;
    }

    return {bestPath, bestCost};
}
/*
1. 초기 경로생성
2. 경로 살짝 바꿔 -> neighbor swap
3. 바꾼 해의 cost delta 계산
4. 좋아지면 수락 , 나빠져도 exp(-delta / T) 확률로 수락
5. T는 점점 감소 (T * alpha) -> 매우 작아질때까지 반복복

*/


// O(K L n)
// n = 도시의 개수, K = 온도 단계 수 (T > T_min), L = 각 온도에서의 반복 횟수수