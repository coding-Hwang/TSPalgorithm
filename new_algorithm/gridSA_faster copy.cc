/*  grid_sa_tsp.cpp  ───────────────────────────────────────────
 *  Adaptive-grid Simulated-Annealing TSP + OpenMP acceleration
 *  - 격자 한 칸에 최대 280개까지
 *  - 격자별 SA를 병렬 수행 (#pragma omp parallel for)
 *  - 최종 메타 MST-2-approx + 전역 2-opt
 *  ----------------------------------------------------------*/

#include "gridSA.h"

#include <bits/stdc++.h>
#include <omp.h>
using namespace std;

/* MST & Euler 헬퍼 (기존 구현 사용) */
using Adj = vector<vector<int>>;
Adj buildMST(const vector<Point>& P);                // 프림/크루스컬 등
void dfs_eulerian(int v, const Adj& adj,
                  vector<vector<bool>>& used,
                  vector<int>& out);                 // Euler DFS

/* ---------- 2-opt ---------- */
double tourLength(const vector<Point>& P,
                  const vector<int>& path)
{
    double len = 0;
    for (int i = 0; i + 1 < (int)path.size(); ++i)
        len += dist(P[path[i]], P[path[i + 1]]);
    return len;
}

bool twoOptOnce(const vector<Point>& P, vector<int>& path)
{
    int n = (int)path.size();
    if (n < 5) return false;
    for (int i = 1; i < n - 2; ++i)
        for (int k = i + 1; k < n - 1; ++k) {
            double delta =
                dist(P[path[i-1]], P[path[k]]) +
                dist(P[path[i  ]], P[path[k+1]]) -
                dist(P[path[i-1]], P[path[i  ]]) -
                dist(P[path[k  ]], P[path[k+1]]);
            if (delta < -1e-9) {
                reverse(path.begin() + i, path.begin() + k + 1);
                return true;
            }
        }
    return false;
}
void twoOptFull(const vector<Point>& P,
                vector<int>& path,
                int maxIter = 6)
{
    for (int it = 0; it < maxIter; ++it)
        if (!twoOptOnce(P, path)) break;
}

/* ---------- 스레드-세이프 Simulated Annealing ---------- */
vector<int> simulatedAnnealing(const vector<Point>& P,
                               const vector<int>& seed,
                               double T0   = 1000.0,
                               double cool = 0.995,
                               double T_min = 1e-6)
{
    if (seed.size() <= 4) return seed;

    /* 스레드별 RNG */
    thread_local std::mt19937 rng(
        std::random_device{}() ^
        (std::hash<std::thread::id>()(std::this_thread::get_id())));
    std::uniform_real_distribution<double> U(0.0, 1.0);

    auto pickIdx = [&](int n) { return 1 + int(rng() % (n - 2)); };

    vector<int> best = seed, cur = seed;
    double bestLen = tourLength(P, best),
           curLen  = bestLen;

    double T = T0;
    while (T > T_min) {
        for (int rep = 0; rep < 50; ++rep) {
            vector<int> nxt = cur;
            int a = pickIdx(nxt.size()), b = pickIdx(nxt.size());
            if (a > b) swap(a, b);
            reverse(nxt.begin() + a, nxt.begin() + b);

            double nxtLen = tourLength(P, nxt);
            double dE = nxtLen - curLen;
            if (dE < 0 || U(rng) < exp(-dE / T)) {
                cur.swap(nxt);  curLen = nxtLen;
                if (curLen < bestLen) { best = cur; bestLen = curLen; }
            }
        }
        T *= cool;
    }
    return best;
}

/* ---------- MST 기반 2-approx ---------- */
vector<int> mst2Approx(const vector<Point>& P)
{
    Adj adj = buildMST(P);
    int n = (int)P.size();
    vector<vector<bool>> used(n, vector<bool>(n,false));
    vector<int> eul;
    dfs_eulerian(0, adj, used, eul);
    reverse(eul.begin(), eul.end());

    vector<int> seen(n), out;
    for (int v : eul) if (!seen[v]) { out.push_back(v); seen[v]=1; }
    out.push_back(0);
    return out;
}

/* ---------- 재귀 격자 분할 ---------- */
void subdivide(const vector<Point>& P, const vector<int>& idx,
               double minX, double minY,
               double maxX, double maxY,
               vector<vector<int>>& cells,
               const int LIMIT = 280)
{
    if (idx.empty()) return;
    if ((int)idx.size() <= LIMIT) { cells.push_back(idx); return; }

    double midX = 0.5*(minX+maxX), midY = 0.5*(minY+maxY);
    vector<int> q[4];
    for (int v : idx) {
        bool right = P[v].x >= midX;
        bool top   = P[v].y >= midY;
        int k = (top?2:0) + (right?1:0);
        q[k].push_back(v);
    }
    subdivide(P, q[0], minX, minY, midX, midY, cells, LIMIT); // 좌하
    subdivide(P, q[1], midX, minY, maxX, midY, cells, LIMIT); // 우하
    subdivide(P, q[2], minX, midY, midX, maxY, cells, LIMIT); // 좌상
    subdivide(P, q[3], midX, midY, maxX, maxY, cells, LIMIT); // 우상
}

pair<vector<vector<int>>, int>
splitGrids(const vector<Point>& P, int LIMIT = 280)
{
    double minX =  numeric_limits<double>::max(),
           minY =  numeric_limits<double>::max(),
           maxX = -numeric_limits<double>::max(),
           maxY = -numeric_limits<double>::max();
    for (auto& p : P) {
        minX = min(minX, p.x); maxX = max(maxX, p.x);
        minY = min(minY, p.y); maxY = max(maxY, p.y);
    }
    vector<int> root(P.size()); iota(root.begin(), root.end(), 0);

    vector<vector<int>> cells;
    subdivide(P, root, minX, minY, maxX, maxY, cells, LIMIT);

    return {cells, (int)cells.size()};
}

/* ---------- 최종 파이프라인 ---------- */
struct TSPResult { vector<int> path; double cost; };

TSPResult gridSA_fast(const vector<Point>& P)
{
    /* 0) 작은 인스턴스는 단일 SA */
    if (P.size() <= 280) {
        vector<int> seed(P.size()+1);
        iota(seed.begin(), seed.end()-1, 0);
        seed.back() = 0;
        auto best = simulatedAnnealing(P, seed);
        return {best, tourLength(P,best)};
    }

    /* 1) 격자 분할 */
    auto [grid, G] = splitGrids(P);
    vector<vector<int>> sub(G);

    /* 2) 격자 내부 SA ─ 병렬 */
#pragma omp parallel for schedule(dynamic) default(none) \
        shared(P, grid, sub, G)
    for (int g = 0; g < G; ++g) {
        const auto& idx = grid[g];
        if (idx.empty()) continue;

        if (idx.size() == 1) { sub[g] = {idx[0]}; continue; }

        vector<int> seed(idx.begin(), idx.end());
        seed.push_back(idx[0]);
        sub[g] = simulatedAnnealing(P, seed);
    }

    /* 3) 메타 그래프(centroid) */
    vector<Point> cent; vector<int> gid;
    for (int g = 0; g < G; ++g)
        if (!sub[g].empty()) {
            double sx=0, sy=0;
            for (int v: grid[g]) { sx+=P[v].x; sy+=P[v].y; }
            cent.push_back({sx/grid[g].size(), sy/grid[g].size()});
            gid.push_back(g);
        }

    /* 4) 메타 MST-2-approx 순서 */
    vector<int> meta = mst2Approx(cent);

    /* 5) 투어 스티칭 */
    vector<int> total;
    for (int m : meta) {
        int g = gid[m];
        const auto& subP = sub[g];
        if (subP.empty()) continue;

        if (total.empty()) total = subP;
        else {
            if (total.back() == subP.front()) total.pop_back();
            bool closed = subP.size()>=2 && subP.front()==subP.back();
            total.insert(total.end(),
                         subP.begin(),
                         closed ? subP.end()-1 : subP.end());
        }
    }
    if (!total.empty() && total.front()!=total.back())
        total.push_back(total.front());

    /* 6) 전역 2-opt */
    twoOptFull(P, total, 8);
    return {total, tourLength(P,total)};
}
