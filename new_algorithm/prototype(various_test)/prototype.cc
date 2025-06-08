#include "me.h"

// tsp_grid_sa.cpp  (C++17)
#include <bits/stdc++.h>
using namespace std;



/* ---------- 2-opt ---------- */
double tourLength(const vector<Point>& P, const vector<int>& path) {
    double len = 0;
    for (int i = 0; i + 1 < (int)path.size(); ++i)
        len += dist(P[path[i]], P[path[i + 1]]);
    return len;
}
bool twoOptOnce(const vector<Point>& P, vector<int>& path) {
    int n = path.size();
    if (n < 5) return false;                 // 안전 가드
    for (int i = 1; i < n - 2; ++i)
        for (int k = i + 1; k < n - 1; ++k) {
            double delta =
                dist(P[path[i - 1]], P[path[k]]) +
                dist(P[path[i]],     P[path[k + 1]]) -
                dist(P[path[i - 1]], P[path[i]]) -
                dist(P[path[k]],     P[path[k + 1]]);
            if (delta < -1e-9) {
                reverse(path.begin() + i, path.begin() + k + 1);
                return true;
            }
        }
    return false;
}
void twoOptFull(const vector<Point>& P, vector<int>& path, int maxIter = 6) {
    for (int iter = 0; iter < maxIter; ++iter)
        if (!twoOptOnce(P, path)) break;
}

/* ---------- 간단 Simulated Annealing ---------- */
vector<int> simulatedAnnealing(const vector<Point>& P, const vector<int>& seed,
                        double T0 = 800, double cool = 0.995,
                        int steps = 1500) {
    if (seed.size() <= 4) return seed;       // 도시 1~2개면 SA 생략

    default_random_engine rng(random_device{}());
    uniform_real_distribution<double> U(0, 1);

    vector<int> best = seed, cur = seed;
    double bestLen = tourLength(P, best), curLen = bestLen;

    auto swapTwo = [&](vector<int>& v) {
        int n = v.size();
        int a = 1 + rng() % (n - 2);
        int b = 1 + rng() % (n - 2);
        if (a > b) swap(a, b);
        reverse(v.begin() + a, v.begin() + b);
    };

    double T = T0;
    for (int s = 0; s < steps; ++s) {
        vector<int> nxt = cur;
        swapTwo(nxt);
        double nxtLen = tourLength(P, nxt);
        double dE = nxtLen - curLen;
        if (dE < 0 || U(rng) < exp(-dE / T)) {
            cur.swap(nxt);
            curLen = nxtLen;
            if (curLen < bestLen) { best = cur; bestLen = curLen; }
        }
        T *= cool;
    }
    return best;
}

/* ---------- MST 기반 2-approx (메타그래프 용) ---------- */

void dfsEuler(int u, const vector<vector<int>>& adj,
              vector<int>& path, vector<vector<int>>& used) {
    for (int v : adj[u]) if (!used[u][v]) {
        used[u][v] = used[v][u] = 1;
        dfsEuler(v, adj, path, used);
    }
    path.push_back(u);
}
vector<int> mst2Approx(const vector<Point>& P) {
    auto adj = buildMST(P);
    int n = P.size();
    vector<vector<int>> used(n, vector<int>(n));
    vector<int> eul; dfsEuler(0, adj, eul, used);
    reverse(eul.begin(), eul.end());
    vector<int> seen(n); vector<int> out;
    for (int v : eul) if (!seen[v]) { out.push_back(v); seen[v] = 1; }
    out.push_back(0);
    return out;
}

/* ---------- 격자 분할 ---------- */
vector<vector<int>> splitGrids(const vector<Point>& P) {
    double minX=1e100,minY=1e100,maxX=-1e100,maxY=-1e100;
    for (auto& p : P) { minX=min(minX,p.x); minY=min(minY,p.y);
                        maxX=max(maxX,p.x); maxY=max(maxY,p.y); }
    double w = (maxX - minX) / 4.0, h = (maxY - minY) / 4.0;
    vector<vector<int>> cell(16);
    for (int i = 0; i < (int)P.size(); ++i) {
        int cx = min(3, int((P[i].x - minX)/w));
        int cy = min(3, int((P[i].y - minY)/h));
        cell[cy*4 + cx].push_back(i);
    }
    return cell;
}

/* ---------------------------------------------------------------
   16-격자 Simulated-Annealing  →  메타 MST-2-approx → 글로벌 2-opt
   --------------------------------------------------------------- */
TSPResult gridSA_TSP(const vector<Point>& P)
{
    /* ───────────────── 0) 격자 분할 ───────────────── */
    const int G = 16;                                 // 4×4
    auto grid = splitGrids(P);                        // 각 격자별 도시 인덱스
    vector<vector<int>> sub(G);                              // 격자별 부분 투어

    /* ───────────────── 1) 격자 내부 SA ────────────── */
    for (int g = 0; g < G; ++g) {
        auto& idx = grid[g];
        if (idx.empty()) continue;

        if (idx.size() == 1) {                        // ★ (도시 1개) → 그대로
            sub[g] = { idx[0] };                      //   {v}
            continue;
        }

        vector<int> seed(idx.begin(), idx.end());
        seed.push_back(idx[0]);                       // 순환 시드
        sub[g] = simulatedAnnealing(P, seed);
    }

    /* ───────────────── 2) 메타 그래프(centroid) ──── */
    vector<Point> cent;   // centroid 좌표
    vector<int>   gid;    // cent[i] ← 어떤 격자?
    for (int g = 0; g < G; ++g)
        if (!sub[g].empty()) {
            double sx = 0, sy = 0;
            for (int v : grid[g]) { sx += P[v].x; sy += P[v].y; }
            cent.push_back({ sx / grid[g].size(), sy / grid[g].size() });
            gid.push_back(g);
        }

    /* ───────────────── 3) 메타 MST-2-approx 순서 ─── */
    vector<int> meta = mst2Approx(cent);                     // 0 … 0 (centroid 인덱스)

    /* ───────────────── 4) 투어 스티칭 ─────────────── */
    vector<int> total;
    for (int mi = 0; mi < (int)meta.size(); ++mi) {   // ★ 마지막 0까지 모두
        int g = gid[ meta[mi] ];                      // 실제 격자 번호
        const vector<int>& subP = sub[g];
        if (subP.empty()) continue;

        if (total.empty()) {                          // 첫 격자
            total = subP;
        } else {
            /* A. 연결점이 중복이면 하나 제거 */
            if (total.back() == subP.front())         // ★ pop 조건 수정
                total.pop_back();

            /* B. subP가 폐회로(첫==끝)인지 체크 */
            bool closed = (subP.size() >= 2 &&
                           subP.front() == subP.back());

            total.insert(total.end(),
                         subP.begin(),
                         closed ? subP.end() - 1      // ★ 마지막 노드 제외
                                : subP.end());
        }
    }

    /* ───────────────── 5) 순환 완성 & 2-opt ──────── */
    if (!total.empty() && total.front() != total.back())
        total.push_back(total.front());

    twoOptFull(P, total, 8);
    return {total, tourLength(P, total)};
}

/* ---------- 실행 ---------- */
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n; if (!(cin >> n)) return 0;
    vector<Point> P(n);
    for (auto& p : P) cin >> p.x >> p.y;

    vector<int> tour = gridSA_TSP(P).path;
    cout << fixed << setprecision(6)
         << "length = " << tourLength(P, tour) << '\n';
    for (int v : tour) cout << v << ' ';
    cout << '\n';
}