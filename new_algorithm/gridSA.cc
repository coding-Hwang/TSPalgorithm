#include <iostream>
#include <vector>
#include <algorithm>
#include <random>

#include "gridSA.h"


using namespace std;

/* ---------- 2-opt ---------- */

// path 길이 계산
double tourLength(const vector<Point>& P, const vector<int>& path) {
    double len = 0;
    for (int i = 0; i + 1 < (int)path.size(); ++i)
        len += dist(P[path[i]], P[path[i + 1]]);
    return len;
}

// 2-opt 한 번 수행 (full에서 호출)
bool twoOptOnce(const vector<Point>& P, vector<int>& path) {
    int n = path.size();
    if (n < 5) return false;  // 안전 가드
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

vector<int> simulatedAnnealing(const vector<Point>& P, const vector<int>& seed) {
    double T = 1000.0;       // 초기 온도
    double T_min = 1e-6;     // 최소 온도 (종료 조건)
    double cool = 0.995;

    if (seed.size() <= 4) return seed;

    srand(time(nullptr));  // rand() 초기화

    vector<int> best = seed, cur = seed;
    double bestLen = tourLength(P, best), curLen = bestLen;

    auto swapTwo = [&](vector<int>& v) {
        int n = v.size();

        int a = 1 + rand() % (n - 2);
        int b = 1 + rand() % (n - 2);
        if (a > b) swap(a, b);
        reverse(v.begin() + a, v.begin() + b);
    };

    int iter = 0;
    while (T > T_min) {
        for (int i = 0; i < 50; ++i) {  // 100번 시도 후 냉각
            vector<int> nxt = cur;
            swapTwo(nxt);
            double nxtLen = tourLength(P, nxt);
            double dE = nxtLen - curLen;

            if (dE < 0 || ((double)rand() / RAND_MAX) < exp(-dE / T)) {
                cur.swap(nxt);
                curLen = nxtLen;
                if (curLen < bestLen) {
                    best = cur;
                    bestLen = curLen;
                }
            }
        }
        T *= cool;
        iter++;
    }

    return best;
}


/* ---------- MST 기반 2-approx (메타그래프 용) ---------- */


vector<int> mst2Approx(const vector<Point>& P) {
    auto adj = buildMST(P);
    int n = P.size();
    vector<vector<bool>> used(n, vector<bool>(n, false));
    vector<int> eul; 
    
    dfs_eulerian(0, adj, used, eul);
    reverse(eul.begin(), eul.end());

    vector<int> seen(n); vector<int> out;

    for (int v : eul) if (!seen[v]) { out.push_back(v); seen[v] = 1; }
    out.push_back(0);
    return out;
}

/* ---------- 격자 분할 ---------- */
// 보조 함수: 재귀 분할
void subdivide(const vector<Point>& P,
               const vector<int>& idx,          // 이 영역에 들어있는 도시 인덱스
               double minX, double minY,
               double maxX, double maxY,
               vector<vector<int>>& cells,      // 완성된 셀들을 여기에 push
               const int LIMIT = 280)           // 한 셀에 허용되는 최대 점 개수
{
    if (idx.empty()) return;

    /* 1) 이미 조건 만족 → 그대로 셀 하나 추가 */
    if (static_cast<int>(idx.size()) <= LIMIT) {
        cells.push_back(idx);
        return;
    }

    /* 2) 4분할 좌표 계산 */
    double midX = (minX + maxX) * 0.5;
    double midY = (minY + maxY) * 0.5;

    vector<int> q[4];   // 0:(좌하) 1:(우하) 2:(좌상) 3:(우상)

    for (int v : idx) {
        bool right = P[v].x >= midX;
        bool top   = P[v].y >= midY;
        int  k     = (top ? 2 : 0) + (right ? 1 : 0);
        q[k].push_back(v);
    }

    /* 3) 각 사분면에 대해 재귀 호출 */
    subdivide(P, q[0], minX, minY, midX, midY, cells, LIMIT);
    subdivide(P, q[1], midX, minY, maxX, midY, cells, LIMIT);
    subdivide(P, q[2], minX, midY, midX, maxY, cells, LIMIT);
    subdivide(P, q[3], midX, midY, maxX, maxY, cells, LIMIT);
}

pair<vector<vector<int>>, int>
splitGrids(const vector<Point>& P)
{
    int LIMIT = 280;
    /* (a) 전체 영역 경계 구하기 */
    double minX =  numeric_limits<double>::max();
    double minY =  numeric_limits<double>::max();
    double maxX = -numeric_limits<double>::max();
    double maxY = -numeric_limits<double>::max();

    for (const auto& p : P) {
        minX = min(minX, p.x);  maxX = max(maxX, p.x);
        minY = min(minY, p.y);  maxY = max(maxY, p.y);
    }

    /* (b) 처음엔 전체 인덱스를 하나의 영역으로 */
    vector<int> rootIdx(P.size());
    iota(rootIdx.begin(), rootIdx.end(), 0);

    vector<vector<int>> cells;
    subdivide(P, rootIdx, minX, minY, maxX, maxY, cells, LIMIT);

    int G = static_cast<int>(cells.size());
    return {cells, G};
}
// 나누는 방식...
// 격자 크기를 데이터 크기 별로 한 격자에 16~280개 정도 들어가도록 조정


/* ---------------------------------------------------------------
   16-격자 Simulated-Annealing  →  메타 MST-2-approx → 글로벌 2-opt
   --------------------------------------------------------------- */
TSPResult gridSA_TSP(const vector<Point>& P)
{
    if (P.size() <= 280) return simulatedAnnealing(P);; // 도시 1~2개면 SA 생략

    /* ───────────────── 0) 격자 분할 ───────────────── */
    auto result = splitGrids(P); // 격자 분할
    auto grid = result.first; 
    int G = result.second;                   // 격자 개수
    cout << "Grid count: " << G << endl;

    vector<vector<int>> sub(G);                              // 격자별 부분 투어

    /* ───────────────── 1) 격자 내부 SA ────────────── */
    for (int g = 0; g < G; ++g) {
        auto& idx = grid[g];
        if (idx.empty()) continue;

        if (idx.size() == 1) {       // (도시 1개) → 그대로
            sub[g] = { idx[0] };     //   {v}
            continue;
        }

        vector<int> seed(idx.begin(), idx.end());
        seed.push_back(idx[0]);   // 순환 시드
        sub[g] = simulatedAnnealing(P, seed);
    }

    /* ───────────────── 2) 메타 그래프(centroid) ──── */
    vector<Point> cent;   // centroid 좌표 이용
    vector<int>   gid;    // cent[i] <- 어떤 격자?
    for (int g = 0; g < G; ++g)
        if (!sub[g].empty()) {
            double sx = 0, sy = 0;
            for (int v : grid[g]) { sx += P[v].x; sy += P[v].y; }
            cent.push_back({ sx / grid[g].size(), sy / grid[g].size() });
            gid.push_back(g);
        }
    // 각 grid 별 cent 찾기
    

    /* ───────────────── 3) 메타 MST-2-approx 순서 ─── */
    vector<int> meta = mst2Approx(cent); // 0 … 0 (centroid 인덱스)
    // 칸 방문 순서 정하기

    /* ───────────────── 4) 투어 스티칭 ─────────────── */
    vector<int> total;
    for (int mi = 0; mi < (int)meta.size(); ++mi) {   // 마지막 0까지 모두
        int g = gid[meta[mi]];                      // 실제 격자 번호
        const vector<int>& subP = sub[g];
        if (subP.empty()) continue;

        if (total.empty()) {   // 첫 격자
            total = subP;
        } else {
            /* A. 연결점이 중복이면 하나 제거 */
            if (total.back() == subP.front())
                total.pop_back();

            /* B. subP가 폐회로(첫==끝)인지 체크 */
            bool closed = (subP.size() >= 2 &&
                           subP.front() == subP.back());

            total.insert(total.end(),
                         subP.begin(),
                         closed ? subP.end() - 1      // 마지막 노드 제외
                                : subP.end());
        }
    }

    /* ───────────────── 5) 순환 완성 & 2-opt ──────── */
    if (!total.empty() && total.front() != total.back())
        total.push_back(total.front());

    twoOptFull(P, total, 8); // 8번 (default)
    return {total, tourLength(P, total)};
}