#pragma once
#include <vector>
#include "__point.h"

vector<vector<int>> buildMST(const vector<Point>& points);
void dfs_eulerian(int u,
                   const vector<vector<int>>& adj,
                   vector<vector<bool>>& used,
                   vector<int>& path);
vector<int> shortcut (int n, vector<int> path);
double calculate_total_distance(const vector<int>& path, const vector<Point>& points);
TSPResult mst_2approx(const vector<Point>& points);
