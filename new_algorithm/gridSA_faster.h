#pragma once
#include <vector>
#include "../existing_TSP/MST-based.h"
#include "../existing_TSP/__point.h"
#include "../existing_TSP/simulated_annealing.h"

TSPResult gridSA_fast(const std::vector<Point>& points);