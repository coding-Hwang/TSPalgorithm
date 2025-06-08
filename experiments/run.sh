#!/bin/bash

# 🔧 1. 컴파일
echo "==> Compiling experiment.cc and all algorithm files..."
g++ -O2 -std=c++17 -o experiment \
    exp.cc \
    ../new_algorithm/me.cc \
    ../existing_TSP/greedy_heuristic.cc \
    ../existing_TSP/Held-Karp.cc \
    ../existing_TSP/insertion_method.cc \
    ../existing_TSP/MST-based.cc \
    ../existing_TSP/simulated_annealing.cc

# 📂 2. 테스트할 TSP 파일 목록
TSP_LIST=(
    "../dataset/UNIstar.tsp"
    "../dataset/a280.tsp"
    "../dataset/xql662.tsp"
    "../dataset/kz9976.tsp"
    "../dataset/mona-lisa100K.tsp"
)

# 🧪 3. 각 TSP 파일에 대해 실험 실행
for tsp_file in "${TSP_LIST[@]}"; do
    echo ""
    echo "----------------------------------------"
    echo "🚀Running experiment on $tsp_file"
    if ./experiment "$tsp_file"; then
    echo "✅ Result saved to: $(basename "$tsp_file" .tsp)_result.csv"
    else
        echo "❌ Experiment failed for: $tsp_file"
    fi

done
