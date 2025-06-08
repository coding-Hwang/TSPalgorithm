#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <algorithm>

// algorithm
#include "../existing_TSP/__point.h"

#include "../existing_TSP/greedy_heuristic.h"
#include "../existing_TSP/Held-Karp.h"
#include "../existing_TSP/insertion_method.h"
#include "../existing_TSP/simulated_annealing.h"
#include "../existing_TSP/MST-based.h"

// #include "../new_algorithm/gridSA.h"
#include "../new_algorithm/gridSA_faster.h"

using namespace std;
using namespace chrono;

// TSP 파일 로드#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cctype>

vector<Point> load_tsp_file(const std::string& fname) {
    std::ifstream in(fname);
    if (!in) { std::cerr << "Cannot open " << fname << '\n'; return {}; }

    std::string line;
    bool reading = false;
    std::vector<Point> pts;

    auto trim  = [](std::string& s) {
        while(!s.empty() && (s.back()=='\r' || s.back()=='\n')) s.pop_back();
    };
    auto upper_nospace = [](std::string s) {
        s.erase(remove_if(s.begin(), s.end(), ::isspace), s.end());
        std::transform(s.begin(), s.end(), s.begin(), ::toupper);
        return s;
    };

    while (std::getline(in, line)) {
        trim(line);
        std::string tag = upper_nospace(line);

        if (tag == "NODE_COORD_SECTION") { reading = true; continue; }
        if (tag == "EOF") break;
        if (!reading) continue;                       // 아직 헤더 영역

        if (line.find_first_not_of(" \t") == std::string::npos) continue;  // 빈 줄

        std::istringstream iss(line);
        int    id;
        double x, y;
        if (iss >> id >> x >> y) {                    // 공백이 몇 개건 상관없음
            pts.push_back({x, y});
        }
    }

    if (pts.empty())
        std::cerr << "⚠️  Could NOT parse coordinates from " << fname << '\n';
    return pts;
}


// 파일 경로에서 csv 출력용 이름 얻기
string get_csv_filename(const string& tsp_path) {
    size_t slash_pos = tsp_path.find_last_of("/\\");
    string filename = tsp_path.substr(slash_pos + 1);  // a280.tsp
    size_t dot_pos = filename.find_last_of('.');
    string name_only = filename.substr(0, dot_pos);    // a280
    return name_only + "_result.csv";                  // a280_result.csv
}

// CSV 헤더 작성
void write_csv_header(const string& filename) {
    ofstream fout(filename, ios::out);
    fout << "TSP_File,Algorithm,Cost,Time_ms,Num_Cities\n";
    fout.close();
}

// CSV 행 추가
void append_csv_row(const string& filename, const string& tsp_file, const string& algo_name,
                    double cost, long long time_ms, int num_cities) {
    ofstream fout(filename, ios::app);
    fout << tsp_file << "," << algo_name << ","
         << fixed << setprecision(2) << cost << "," << time_ms << "," << num_cities << "\n";
    fout.close();
}

template <typename AlgoFunc>
void run_and_record(const string& tsp_file, const string& name, AlgoFunc algo,
                    const vector<Point>& points, const string& csv_file) {
    auto start = high_resolution_clock::now();
    TSPResult result = algo(points);
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start).count();

    int n = points.size();

    cout << "[" << name << "] "
         << "Cost: " << result.total_cost
         << ", Time: " << duration << " ms" << endl;

    append_csv_row(csv_file, tsp_file, name, result.total_cost, duration, n);
}


// 실험 실행
void run_exp(const string& tsp_file) {
    cout << "==== Running TSP algorithms on: " << tsp_file << " ====" << endl;
    vector<Point> points = load_tsp_file(tsp_file);
    cout << "Loaded points: " << points.size() << endl;

    int n = points.size();

    string csv_file = get_csv_filename(tsp_file);
    write_csv_header(csv_file);

    // run_and_record(tsp_file, "Grid SA", gridSA_TSP, points, csv_file); // O(n^2 log n) (heuristic)
    // run_and_record(tsp_file, "Grid SA fast", gridSA_fast, points, csv_file); // O(n^2 log n) (heuristic)

    // run_and_record(tsp_file, "Greedy Heuristic", greedy, points, csv_file); // O(n^2)
    // run_and_record(tsp_file, "Held-Karp", held_karp, points, csv_file); // O(n^2 * 2^n)
    run_and_record(tsp_file, "MST-based 2-approx", mst_2approx, points, csv_file); // O(n^2 log n)
    run_and_record(tsp_file, "Simulated Annealing", simulatedAnnealing, points, csv_file); // O(n^2 log n) (heuristic)
    // run_and_record(tsp_file, "Insertion Method", nearestInsertion, points, csv_file); // O(n^2)

}

// main
int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " [tsp_file_path]" << endl;
        return 1;
    }

    string tsp_file = argv[1];
    run_exp(tsp_file);
    return 0;
}
