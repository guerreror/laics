#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include "json.hpp"
#include "migprob.h"

using json = nlohmann::json;

static bool load_json(const std::string& filename, json& out) {
    std::ifstream f(filename);
    if (!f) {
        std::cerr << "[migprob] ERROR: cannot open " << filename << "\n";
        return false;
    }
    try { f >> out; }
    catch (const std::exception& e) {
        std::cerr << "[migprob] ERROR: bad JSON: " << e.what() << "\n";
        return false;
    }
    return true;
}

Matrix readMigrationMatrix(const std::string& filename, const std::string& selected_time)
{
    json j;
    if (!load_json(filename, j)) return {};
    if (!j.contains(selected_time)) {
        std::cerr << "[migprob] ERROR: time key " << selected_time << " not found in " << filename << "\n";
        return {};
    }
    return j[selected_time].get<Matrix>();
}

std::vector<std::pair<double, Matrix>> readMigrationSchedule(const std::string& filename)
{
    json j;
    std::vector<std::pair<double, Matrix>> schedule;
    if (!load_json(filename, j)) return schedule;

    for (auto it = j.begin(); it != j.end(); ++it) {
        try {
            double t = std::stod(it.key());
            Matrix M = it.value().get<Matrix>();
            schedule.emplace_back(t, std::move(M));
        } catch (const std::exception& e) {
            std::cerr << "[migprob] WARNING: skipping key '" << it.key() << "': " << e.what() << "\n";
        }
    }
    std::sort(schedule.begin(), schedule.end(),
              [](auto& a, auto& b){ return a.first < b.first; });
    return schedule;
}
