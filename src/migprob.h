#ifndef MIGPROB_H
#define MIGPROB_H

#include <vector>
#include <string>
#include <map>

using Matrix = std::vector<std::vector<double>>;

Matrix readMigrationMatrix(const std::string& filename, const std::string& selected_time);

std::vector<std::pair<double, Matrix>> readMigrationSchedule(const std::string& filename);

#endif 
