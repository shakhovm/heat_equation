//
// Created by shakhov on 4/4/20.
//

#include "../../include/conf_handler/confhandler.h"
#include <sstream>
#include <regex>
ConfHandler::ConfHandler(const std::string& filename) {
    std::ifstream conf_file(filename);
    if (!conf_file.is_open()) {
        throw std::runtime_error("File " + filename +  " is not found!");
    }

    content = dynamic_cast<std::ostringstream&>(
            std::ostringstream{} << conf_file.rdbuf()).str();
    conf_file_handler();
}

std::string ConfHandler::file_pattern(const std::string& pattern) {
    std::smatch match;
    std::regex r(pattern);

    regex_search(content, match, r);
    if (match.empty()) {
        throw std::runtime_error("Cannot handle conf file!");
    }
    std::string new_string = match.str();;
    return new_string.substr(new_string.find_first_of('=') + 1);
}

void ConfHandler::conf_file_handler() {

    conf_params.matrix_file = file_pattern(R"(matrix_file=.+\..+)");
    conf_params.out_file = file_pattern(R"(out_file=.+\..+)");
    std::stringstream(file_pattern( R"(delta_t=[\d\.]+)")) >> conf_params.delta_t;
    std::stringstream(file_pattern( R"(delta_x=[\d\.]+)")) >> conf_params.delta_x;
    std::stringstream(file_pattern( R"(delta_y=[\d\.]+)")) >> conf_params.delta_y;
    std::stringstream(file_pattern( R"(density=[\d\.]+)")) >> conf_params.density;
    std::stringstream(file_pattern( R"(temp_capacity=[\d\.]+)")) >> conf_params.temp_capacity;
    std::stringstream(file_pattern( R"(temp_conduct=[\d\.]+)")) >> conf_params.temp_conduct;
    std::stringstream(file_pattern( R"(time_to_save=[\d]+)")) >> conf_params.time_to_save;
    std::stringstream(file_pattern( R"(max_time=[\d]+)")) >> conf_params.max_time;
}
