#ifndef CONFHANDLER_H
#define CONFHANDLER_H

#include <fstream>
class ConfHandler {
public:
    struct ConfParams {
        std::string matrix_file;
        std::string out_file;
        double delta_t;
        double delta_x;
        double delta_y;
        double density;
        double temp_capacity;
        double temp_conduct;
        size_t time_to_save;
        size_t max_time;
    };
private:
    std::string content;

    ConfParams conf_params;

public:
    ConfHandler() = default;
    explicit ConfHandler(const std::string& filename);
    ~ConfHandler() = default;
    std::string file_pattern(const std::string& pattern);
    void conf_file_handler();
    inline ConfParams& getConfParams() { return conf_params; }
};
#endif // CONFHANDLER_H
