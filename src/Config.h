#pragma once

#include <string>

struct AppConfig
{
    bool filter_experimental = false;
    std::string library_path;
    std::string experimental_path;
    std::string output_path;
    float resolution;
    double cutoff;
};