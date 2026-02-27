#pragma once

#include <string>

struct AppConfig
{
    std::string library_path;
    std::string experimental_path;
    std::string output_path;
    float resolution;
    double cutoff;
};