#include "src/SpectrumBitSet.h"
#include "src/Config.h"
#include "include/CLI11.hpp"

#include <iostream>
#include <chrono>
#include <string>

int main(int argc, char** argv)
{
    /*
        Option handling
    */
    CLI::App app{"Protein Bit Filter"};
    argv = app.ensure_utf8(argv);

    AppConfig config;

    app.add_option("-l, --library", config.library_path, "Path to spectral library file")->required();
    app.add_option("-e, --experimental", config.experimental_path, "Path to experimental file or directory")->required();
    app.add_option("-o, --out", config.output_path, "Output directory")->required();
    app.add_option("-r, --resolution", config.resolution, "Resolution fo the bit sets")->required();
    app.add_option("-c, --cutoff", config.cutoff, "Cutoff for overlap coefficient during filtering")->required();
    app.add_flag("-f, --filter_experimental", config.filter_experimental, "Filter experimental spectras instead of spectras in the seach library");

    app.parse(argc, argv);

    std::cout << "\033[1;32m"; 
    std::cout << "╔════════════════════╗" << std::endl;
    std::cout << "║ Protein Bit Filter ║" << std::endl;
    std::cout << "╚════════════════════╝" << std::endl;
    std::cout << "\033[0m";

    std::cout << "-> Constructing spectrum bit set..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    SpectrumBitSet sbs(config);
    config.filter_experimental ? sbs.loadFile(config.library_path) : sbs.loadFile(config.experimental_path);
    auto stop = std::chrono::high_resolution_clock::now();
    auto durration = duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "-> Done! Took " << durration.count() << " seconds! Loaded " << sbs.loaded() << " spectras!" <<std::endl;
    std::cout << " " << std::endl;

    std::cout << "-> Loading and filtering spectras..." << std::endl;
    start = std::chrono::high_resolution_clock::now();
    config.filter_experimental ? sbs.loadFile(config.experimental_path) : sbs.loadFile(config.library_path);
    stop = std::chrono::high_resolution_clock::now();
    durration = duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "-> Done! Took " << durration.count() << " seconds! Filtered" << sbs.filtered() << " spectras!" << std::endl;
    std::cout << " " << std::endl;

    std::cout << "\033[1;32m"; 
    std::cout << "Run completed! Thank you for using ProteinBiFi!"<< std::endl;
    std::cout << "\033[0m";
    return 0;
}