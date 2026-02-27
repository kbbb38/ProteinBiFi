#include "src/SpectrumBitSet.h"
#include "src/ExperimentalBitSets.h"
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

    app.parse(argc, argv);

    std::cout << "\033[1;32m"; 
    std::cout << "╔════════════════════╗" << std::endl;
    std::cout << "║ Protein Bit Filter ║" << std::endl;
    std::cout << "╚════════════════════╝" << std::endl;
    std::cout << "\033[0m";

    std::cout << "-> Constructing spectrum bit set..." << std::endl;
    auto start_lib = std::chrono::high_resolution_clock::now();
    SpectrumBitSet sbs(config.resolution);
    sbs.loadBitSet(config.library_path);
    auto stop_lib = std::chrono::high_resolution_clock::now();
    auto durration_lib = duration_cast<std::chrono::seconds>(stop_lib - start_lib);
    std::cout << "...Done! Took " << durration_lib.count() << " seconds" << std::endl;;
    std::cout << "═════════════════════" << std::endl;

    std::cout << "-> Loading and filtering experimental spectras..." << std::endl;
    start_lib = std::chrono::high_resolution_clock::now();
    ExperimentalBitSets ebs(config, sbs);
    ebs.loadExperimentalBitSets(config.experimental_path);
    stop_lib = std::chrono::high_resolution_clock::now();
    durration_lib = duration_cast<std::chrono::seconds>(stop_lib - start_lib);
    std::cout << "...Done! Took " << durration_lib.count() << " seconds" << std::endl;;
    std::cout << "═════════════════════" << std::endl;

    std::cout << "\033[1;32m"; 
    std::cout << "Run completed, thank you for using ProteinBiFi!" << std::endl;
    std::cout << "\033[0m";
    return 0;
}