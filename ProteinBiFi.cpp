#include "src/SpectrumBitSet.h"
#include "src/Config.h"
#include "src/CutoffAnalysis.h"
#include "include/CLI11.hpp"

#include <iostream>
#include <chrono>
#include <string>
#include <filesystem>

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
    app.add_option("-g, --ground_truth", config.ground_truth_path, "Path to ground truth file for cutoff analysis");
    app.add_flag("-f, --filter_experimental", config.filter_experimental, "Filter experimental spectras instead of spectras in the seach library");

    app.parse(argc, argv);

    // Create output directory if it doesn't exist
    std::filesystem::create_directories(config.output_path);

    std::cout << "\033[1;32m"; 
    std::cout << "╔════════════════════╗" << std::endl;
    std::cout << "║ Protein Bit Filter ║" << std::endl;
    std::cout << "╚════════════════════╝" << std::endl;
    std::cout << "\033[0m";

    std::cout << "-> Loading experimental spectras..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    SpectrumBitSet sbs(config);
    sbs.loadFile(config.experimental_path);
    auto stop = std::chrono::high_resolution_clock::now();
    auto durration = duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "-> Done! Took " << durration.count() << " seconds! Loaded " << sbs.loaded() << " spectras!" <<std::endl;
    std::cout << " " << std::endl;

    std::cout << "-> Loading library spectras..." << std::endl;
    start = std::chrono::high_resolution_clock::now();
    sbs.loadFile(config.library_path);
    stop = std::chrono::high_resolution_clock::now();
    durration = duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "-> Done! Took " << durration.count() << " seconds! Loaded " << sbs.loaded() << " spectras!" <<std::endl;
    std::cout << " " << std::endl;

    float cutoff = 0.0f;
    if (!config.ground_truth_path.empty())
    {
        std::cout << "-> Starting cutoff analysis..." << std::endl;
        start = std::chrono::high_resolution_clock::now();
        CutoffAnalysis ca(sbs, config.ground_truth_path, config.output_path, config.resolution);
        stop = std::chrono::high_resolution_clock::now();
        durration = duration_cast<std::chrono::seconds>(stop - start);
        std::cout << "-> Done! Took " << durration.count() << " seconds!" << std::endl;
        std::cout << " " << std::endl;
        cutoff = ca.getMean();
    }

    std::cout << "-> Finding matches..." << std::endl;
    start = std::chrono::high_resolution_clock::now();
    sbs.matchSpectras();
    stop = std::chrono::high_resolution_clock::now();
    durration = duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "-> Done! Took " << durration.count() << " seconds!" << std::endl;
    std::cout << " " << std::endl;

    std::cout << "-> Filtering and writing output..." << std::endl;
    start = std::chrono::high_resolution_clock::now();
    sbs.writeOutput(config.output_path, config.experimental_path, cutoff);
    stop = std::chrono::high_resolution_clock::now();
    durration = duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "-> Done! Took " << durration.count() << " seconds! Filtered " << sbs.filtered() << " spectras!" << std::endl;
    std::cout << " " << std::endl;

    std::cout << "\033[1;32m"; 
    std::cout << "Run completed! Thank you for using ProteinBiFi!"<< std::endl;
    std::cout << "\033[0m";
    return 0;
}