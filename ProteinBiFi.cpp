#include "src/SpectrumBitSet.h"
#include "src/Config.h"
#include "src/CutoffAnalysis.h"
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

    // Options 
    app.add_option("-l, --library", config.library_path, "Path to spectral library file or directory")->required();
    app.add_option("-e, --experimental", config.experimental_path, "Path to experimental file or directory")->required();
    app.add_option("-o, --out", config.output_path, "Output directory")->required();
    app.add_option("-r, --resolution", config.resolution, "Resolution fo the bit sets")->required();
    app.add_option("-g, --ground-truth", config.ground_truth_path, "Path to ground truth file")->required();

    app.parse(argc, argv);

    std::cout << "\033[1;32m"; 
    std::cout << "╔════════════════════╗" << std::endl;
    std::cout << "║ Protein Bit Filter ║" << std::endl;
    std::cout << "╚════════════════════╝" << std::endl;
    std::cout << "\033[0m";

    // Load experimental file (mgf)
    std::cout << "-> Loading experimental spectras..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    SpectrumBitSet sbs(config);
    sbs.loadFile(config.experimental_path);
    auto stop = std::chrono::high_resolution_clock::now();
    auto durration = duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "-> Done! Took " << durration.count() << " seconds! Loaded " << sbs.loaded() << " spectras!" <<std::endl;
    std::cout << " " << std::endl;

    // Load library file (msp)
    std::cout << "-> Loading library spectras..." << std::endl;
    start = std::chrono::high_resolution_clock::now();
    sbs.loadFile(config.library_path);
    stop = std::chrono::high_resolution_clock::now();
    durration = duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "-> Done! Took " << durration.count() << " seconds! Loaded " << sbs.loaded() << " spectras!" <<std::endl;
    std::cout << " " << std::endl;

    // Analyse the cutoff of ground truth matches
    std::cout << "-> Starting cutoff analysis..." << std::endl;
    start = std::chrono::high_resolution_clock::now();
    CutoffAnalysis ca(sbs, config.ground_truth_path, config.output_path);
    stop = std::chrono::high_resolution_clock::now();
    durration = duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "-> Done! Took " << durration.count() << " seconds!" << std::endl;
    std::cout << " " << std::endl;

    // Find best tanimoto score for all library spectras
    std::cout << "-> Finding matches..." << std::endl;
    start = std::chrono::high_resolution_clock::now();
    sbs.matchSpectras();
    stop = std::chrono::high_resolution_clock::now();
    durration = duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "-> Done! Took " << durration.count() << " seconds!" << std::endl;
    std::cout << " " << std::endl;

    // Write tanimoto scores and spectras above cutoff into corresponding files
    std::cout << "-> Filtering and writing output..." << std::endl;
    start = std::chrono::high_resolution_clock::now();
    sbs.writeOutput(config.output_path);
    stop = std::chrono::high_resolution_clock::now();
    durration = duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "-> Done! Took " << durration.count() << " seconds! Filtered " << sbs.filtered() << " spectras!" << std::endl;
    std::cout << " " << std::endl;

    // Programm is done!
    std::cout << "\033[1;32m"; 
    std::cout << "Run completed! Thank you for using ProteinBiFi!"<< std::endl;
    std::cout << "\033[0m";
    return 0;
}