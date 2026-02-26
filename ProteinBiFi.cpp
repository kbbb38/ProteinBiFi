#include "src/SpectrumBitSet.h"
#include "src/ExperimentalBitSets.h"

#include <iostream>
#include <thread>
#include <chrono>

int main()
{
    std::cout << "\033[1;32m"; 
    std::cout << "╔════════════════════╗" << std::endl;
    std::cout << "║ Protein Bit Filter ║" << std::endl;
    std::cout << "╚════════════════════╝" << std::endl;
    std::cout << "\033[0m";

    std::cout << "-> Constructing spectrum bit set..." << std::endl;
    SpectrumBitSet sbs(0.1);
    sbs.loadBitSet("/home/malekk/scratch/MistleRun/data/human_consensus_final_true_lib.msp");
    std::cout << "Done!" << std::endl;;
    std::cout << "═════════════════════" << std::endl;

    std::cout << "-> Loading experimental spectras..." << std::endl;
    ExperimentalBitSets ebs(0.1);
    ebs.loadExperimentalBitSets("/home/malekk/scratch/MistleRun/data/PXD001197_MFG/OrbiElite04549.mgf");
    std::cout << "Done!" << std::endl;
    std::cout << "═════════════════════" << std::endl;

    std::cout << "-> Filtering spectras..." << std::endl;
    std::cout << "Done!" << std::endl;
    std::cout << "═════════════════════" << std::endl;
    return 0;
}