#include "ExperimentalSpectra.h"
#include "Constants.h"

#include <vector>
#include <string>
#include <sstream>
#include <bit>

ExperimentalSpectra::ExperimentalSpectra(const std::string& b, const AppConfig& config) : config_(config)
{
    std::string line, value;
    std::stringstream ss(b);

    getline(ss, line);
    getline(ss, line);

    size_t colon_pos = line.find(':');
    value = line.substr(6, std::string::npos);
    name_ = value;

    while (getline(ss, line))
    {
        if (line.starts_with("CHARGE=")) break;
    }

    while(getline(ss, value, ' ') && value != "END")
    {
        peak_positions_.push_back(stof(value));
        getline(ss, value, '\n');
        intensities_.push_back(stof(value));
    }
    
    createBitSet();
}

void ExperimentalSpectra::createBitSet()
{
    bitset_.resize(int(((BIN_MAX_MZ - BIN_MIN_MZ) / config_.resolution) / 64) + 1);
    for (float p : peak_positions_)
    {
        if (p < BIN_MIN_MZ || p > BIN_MAX_MZ) continue;
        int index = int((p - BIN_MIN_MZ) / config_.resolution);
        int word_index = index / 64;
        int bit_index = index % 64;
        bitset_[word_index] |= (1ULL << bit_index);
    }
    for(uint64_t subset : bitset_)
    {
        bit_count_ += std::popcount(subset);
    }
}
