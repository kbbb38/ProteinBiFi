#include "LibrarySpectra.h"
#include "Constants.h"

#include <string>
#include <sstream>
#include <bit>

LibrarySpectra::LibrarySpectra(const std::string& b, const AppConfig& config) : config_(config)
{
    std::string line, value;
    std::stringstream ss(b);

    getline(ss, line);
    peptide_ = line.substr(6, std::string::npos);

    while (getline(ss, line))
    {
        if (line.starts_with("Num peaks:")) break;
    }

    size_t colon_pos = line.find(':');
    value = line.substr(colon_pos + 2, std::string::npos);
    uint16_t num_peaks = stoi(value);

    std::vector<float> tmp_peaks;

    for (uint16_t i = 0; i < num_peaks; ++i)
    {
        getline(ss, value, '\t');
        peak_positions_.push_back(stof(value));
        getline(ss, value, '\n');
        intensities_.push_back(stof(value));
    }

    createBitSet();
}

void LibrarySpectra::createBitSet()
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