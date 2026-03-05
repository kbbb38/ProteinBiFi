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

    std::vector<float> tmp_peaks;
    while(getline(ss, value, ' ') && value != "END")
    {
        float peak = stof(value);
        tmp_peaks.push_back(peak);
        getline(ss, value, '\n');
    }
    
    createBitSet(tmp_peaks);
}

void ExperimentalSpectra::createBitSet(const std::vector<float>& tmp_peaks)
{
    bitset_.resize(int(((BIN_MAX_MZ - BIN_MIN_MZ) / config_.resolution) / 64) + 1);
    for (float p : tmp_peaks)
    {
        if (p < BIN_MIN_MZ || p > BIN_MAX_MZ) continue;
        int index = int((p - BIN_MIN_MZ) / config_.resolution);
        int word_index = index / 64;
        int bit_index = index % 64;
        bitset_[word_index] |= (1ULL << bit_index);
    }
    for(size_t i = 0; i < bitset_.size(); ++i)
    {
        bit_count_ += std::popcount(bitset_[i]);
    }
}
