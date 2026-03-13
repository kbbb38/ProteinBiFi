#include "ExperimentalSpectra.h"
#include "Constants.h"

#include <vector>
#include <string>
#include <sstream>
#include <bit>
#include <cmath>

#include <iostream>

ExperimentalSpectra::ExperimentalSpectra(const std::string& b, const AppConfig& config) : config_(config)
{
    std::string line, value;
    std::stringstream ss(b);

    getline(ss, line);
    getline(ss, line);

    value = line.substr(6, std::string::npos);
    name_ = value;

    while (getline(ss, line))
    {
        if (line.starts_with("PEPMASS="))
        {
            pepmass_ = stof(line.substr(8, std::string::npos));
        }
        if (line.starts_with("CHARGE="))
        {
            charge_ = stoi(line.substr(7,1));
            break;
        }
    }

    while(getline(ss, value, ' ') && value != "END")
    {
        peak_positions_.push_back(stof(value));
        getline(ss, value, '\n');
        intensities_.push_back(stof(value));
    }
    
    createBitSet();
    normalizeAndScaleIntensities();
    binIntensities();
}

void ExperimentalSpectra::createBitSet()
{
    bit_count_ = 0;
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

void ExperimentalSpectra::binIntensities()
{
    size_t num_bins = bitset_.size() * 64;
    binned_intensities_.resize(num_bins);

    std::vector<double> bin_squares(num_bins, 0.0);
    for (size_t i = 0; i < peak_positions_.size(); ++i) 
    {
        float intensity = intensities_[i];

        size_t bin_index = size_t((peak_positions_[i] - BIN_MIN_MZ) / config_.resolution);

        if (bin_index < num_bins)
        {
            bin_squares[bin_index] += (float(intensity) * intensity);
        }
    }

    for (size_t i = 0; i < num_bins; ++i) {
        binned_intensities_[i] = static_cast<float>(std::sqrt(bin_squares[i]));
    }
}

void ExperimentalSpectra::normalizeAndScaleIntensities()
{
    float magnitude = 0;
    for (float &i : intensities_)
    {
        i = std::sqrt(i);
        magnitude += i * i;
    }

    magnitude = std::sqrt(magnitude);

    for (float &i : intensities_)
    {
        i = i / magnitude;
    }
}