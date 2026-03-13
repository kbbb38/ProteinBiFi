#include "LibrarySpectra.h"
#include "Constants.h"

#include <string>
#include <sstream>
#include <bit>
#include <cmath>

#include <iostream>

LibrarySpectra::LibrarySpectra(const std::string& b, const AppConfig& config) : config_(config)
{
    std::string line, value;
    std::stringstream ss(b);

    getline(ss, line);

    size_t slash_pos = line.find('/');

    peptide_ = line.substr(6, slash_pos - 6);
    charge_ = stoi(line.substr(slash_pos + 1, 1));
    
    if (!peptide_.empty() && peptide_.back() == '\r') 
    {
        peptide_.pop_back();  
    }

    getline(ss, line);

    pepmass_ = stof(line.substr(3, std::string::npos));

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
    normalizeAndScaleIntensities();
    binIntensities();
}

void LibrarySpectra::createBitSet()
{
    bitset_.resize(int(((BIN_MAX_MZ - BIN_MIN_MZ) / config_.resolution) / 64) + 1);
    bit_count_ = 0;
    
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

void LibrarySpectra::binIntensities()
{
    size_t num_bins = bitset_.size() * 64 + 1;
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

void LibrarySpectra::normalizeAndScaleIntensities()
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