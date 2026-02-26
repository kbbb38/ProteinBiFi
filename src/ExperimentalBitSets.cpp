#include "ExperimentalBitSets.h"
#include "Constants.h"

#include <filesystem>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <bit>

//tmp
#include <bitset>

ExperimentalBitSets::ExperimentalBitSets(const float resolution) : resolution_(resolution) { }

void ExperimentalBitSets::loadExperimentalBitSets(const std::string& path_string)
{
    const std::filesystem::path path(path_string);

    if (std::filesystem::exists(path))
    {
        if(std::filesystem::is_directory(path))
        {
            loadFromDirectory(path_string);
        }
        else if (std::filesystem::is_regular_file(path))
        {
            loadSingleFile(path_string);
        }
    }
}

void ExperimentalBitSets::loadFromDirectory(const std::string& path_string)
{
    std::filesystem::path path(path_string);
    size_t file_count = 0;

    for (const auto& entry : std::filesystem::directory_iterator(path))
    {
        if (entry.is_regular_file()) 
        {
            ++file_count;
        }
        
        
    }
    std::cout << "Directory detected! Loading " << file_count << " files..." << std::endl;

    for (const auto& entry : std::filesystem::directory_iterator(path)) 
    {
        if (entry.is_regular_file())
        {
            loadSingleFile(entry.path());
        }
    }
}

void ExperimentalBitSets::loadSingleFile(const std::string& path)
{
    std::cout << "Loading file: " << '"' << path << '"' << std::endl;

    std::ifstream f(path);
    std::string buffer;

    std::vector<float> tmp_peaks;
    while (readEntryIntoBuffer(f, buffer))
    {
        tmp_peaks = readPeaksFromBuffer(buffer);
        std::vector<uint64_t> tmp_bitset(int(((BIN_MAX_MZ - BIN_MIN_MZ) / resolution_) / 64) + 1);
        for (float p : tmp_peaks)
        {
            if (p < BIN_MIN_MZ || p > BIN_MAX_MZ) continue;
            int index = int((p - BIN_MIN_MZ) / resolution_);
            int word_index = index / 64;
            int bit_index = index % 64;
            tmp_bitset[word_index] |= (1ULL << bit_index);
        }
        experimental_bitsets_.push_back(std::move(tmp_bitset));
    }
    std::cout << "      " << experimental_bitsets_.size() << " spectras loaded" << std::endl;
}

bool ExperimentalBitSets::readEntryIntoBuffer(std::ifstream& f, std::string& buffer) const
{
    buffer.clear();

    std::string line;

    if (!getline(f, line) || line != "BEGIN IONS") return false;
    
    buffer += line;
    buffer += "\n";

    while(getline(f, line))
    {
        if (line == "END IONS") break;
        buffer += line;
        buffer += "\n";
    }
    return true;
}

std::vector<float> ExperimentalBitSets::readPeaksFromBuffer(const std::string& buffer) const
{
    std::string line, value;
    std::stringstream ss(buffer);

    while (getline(ss, line))
    {
        if (line.starts_with("CHARGE=")) break;
    }

    std::vector<float> tmp_peaks;
    while(getline(ss, value, ' '))
    {
        float peak = stof(value);
        tmp_peaks.push_back(peak);
        getline(ss, value, '\n');
    }
    return tmp_peaks;
}

void ExperimentalBitSets::filterExperimentalSpectra(const SpectrumBitSet& sbs)
{
    post_filter_indicies_.clear();
    const std::vector<uint64_t>& library_bitset = sbs.bitset(); 

    for(uint i = 0; i < experimental_bitsets_.size(); ++i)
    {
        if (popCountBitsets(library_bitset, experimental_bitsets_[i]))
        {
            post_filter_indicies_.push_back(i);
        }
    }
    std::cout << experimental_bitsets_.size() - post_filter_indicies_.size() << " out of " << experimental_bitsets_.size() << " spectras removed by bitset filter. " << std::endl;
    std::cout << post_filter_indicies_.size() << " spectras remaining" << std::endl;
}   

bool ExperimentalBitSets::popCountBitsets(const std::vector<uint64_t>& lib, const std::vector<uint64_t>& exp) const
{
    uint64_t exp_bit_count = 0;
    uint64_t intersection_bit_count = 0;

    for (uint i = 0; i < exp.size(); ++i)
    {
        uint64_t exp_bits = exp[i];
        uint64_t lib_bits = lib[i];

        if ((exp_bits & lib_bits) != exp_bits) return false;
    }
    return true;
}