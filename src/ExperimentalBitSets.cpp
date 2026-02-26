#include "ExperimentalBitSets.h"
#include "Constants.h"

#include <filesystem>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>

//tmp includes remove later
#include <bitset>

ExperimentalBitSets::ExperimentalBitSets(const float resolution) : resolution_(resolution) { }

void ExperimentalBitSets::loadExperimentalBitSets(const std::string& path_string)
{
    const std::filesystem::path path(path_string);

    if (std::filesystem::exists(path))
    {
        if(std::filesystem::is_directory(path))
        {
            
        }
        else if (std::filesystem::is_regular_file(path))
        {
            loadSingleFile(path_string);
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