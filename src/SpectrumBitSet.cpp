#include "SpectrumBitSet.h"
#include "Constants.h"

#include <fstream>
#include <iostream>
#include <sstream>

SpectrumBitSet::SpectrumBitSet(const float resolution) : resolution_(resolution)
{
    const uint32_t vector_size = int(((BIN_MAX_MZ - BIN_MIN_MZ) / resolution_) / 64) + 1;
    bitset_.resize(vector_size);
}

void SpectrumBitSet::loadFromFile(const std::string& path)
{
    std::ifstream f(path);
    std::string buffer;

    std::vector<float> tmp_peaks;
    while (readEntryIntoBuffer(f, buffer))
    {
        tmp_peaks = readPeaksFromBuffer(buffer);
        
        for (float p : tmp_peaks)
        {
            if (p < BIN_MIN_MZ || p > BIN_MAX_MZ) continue;
            int index = int((p - BIN_MIN_MZ) / resolution_);
            int word_index = index / 64;
            int bit_index = index % 64;
            bitset_[word_index] |= (1ULL << bit_index);
        }
    }

    for (auto i : bitset_)
    {
        std::cout << i << std::endl;
    }

    f.close();
}

bool SpectrumBitSet::readEntryIntoBuffer(std::ifstream& f, std::string& buffer) const
{
    buffer.clear();

    std::string line;
    
    if (!getline(f, line)) return false;
    if (line.compare(0, 5, "Name:") != 0) return false;

    buffer += line;
    buffer += "\n";

    while (getline(f, line))
    {
        if (line.compare(0, 5, "Name:") == 0)
        {
            f.clear();
            f.seekg(-(static_cast<std::streamoff>(line.size()+1)), std::ios::cur);
            break;
        }
        buffer += line;
        buffer += "\n";
    }
    return true;
}

std::vector<float> SpectrumBitSet::readPeaksFromBuffer(const std::string& buffer) const
{
    std::string line, value;
    std::stringstream ss(buffer);

    while (getline(ss, line))
    {
        if (line.starts_with("Num peaks:")) break;
    }
    
    size_t colon_pos = line.find(':');
    value = line.substr(colon_pos + 2, std::string::npos);
    uint16_t num_peaks = stoi(value);

    std::vector<float> tmp_peaks;
    tmp_peaks.reserve(num_peaks);

    for (uint16_t i = 0; i < num_peaks; ++i)
    {
        getline(ss, value, '\t');
        float peak = stof(value);
        tmp_peaks.push_back(peak);
        getline(ss, value, '\n');
    }
    return tmp_peaks;
}