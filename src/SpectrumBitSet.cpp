#include "SpectrumBitSet.h"
#include "Constants.h"

#include <string>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <bit>

namespace fs = std::filesystem;

SpectrumBitSet::SpectrumBitSet(AppConfig config) : config_(config)
{
}

void SpectrumBitSet::loadFile(const std::string& path_string)
{
    const fs::path path(path_string);

    if (fs::exists(path))
    {
        if(fs::is_directory(path))
        {
            loadFromDirectory(path_string);
        }
        else if (std::filesystem::is_regular_file(path))
        {
            loadSingleFile(path_string);
        }
    }

    bitset_complete_ = true;
}

void SpectrumBitSet::loadFromDirectory(const std::string& path_string)
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

void SpectrumBitSet::loadSingleFile(const std::string& path_string)
{
    std::cout << "Loading file: " << '"' << path_string << '"' << std::endl;
    std::ifstream f(path_string);
    std::string buffer;

    std::vector<float> tmp_peaks;
    if (!bitset_complete_)
    {
        bitset_.resize(int(((BIN_MAX_MZ - BIN_MIN_MZ) / config_.resolution) / 64) + 1);
        int load_count = 0;

        if (!config_.filter_experimental)
        {
            while (readEntryIntoBufferExp(f, buffer))
            {
                tmp_peaks = readPeaksFromBufferExp(buffer);
                addPeaksToBitset(tmp_peaks);
                load_count += 1;
            }
            std::cout << "  Loaded " << load_count << " spectras" << std::endl;
        }
        else
        {
            while (readEntryIntoBufferLib(f, buffer))
            {
                tmp_peaks = readPeaksFromBufferLib(buffer);
                addPeaksToBitset(tmp_peaks);
                load_count += 1;
            }
            std::cout << "  Loaded " << load_count << " spectras" << std::endl;
        }
        total_loaded_ += load_count;
    }
    else
    {
        std::filesystem::path path(path_string);
        path = config_.output_path / path.filename();
        std::ofstream o(path.string());

        int file_count = 0;
        int nfiltered_count = 0;

        std::vector<uint64_t> tmp_bitset(int(((BIN_MAX_MZ - BIN_MIN_MZ) / config_.resolution) / 64) + 1);
        if (!config_.filter_experimental)
        {
            while (readEntryIntoBufferLib(f, buffer))
            {
                file_count += 1;
                tmp_peaks = readPeaksFromBufferLib(buffer);
                tmp_bitset = readPeaksIntoBitset(tmp_peaks);

                if(popCountBitsets(tmp_bitset))
                {
                    o << buffer;
                    nfiltered_count += 1;
                }
            }
            total_filtered_ += file_count - nfiltered_count;
            std::cout << "  " << file_count << " spectras loaded; " << file_count - nfiltered_count << " spectras filtered; " << nfiltered_count << " spectras remaining" << std::endl;
            std::cout << "  Remaining spectras written to: " << path.string() << std::endl;
        }
        else
        {
            while (readEntryIntoBufferExp(f, buffer))
            {
                file_count += 1;
                tmp_peaks = readPeaksFromBufferExp(buffer);
                tmp_bitset = readPeaksIntoBitset(tmp_peaks);

                if(popCountBitsets(tmp_bitset))
                {
                    o << buffer;
                    nfiltered_count += 1;
                }
            }
            total_filtered_ += file_count - nfiltered_count;
            std::cout << "  " << file_count << " spectras loaded; " << file_count - nfiltered_count << " spectras filtered; " << nfiltered_count << " spectras remaining" << std::endl;
            std::cout << "  Remaining spectras written to: " << path.string() << std::endl;
        }
    }
}

bool SpectrumBitSet::readEntryIntoBufferExp(std::ifstream& f, std::string& buffer) const
{
    buffer.clear();

    std::string line;

    if (!getline(f, line) || line != "BEGIN IONS") return false;
    
    buffer += line;
    buffer += "\n";

    while(getline(f, line))
    {
        buffer += line;
        buffer += "\n";
        if (line == "END IONS") break;
    }
    return true;
}

std::vector<float> SpectrumBitSet::readPeaksFromBufferExp(const std::string& buffer) const
{
    std::string line, value;
    std::stringstream ss(buffer);

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
    return tmp_peaks;    
}

bool SpectrumBitSet::readEntryIntoBufferLib(std::ifstream& f, std::string& buffer) const
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

std::vector<float> SpectrumBitSet::readPeaksFromBufferLib(const std::string& buffer) const
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

void SpectrumBitSet::addPeaksToBitset(const std::vector<float>& tmp_peaks)
{
    for (float p : tmp_peaks)
    {
        if (p < BIN_MIN_MZ || p > BIN_MAX_MZ) continue;
        int index = int((p - BIN_MIN_MZ) / config_.resolution);
        int word_index = index / 64;
        int bit_index = index % 64;
        bitset_[word_index] |= (1ULL << bit_index);
    }
}

std::vector<uint64_t> SpectrumBitSet::readPeaksIntoBitset(const std::vector<float>& tmp_peaks) const
{
    std::vector<uint64_t> tmp_bitset(int(((BIN_MAX_MZ - BIN_MIN_MZ) / config_.resolution) / 64) + 1);
    for (float p : tmp_peaks)
    {
        if (p < BIN_MIN_MZ || p > BIN_MAX_MZ) continue;
        int index = int((p - BIN_MIN_MZ) / config_.resolution);
        int word_index = index / 64;
        int bit_index = index % 64;
        tmp_bitset[word_index] |= (1ULL << bit_index);
    }
    return tmp_bitset;
}

bool SpectrumBitSet::popCountBitsets(const std::vector<uint64_t>& tmp_bitset) const
{
    uint64_t exp_bit_count = 0;
    uint64_t intersection_bit_count = 0;

    for (uint i = 0; i < tmp_bitset.size(); ++i)
    {
        uint64_t exp_bits = tmp_bitset[i];
        uint64_t lib_bits = bitset_[i];

        exp_bit_count += std::popcount(exp_bits);
        intersection_bit_count += std::popcount(exp_bits & lib_bits);
    }

    double overlap_coefficient = 0.0;
    std::cout << overlap_coefficient << " ";
    if (exp_bit_count != 0) overlap_coefficient = static_cast<double>(intersection_bit_count) / exp_bit_count;

    if (overlap_coefficient > config_.cutoff) return true;
    return false;
}