#include "ExperimentalBitSets.h"
#include "Constants.h"

#include <filesystem>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <bit>

ExperimentalBitSets::ExperimentalBitSets(const AppConfig& config, const SpectrumBitSet& sbs) : config_(config), sbs_(sbs)
{ 
}

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

void ExperimentalBitSets::loadSingleFile(const std::string& path_string)
{
    std::cout << "Loading file: " << '"' << path_string << '"' << std::endl;

    std::filesystem::path path(path_string);
    path = config_.output_path / path.filename();

    std::ifstream i(path_string);
    std::ofstream o(path.string());
    std::string buffer;

    std::vector<float> tmp_peaks;

    int count_files = 0;
    int count_filtered = 0;

    while (readEntryIntoBuffer(i, buffer))
    {
        count_files += 1;
        tmp_peaks = readPeaksFromBuffer(buffer);
        std::vector<uint64_t> tmp_bitset(int(((BIN_MAX_MZ - BIN_MIN_MZ) / config_.resolution) / 64) + 1);
        for (float p : tmp_peaks)
        {
            if (p < BIN_MIN_MZ || p > BIN_MAX_MZ) continue;
            int index = int((p - BIN_MIN_MZ) / config_.resolution);
            int word_index = index / 64;
            int bit_index = index % 64;
            tmp_bitset[word_index] |= (1ULL << bit_index);
        }
        if ((popCountBitsets(sbs_.bitset(), tmp_bitset, config_.cutoff)))
        {
            o << buffer;
            count_filtered += 1;
        }
    }
    std::cout << "  " << count_files << " spectras loaded; " << count_filtered << " spectras filtered; " << count_files - count_filtered << " spectras remaining" << std::endl;
    std::cout << "Remaining spectras written to " << path.string() << std::endl;
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
        buffer += line;
        buffer += "\n";
        if (line == "END IONS") break;
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
    while(getline(ss, value, ' ') && value != "END")
    {
        float peak = stof(value);
        tmp_peaks.push_back(peak);
        getline(ss, value, '\n');
    }
    return tmp_peaks;
}

bool ExperimentalBitSets::popCountBitsets(const std::vector<uint64_t>& lib, const std::vector<uint64_t>& exp, const double cutoff) const
{
    uint64_t exp_bit_count = 0;
    uint64_t intersection_bit_count = 0;

    for (uint i = 0; i < exp.size(); ++i)
    {
        uint64_t exp_bits = exp[i];
        uint64_t lib_bits = lib[i];

        exp_bit_count += std::popcount(exp_bits);
        intersection_bit_count += std::popcount(exp_bits & lib_bits);
    }

    double overlap_coefficient = 0.0;
    if (exp_bit_count != 0) overlap_coefficient = static_cast<double>(intersection_bit_count) / exp_bit_count;

    if (overlap_coefficient > cutoff) return true;
    return false;
}