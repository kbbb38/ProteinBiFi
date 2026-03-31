#include "SpectrumBitSet.h"
#include "Constants.h"

#include <string>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <bit>
#include <algorithm>

namespace fs = std::filesystem;

SpectrumBitSet::SpectrumBitSet(AppConfig config) : config_(config)
{
}

void SpectrumBitSet::loadFile(const std::string& path_string)
{
    total_loaded_ = 0;
    is_dir_ = false;
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
}

void SpectrumBitSet::loadFromDirectory(const std::string& path_string)
{
    is_dir_ = true;
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
    if (path_string.ends_with(".mgf"))
    {
        std::cout << "Loading file: " << '"' << path_string << '"' << std::endl;
        std::ifstream f(path_string);
        std::string buffer;

        int load_count = 0;

        while (readEntryIntoBufferMgf(f, buffer))
        {
            experimental_spectra_.emplace_back(buffer, config_);
            load_count += 1;
        }
        total_loaded_ += load_count;
        if (is_dir_) std::cout << "      +" << load_count << " spectra" << std::endl;
        f.close();
    }
    if (path_string.ends_with(".msp"))
    {
        std::cout << "Loading file: " << '"' << path_string << '"' << std::endl;
        std::ifstream f(path_string);
        std::string buffer;

        int load_count = 0;

        while (readEntryIntoBufferMsp(f, buffer))
        {
            library_spectra_.emplace_back(buffer, config_);
            load_count += 1;
        }
        total_loaded_ += load_count;
        if (is_dir_) std::cout << "      +" << load_count << " spectra" << std::endl;
        f.close();
    }
}

bool SpectrumBitSet::readEntryIntoBufferMgf(std::ifstream& f, std::string& buffer) const
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

bool SpectrumBitSet::readEntryIntoBufferMsp(std::ifstream& f, std::string& buffer) const
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

void SpectrumBitSet::matchSpectras()
{
    for (ExperimentalSpectra& e_spec : experimental_spectra_)
    {
        size_t index_count = 0;

        float highest_t = 0;
        int highest_id_t = 0;

        for (LibrarySpectra& l_spec : library_spectra_)
        {
            if(l_spec.getCharge() == e_spec.getCharge() && abs(l_spec.getMass() - e_spec.getMass() < 10))
            {
                float score_t = calculateTanimotoScore(e_spec.getBitset(), e_spec.getBitCount(), l_spec.getBitset(), l_spec.getBitCount());

                if (score_t > highest_t)
                {
                    highest_t = score_t;
                    highest_id_t = index_count;
                }
            }
            ++index_count;
        }
        library_spectra_[highest_id_t].setIfMatch();
        e_spec.setMatch(Match(highest_t, library_spectra_[highest_id_t].getPeptide(), highest_id_t, true));
    }
}

float SpectrumBitSet::calculateTanimotoScore(const std::vector<uint64_t>& e_spec, const uint64_t e_count, const std::vector<uint64_t>& l_spec, const uint64_t l_count) const
{
    uint64_t count_intersection = 0;
    const uint64_t* e_ptr = e_spec.data();
    const uint64_t* l_ptr = l_spec.data();
    const size_t size = e_spec.size();

    for(size_t i = 0; i < size; ++i)
    {
        count_intersection += std::popcount(e_ptr[i] & l_ptr[i]);
    }

    if (e_count == 0 && l_count == 0) return 0.0;
    uint64_t count_union = e_count + l_count - count_intersection;
    if (count_union == 0) return 0.0;

    return float(count_intersection) / float(count_union);
}

float SpectrumBitSet::calculateOverlapCoefficient(const std::vector<uint64_t>& e_spec, const uint64_t e_count, const std::vector<uint64_t>& l_spec, const uint64_t l_count) const
{
    uint64_t count_intersection = 0;
    const uint64_t* e_ptr = e_spec.data();
    const uint64_t* l_ptr = l_spec.data();
    const size_t size = e_spec.size();

    for(size_t i = 0; i < size; ++i)
    {
        count_intersection += std::popcount(e_ptr[i] & l_ptr[i]);
    }

    if (e_count == 0 || l_count == 0) return 0.0;
    uint64_t min_count = std::min(e_count, l_count);
    return float(count_intersection) / float(min_count);
}

float SpectrumBitSet::calculateDotProductScore(const std::vector<uint64_t>& e_spec, const std::vector<float>& e_intensities, const std::vector<uint64_t>& l_spec, const std::vector<float>& l_intensities) const
{
    float dot = 0.0;
    for (size_t i = 0; i < e_spec.size(); ++i)
    {
        uint64_t overlap = e_spec[i] & l_spec[i];

        while (overlap != 0)
        {
            size_t bit_pos = std::countr_zero(overlap);
            size_t bin_index = (i * 64) + bit_pos;

            dot += e_intensities[bin_index] * l_intensities[bin_index];

            overlap &= (overlap - 1);
        }
    }
    return dot;
}

void SpectrumBitSet::writeOutput(const std::string& path_string_out, std::string& path_string_in, const float cutoff)
{
    std::ifstream f(path_string_in);
    
    // Create output filename based on resolution and cutoff value
    std::string bin_size_str = std::to_string(config_.resolution);
    // Remove trailing zeros from bin_size string
    bin_size_str.erase(bin_size_str.find_last_not_of('0') + 1, std::string::npos);
    if (bin_size_str.back() == '.') bin_size_str.pop_back();
    
    std::string filtered_filename = "filtered_spectra_bin_" + bin_size_str + "_cutoff_" + std::to_string(cutoff) + ".mgf";
    std::string all_scores_filename = "all_tanimoto_scores_bin_" + bin_size_str + ".txt";
    std::string not_gt_scores_filename = "not_ground_truth_scores_bin_" + bin_size_str + ".txt";
    
    std::ofstream of(std::filesystem::path(path_string_out) / filtered_filename);
    std::ofstream txt1(std::filesystem::path(path_string_out) / all_scores_filename);
    std::ofstream txt2(std::filesystem::path(path_string_out) / not_gt_scores_filename);
    
    std::string buffer;
    for (size_t i = 0; i < experimental_spectra_.size(); ++i)
    {
        readEntryIntoBufferMgf(f, buffer);

        if(experimental_spectra_[i].getMatch().is_initialized_m)
        {
            txt1 << experimental_spectra_[i].getMatch().tanimoto_m << "\n";
            if (experimental_spectra_[i].is_gt) txt2 << experimental_spectra_[i].getMatch().tanimoto_m << "\n";
        }
        if(experimental_spectra_[i].getMatch().is_initialized_m && experimental_spectra_[i].getMatch().tanimoto_m > cutoff)
        {
            of << buffer;
        }
        else total_filtered_ += 1;
    }

    f.close();
    of.close();
    txt1.close();
    txt2.close();
}