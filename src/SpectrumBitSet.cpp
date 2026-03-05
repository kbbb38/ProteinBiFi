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
    size_t total = experimental_spectra_.size();
    size_t current = 0;
    for (ExperimentalSpectra& e_spec : experimental_spectra_)
    {
        size_t index_count = 0;
        float highest_tanimoto = 0;
        int highest_tanimoto_id = 0;
        for (LibrarySpectra& l_spec : library_spectra_)
        {
            float tanimoto = calculateTanimotoScore(e_spec.getBitset(), e_spec.getBitCount(), l_spec.getBitset(), l_spec.getBitCount());
            if (tanimoto > highest_tanimoto)
            {
                highest_tanimoto = tanimoto;
                highest_tanimoto_id = index_count;
            }
            ++index_count;
        }
        library_spectra_[highest_tanimoto_id].setIfMatch();
        e_spec.setMatch(Match(highest_tanimoto, 0, library_spectra_[highest_tanimoto_id].getPeptide(), highest_tanimoto_id));
        ++current;
        int percent = (current * 100) / total;
        std::cout << "\r[" << std::string(percent / 2, '=') << std::string(50 - percent / 2, ' ') << "] " << percent << "% (" << current << "/" << total << ")" << std::flush;
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

    if (e_count == 0 && l_count == 0) return 1.0;
    uint64_t count_union = e_count + l_count - count_intersection;
    if (count_union == 0) return 0.0;

    return float(count_intersection) / float(count_union);
}

/* float SpectrumBitSet::calculateDotScore(const std::vector<uint64_t>& e_spec, const std::vector<uint64_t>& l_spec)
{

} */

void SpectrumBitSet::writeOutput(const std::string& path_string) const
{
    std::ofstream f(path_string);
    for (const ExperimentalSpectra& es : experimental_spectra_)
    {
        Match m = es.getMatch();
        f << es.getName() << "\t" << m.peptide_m << "\t" << m.tanimoto_m << "\t" << m.dot_product_m << "\n"; 
    }
    f.close();
}
