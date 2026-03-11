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

/* void SpectrumBitSet::matchSpectras()
{
    size_t total = experimental_spectra_.size();
    for (ExperimentalSpectra& e_spec : experimental_spectra_)
    {
        size_t index_count = 0;

        float highest_o = 0;
        int highest_id_o = 0;

        float highest_t = 0;
        int highest_id_t = 0;

        float highest_d = 0;
        int highest_id_d = 0;

        for (LibrarySpectra& l_spec : library_spectra_)
        {
            float score_t = calculateTanimotoScore(e_spec.getBitset(), e_spec.getBitCount(), l_spec.getBitset(), l_spec.getBitCount());
            float score_d = calculateDotProductScore(e_spec.getBitset(), e_spec.getIntensities(), l_spec.getBitset(), l_spec.getIntensities());
            float score_o = calculateOverlapCoefficient(e_spec.getBitset(), e_spec.getBitCount(), l_spec.getBitset(), l_spec.getBitCount());

            if (score_t > highest_t)
            {
                highest_t = score_t;
                highest_id_t = index_count;
            }

            if (score_o > highest_o)
            {
                highest_o = score_o;
                highest_id_o = index_count;
            }

            if (score_d > highest_d)
            {
                highest_d = score_d;
                highest_id_d = index_count;
            }

            ++index_count;
        }
        std::vector<std::string> matched_peptides = {library_spectra_[highest_id_t].getPeptide(), library_spectra_[highest_id_d].getPeptide(), library_spectra_[highest_id_o].getPeptide()};
        library_spectra_[highest_id_d].setIfMatch();
        e_spec.setMatch(Match(highest_t, highest_d, highest_o,
                                matched_peptides, 
                                highest_id_d));
    }
} */

void SpectrumBitSet::matchSpectras()
{
    size_t total = experimental_spectra_.size();
    
    for (ExperimentalSpectra& e_spec : experimental_spectra_)
    {
        // Temporary storage for all candidates for this experimental spectrum
        std::vector<Hit> all_candidates;
        all_candidates.reserve(library_spectra_.size());

        size_t index_count = 0;

        for (LibrarySpectra& l_spec : library_spectra_)
        {
            // Calculate all three scores for every pair
            float score_t = calculateTanimotoScore(e_spec.getBitset(), e_spec.getBitCount(), l_spec.getBitset(), l_spec.getBitCount());
            float score_d = calculateDotProductScore(e_spec.getBitset(), e_spec.getIntensities(), l_spec.getBitset(), l_spec.getIntensities());
            float score_o = calculateOverlapCoefficient(e_spec.getBitset(), e_spec.getBitCount(), l_spec.getBitset(), l_spec.getBitCount());

            Hit hit;
            hit.tanimoto_m = score_t;
            hit.dot_product_m = score_d;
            hit.overlap_coefficient_m = score_o;
            hit.peptide_m = l_spec.getPeptide();
            hit.library_id = index_count;

            all_candidates.push_back(hit);
            ++index_count;
        }

        // Prepare the Match object
        Match match_result;

        // Helper lambda to sort and trim to top 5
        auto processTop5 = [](std::vector<Hit>& candidates, auto compareFunc) {
            std::sort(candidates.begin(), candidates.end(), compareFunc);
            
            if (candidates.size() > 5) {
                candidates.resize(5);
            }
        };

        // 1. Process Tanimoto Top 5
        match_result.hits_tanimoto = all_candidates;
        processTop5(match_result.hits_tanimoto, [](const Hit& a, const Hit& b) {
            return a.tanimoto_m > b.tanimoto_m;
        });

        // 2. Process Overlap Top 5
        match_result.hits_overlap = all_candidates;
        processTop5(match_result.hits_overlap, [](const Hit& a, const Hit& b) {
            return a.overlap_coefficient_m > b.overlap_coefficient_m;
        });

        // 3. Process Dot Product Top 5
        match_result.dot_product = all_candidates;
        processTop5(match_result.dot_product, [](const Hit& a, const Hit& b) {
            return a.dot_product_m > b.dot_product_m;
        });

        // Maintain original side-effect logic: Mark the best Dot Product match as 'matched'
        // We assume the first element in the sorted dot_product vector is the best
        if (!match_result.dot_product.empty()) {
            size_t best_dot_id = match_result.dot_product[0].library_id;
            if (best_dot_id < library_spectra_.size()) {
                library_spectra_[best_dot_id].setIfMatch();
            }
        }

        // Set the match object to the experimental spectrum
        e_spec.setMatch(match_result);
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

/* void SpectrumBitSet::writeOutput(const std::string& path_string) const
{
    std::ofstream f(path_string);
    for (const ExperimentalSpectra& es : experimental_spectra_)
    {
        Match m = es.getMatch();
        f << es.getName() << "\t" << m.peptide_m[0] << "\t" << m.tanimoto_m << "\t" << m.peptide_m[1] << "\t" << m.overlap_coefficient_m << "\t" << m.peptide_m[2] << "\t" << m.dot_product_m << "\n"; 
    }
    f.close();
} */

void SpectrumBitSet::writeOutput(const std::string& path_string) const
{
    std::ofstream f(path_string);
    
    // Write header
    f << "Experimental_Spectrum\tMatch_Type\tRank\tPeptide\tTanimoto\tOverlap\tDot_Product\n";
    
    for (const ExperimentalSpectra& es : experimental_spectra_)
    {
        f << "START MATCHES" << "\n";

        Match m = es.getMatch();
        std::string spec_name = es.getName();
        
        // Output top 5 Tanimoto matches
        for (size_t i = 0; i < m.hits_tanimoto.size(); ++i)
        {
            const Hit& hit = m.hits_tanimoto[i];
            f << spec_name << "\t" 
              << "Tanimoto" << "\t" 
              << (i + 1) << "\t" 
              << hit.peptide_m << "\t" 
              << hit.tanimoto_m << "\t" 
              << hit.overlap_coefficient_m << "\t" 
              << hit.dot_product_m << "\n";
        }

        f << "\n";
        
        // Output top 5 Overlap matches
        for (size_t i = 0; i < m.hits_overlap.size(); ++i)
        {
            const Hit& hit = m.hits_overlap[i];
            f << spec_name << "\t" 
              << "Overlap" << "\t" 
              << (i + 1) << "\t" 
              << hit.peptide_m << "\t" 
              << hit.tanimoto_m << "\t" 
              << hit.overlap_coefficient_m << "\t" 
              << hit.dot_product_m << "\n";
        }
        
        f << "\n";

        // Output top 5 Dot Product matches
        for (size_t i = 0; i < m.dot_product.size(); ++i)
        {
            const Hit& hit = m.dot_product[i];
            f << spec_name << "\t" 
              << "Dot_Product" << "\t" 
              << (i + 1) << "\t" 
              << hit.peptide_m << "\t" 
              << hit.tanimoto_m << "\t" 
              << hit.overlap_coefficient_m << "\t" 
              << hit.dot_product_m << "\n";
        }
        f << "End MATCHES" << "\n";
    }
    
    f.close();
}