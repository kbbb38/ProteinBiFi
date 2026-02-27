#pragma once

#include "SpectrumBitSet.h"
#include "config.h"

#include <string>
#include <vector>

class ExperimentalBitSets 
{
    public:
        explicit ExperimentalBitSets(const AppConfig& config);

        ~ExperimentalBitSets() = default;

        ExperimentalBitSets(const ExperimentalBitSets&) = default;
        ExperimentalBitSets& operator=(const ExperimentalBitSets&) = default;

        ExperimentalBitSets(ExperimentalBitSets&& other) = default;
        ExperimentalBitSets& operator=(ExperimentalBitSets&& other) = default;

        void loadExperimentalBitSets(const std::string& path_string);
        void filterExperimentalSpectra(const SpectrumBitSet& sbs);
        void saveFilteredExperimentalSpectras(const std::string& path_string);

    private:
        AppConfig config_;
        std::vector<std::vector<uint64_t>> experimental_bitsets_;
        std::vector<int> post_filter_indicies_;


        void loadFromDirectory(const std::string& path_string);
        void loadSingleFile(const std::string& path_string);
        bool readEntryIntoBuffer(std::ifstream& f, std::string& buffer) const;
        std::vector<float> readPeaksFromBuffer(const std::string& buffer) const;

        bool popCountBitsets(const std::vector<uint64_t>& exp, const std::vector<uint64_t>& lib, const double cutoff) const;

        void saveDirectory(const std::string& path_string);
        void saveSingleFile(const std::string& path_string)
};