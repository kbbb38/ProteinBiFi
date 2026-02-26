#pragma once

#include "SpectrumBitSet.h"

#include <string>
#include <vector>

class ExperimentalBitSets 
{
    public:
        explicit ExperimentalBitSets(const float resolution);

        ~ExperimentalBitSets() = default;

        ExperimentalBitSets(const ExperimentalBitSets&) = default;
        ExperimentalBitSets& operator=(const ExperimentalBitSets&) = default;

        ExperimentalBitSets(ExperimentalBitSets&& other) = default;
        ExperimentalBitSets& operator=(ExperimentalBitSets&& other) = default;

        void loadExperimentalBitSets(const std::string& path);
        void filterExperimentalSpectra(const SpectrumBitSet& sbs);

    private:
        float resolution_;
        std::vector<std::vector<uint64_t>> experimental_bitsets_;
        std::vector<int> post_filter_indicies_;

        void loadFromDirectory(const std::string& path);
        void loadSingleFile(const std::string& path);
        bool readEntryIntoBuffer(std::ifstream& f, std::string& buffer) const;
        std::vector<float> readPeaksFromBuffer(const std::string& buffer) const;

        bool popCountBitsets(const std::vector<uint64_t>& exp, const std::vector<uint64_t>& lib) const;
};