#pragma once

#include<cstdint>
#include<vector>
#include<string>

#include "Config.h"

class SpectrumBitSet 
{
    public:
        explicit SpectrumBitSet(const float resolution);

        ~SpectrumBitSet() = default;

        SpectrumBitSet(const SpectrumBitSet&) = default;
        SpectrumBitSet& operator=(const SpectrumBitSet&) = default;

        SpectrumBitSet(SpectrumBitSet&& other) = default;
        SpectrumBitSet& operator=(SpectrumBitSet&& other) = default;

        void loadFile(const std::string& path);
        void filterSpectra();
        int filtered() { return total_filtered_; };

        const std::vector<uint64_t>& bitset() const { return bitset_; }
        
    private:
        std::vector<uint64_t> bitset_;
        int total_filtered_ = 0;
        AppConfig config_;

        void loadFromDirectory(const std::string& path_string);
        void loadSingleFile(const std::string& path_string);

        bool readEntryIntoBufferExp(std::ifstream& f, std::string& buffer) const;
        std::vector<float> readPeaksFromBufferExp(const std::string& buffer) const;

        bool readEntryIntoBufferLib(std::ifstream& f, std::string& buffer) const;
        std::vector<float> readPeaksFromBufferLib(const std::string& buffer) const;

        void addPeaksToBitset(const std::vector<float>& tmp_peaks);
        std::vector<uint64_t> readPeaksIntoBitset(const std::vector<float>& tmp_peak) const;

        bool popCountBitsets(const std::vector<uint64_t>& tmp_bitset) const;
};