#pragma once

#include<cstdint>
#include<vector>
#include<string>

class SpectrumBitSet 
{
    public:
        explicit SpectrumBitSet(const float resolution);

        ~SpectrumBitSet() = default;

        SpectrumBitSet(const SpectrumBitSet&) = default;
        SpectrumBitSet& operator=(const SpectrumBitSet&) = default;

        SpectrumBitSet(SpectrumBitSet&& other) = default;
        SpectrumBitSet& operator=(SpectrumBitSet&& other) = default;

        void loadBitSet(const std::string& path);

    private:
        const float resolution_;
        std::vector<uint64_t> bitset_;

        bool readEntryIntoBuffer(std::ifstream& f, std::string& buffer) const;
        std::vector<float> readPeaksFromBuffer(const std::string& buffer) const;
};