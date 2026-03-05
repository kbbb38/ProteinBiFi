#pragma once

#include<cstdint>
#include<vector>
#include<string>

#include "Config.h"
#include "ExperimentalSpectra.h"
#include "LibrarySpectra.h"

class SpectrumBitSet 
{
    public:
        explicit SpectrumBitSet(AppConfig config);

        ~SpectrumBitSet() = default;

        SpectrumBitSet(const SpectrumBitSet&) = default;
        SpectrumBitSet& operator=(const SpectrumBitSet&) = default;

        SpectrumBitSet(SpectrumBitSet&& other) = default;
        SpectrumBitSet& operator=(SpectrumBitSet&& other) = default;

        void loadFile(const std::string& path_string);
        void matchSpectras();
        void writeOutput(const std::string& path_string) const;

        int loaded() { return total_loaded_; }
        
    private:
        std::vector<ExperimentalSpectra> experimental_spectra_;
        std::vector<LibrarySpectra> library_spectra_;
        int total_filtered_ = 0;
        int total_loaded_ = 0;
        bool is_dir_ = false;
        const AppConfig config_;

        void loadFromDirectory(const std::string& path_string);
        void loadSingleFile(const std::string& path_string);
        bool readEntryIntoBufferMgf(std::ifstream& f, std::string& buffer) const;
        bool readEntryIntoBufferMsp(std::ifstream& f, std::string& buffer) const;

        float calculateTanimotoScore(const std::vector<uint64_t>& e_spec, const uint64_t e_count, const std::vector<uint64_t>& l_spec, const uint64_t l_count) const;
};