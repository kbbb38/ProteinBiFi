#pragma once

#include "Config.h"

#include <string>
#include <vector>

class LibrarySpectra {
    private:
        AppConfig config_;

        std::vector<float> peak_positions_;
        std::vector<float> intensities_;

        std::vector<uint64_t> bitset_;
        std::vector<float> binned_intensities_; 

        std::string peptide_;
        u_int64_t bit_count_;
        bool is_a_match_ = false;
        int charge_;
        float pepmass_;

        void createBitSet();
        void binIntensities();
        void normalizeAndScaleIntensities();

    public:
        LibrarySpectra() = default;
        LibrarySpectra(const std::string& b, const AppConfig& config);

        const std::string& getPeptide() const { return peptide_; }
        const std::vector<uint64_t>& getBitset() const { return bitset_; }
        const uint64_t getBitCount() const { return bit_count_; }
        const bool getIfMatch() const { return is_a_match_; }
        const std::vector<float>& getIntensities() const { return binned_intensities_; }
        const std::vector<float>& getPeakPositions() const { return peak_positions_; }
        const int getCharge() const {return charge_; };
        const float getMass() const { return pepmass_; };

        void setPeptide(const std::string& n) { peptide_ = n; }
        void setBitset(const std::vector<uint64_t>& bs) { bitset_ = bs; }
        void setIfMatch() { is_a_match_ = true; }
};