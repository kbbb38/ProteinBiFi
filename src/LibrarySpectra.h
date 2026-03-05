#pragma once

#include "Config.h"

#include <string>
#include <vector>

class LibrarySpectra {
    private:
        AppConfig config_;

        std::string peptide_;
        std::vector<uint64_t> bitset_;
        u_int64_t bit_count_;
        bool is_a_match_ = false;

        void createBitSet(const std::vector<float>& tmp_peaks);

    public:
        LibrarySpectra() = default;
        LibrarySpectra(const std::string& b, const AppConfig& config);

        const std::string& getPeptide() const { return peptide_; }
        const std::vector<uint64_t>& getBitset() const { return bitset_; }
        const uint64_t getBitCount() const { return bit_count_; }
        const bool getIfMatch() const { return is_a_match_; }

        void setPeptide(const std::string& n) { peptide_ = n; }
        void setBitset(const std::vector<uint64_t>& bs) { bitset_ = bs; }
        void setIfMatch() { is_a_match_ = true; }
};