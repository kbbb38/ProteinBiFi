#pragma once

#include "Config.h"

#include <string>
#include <vector>

struct Match
{
    float tanimoto_m;
    float dot_product_m;
    float overlap_coefficient_m;
    std::string peptide_m;
    int id_m = 0;
};

class ExperimentalSpectra {
private:
    AppConfig config_;

    std::vector<float> peak_positions_;
    std::vector<float> intensities_;

    std::string name_;
    std::vector<uint64_t> bitset_;
    u_int64_t bit_count_;

    Match match_;

    void createBitSet();

public:
    ExperimentalSpectra() = default;
    ExperimentalSpectra(const std::string& b, const AppConfig& config);

    const std::string& getName() const { return name_; }
    const std::vector<uint64_t>& getBitset() const { return bitset_; }
    const uint64_t getBitCount() const { return bit_count_; }
    const Match getMatch() const { return match_; }

    void setName(const std::string& n) { name_ = n; }
    void setBitset(const std::vector<uint64_t>& bs) { bitset_ = bs; }
    void setMatch(const Match& m) {match_ = m; }
};