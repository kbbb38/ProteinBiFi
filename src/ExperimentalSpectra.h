#pragma once

#include "Config.h"

#include <string>
#include <vector>

struct Hit
{
    float tanimoto_m;
    float dot_product_m;
    float overlap_coefficient_m;
    std::string peptide_m;
    size_t library_id;
};


struct Match
{
    std::vector<Hit> hits_tanimoto;
    std::vector<Hit> hits_overlap;
    std::vector<Hit> dot_product;
};

class ExperimentalSpectra 
{
private:
    AppConfig config_;

    std::vector<float> peak_positions_;
    std::vector<float> intensities_;

    std::string name_;
    std::vector<uint64_t> bitset_;
    std::vector<float> binned_intensities_; 
    u_int64_t bit_count_;
    int charge_;
    float pepmass_;

    Match match_;

    void createBitSet();
    void binIntensities();
    void normalizeAndScaleIntensities();
    
public:
    ExperimentalSpectra() = default;
    ExperimentalSpectra(const std::string& b, const AppConfig& config);

    const std::string& getName() const { return name_; }
    const std::vector<uint64_t>& getBitset() const { return bitset_; }
    const uint64_t getBitCount() const { return bit_count_; }
    const Match getMatch() const { return match_; }
    const std::vector<float>& getIntensities() const { return binned_intensities_; }
    const std::vector<float>& getPeakPositions() const { return peak_positions_; }
    const int getCharge() const { return charge_; };
    const float getMass() const { return pepmass_; };

    void setName(const std::string& n) { name_ = n; }
    void setBitset(const std::vector<uint64_t>& bs) { bitset_ = bs; }
    void setMatch(const Match& m) {match_ = m; }
};