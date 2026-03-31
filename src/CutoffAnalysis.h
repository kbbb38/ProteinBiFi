#include "SpectrumBitSet.h"

#include <vector>
#include <string>

class CutoffAnalysis
{
  private:
    float lowest_score_ = 1;
    float mean_score_ = 0;

    float calculateTanimotoScore(const std::vector<uint64_t>& e_spec, const uint64_t e_count, const std::vector<uint64_t>& l_spec, const uint64_t l_count) const;
    
  public:
    CutoffAnalysis(SpectrumBitSet& sbs, const std::string& path_string, const std::string& output_path, const float resolution);
    void analyseWithGroundTruth();
    const float getLowest() {return lowest_score_; }
    const float getMean() { return mean_score_; }
};