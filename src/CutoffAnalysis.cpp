#include "CutoffAnalysis.h"
#include "ExperimentalSpectra.h"
#include "LibrarySpectra.h"

#include <iostream>
#include <fstream>
#include <bit>

CutoffAnalysis::CutoffAnalysis(const SpectrumBitSet& sbs, const std::string& path_string)
{
  std::vector<ExperimentalSpectra> exp = sbs.getExpSpec();
  std::vector<LibrarySpectra> lib = sbs.getLibSpec();

  std::ifstream f(path_string);
  std::string line;

  getline(f, line);
  std::string id;
  std::string sequence;
  while(getline(f, line, '\t'))
  {
    id = line;
    getline(f, line, '\n');
    sequence = line;

    std::vector<uint64_t> bitset_e;
    std::vector<uint64_t> bitset_l;

    std::vector<uint64_t> tmp_bitset_l;
    int tmp_bitcount_l;
    for (auto l : lib)
    {
      if (l.getPeptide() == sequence)
      {
        for (auto e : exp)
        {
          if (e.getName() == id)
          {
            float tmp_score = calculateTanimotoScore(e.getBitset(), e.getBitCount(), tmp_bitset_l = l.getBitset(), l.getBitCount());
            if (tmp_score < lowest_score_) lowest_score_ = tmp_score;
            break;
          }
        }
        break;
      }
    }
  }
  std::cout << "Lowest Tanimoto Score: " << lowest_score_ << std::endl;
  f.close();
}

float CutoffAnalysis::calculateTanimotoScore(const std::vector<uint64_t>& e_spec, const uint64_t e_count, const std::vector<uint64_t>& l_spec, const uint64_t l_count) const
{
    uint64_t count_intersection = 0;
    const uint64_t* e_ptr = e_spec.data();
    const uint64_t* l_ptr = l_spec.data();
    const size_t size = e_spec.size();

    for(size_t i = 0; i < size; ++i)
    {
        count_intersection += std::popcount(e_ptr[i] & l_ptr[i]);
    }

    if (e_count == 0 && l_count == 0) return 0.0;
    uint64_t count_union = e_count + l_count - count_intersection;
    if (count_union == 0) return 0.0;

    return float(count_intersection) / float(count_union);
}