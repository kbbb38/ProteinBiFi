#include "CutoffAnalysis.h"
#include "ExperimentalSpectra.h"
#include "LibrarySpectra.h"

#include <iostream>
#include <fstream>
#include <bit>
#include <unordered_map>
#include <filesystem>

CutoffAnalysis::CutoffAnalysis(SpectrumBitSet& sbs, const std::string& path_string, const std::string& output_path, const float resolution)
{
  std::vector<ExperimentalSpectra>& exp = sbs.getExpSpec();
  std::vector<LibrarySpectra>& lib = sbs.getLibSpec();

  std::unordered_map<std::string, ExperimentalSpectra*> exp_m;
  std::unordered_map<std::string, LibrarySpectra*> lib_m;

  for (auto& e : exp) exp_m[e.getName()] = &e;
  for (auto& l : lib) lib_m[l.getPeptide()] = &l;

  std::ifstream f(path_string);
  
  // Create output filename based on resolution
  std::string scores_filename = "tanimoto_scores_res_" + std::to_string(resolution) + ".txt";
  std::ofstream of(std::filesystem::path(output_path) / scores_filename);
  
  std::string line;

  getline(f, line);
  std::string id;
  std::string sequence;

  int count = 0;
  while(getline(f, line, '\t'))
  {
    id = line;
    getline(f, line, '\n');
    sequence = line;

    if (exp_m.contains(id) && lib_m.contains(sequence))
    {
      float tmp_score = calculateTanimotoScore(exp_m[id]->getBitset(), exp_m[id]->getBitCount(), lib_m[sequence]->getBitset(), lib_m[sequence]->getBitCount());
      mean_score_ += tmp_score;
      ++count; 
      if (tmp_score < lowest_score_) lowest_score_ = tmp_score;
      of << tmp_score << "\n";
      exp_m[id]->is_gt = true;
    }
  }
  mean_score_ /= count;
  std::cout << "Lowest Tanimoto Score: " << lowest_score_ << std::endl;
  std::cout << "Mean Tanimoto Score: " << mean_score_ << std::endl;
  f.close();
  of.close();
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