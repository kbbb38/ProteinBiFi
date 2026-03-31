#include "CutoffAnalysis.h"
#include "ExperimentalSpectra.h"
#include "LibrarySpectra.h"

#include <iostream>
#include <fstream>
#include <bit>
#include <unordered_map>

CutoffAnalysis::CutoffAnalysis(SpectrumBitSet& sbs, const std::string& ground_truth_path, const std::string& output_path)
{
  std::vector<ExperimentalSpectra>& exp = sbs.getExpSpec();
  std::vector<LibrarySpectra>& lib = sbs.getLibSpec();

  std::unordered_map<std::string, ExperimentalSpectra*> exp_m;
  std::unordered_map<std::string, LibrarySpectra*> lib_m;

  for (auto& e : exp) exp_m[e.getName()] = &e;
  for (auto& l : lib) lib_m[l.getPeptide()] = &l;

  std::ifstream ground_truth(ground_truth_path);
  std::ofstream ground_truth_output(output_path + "ground_truth_tanimotos.txt");
  std::string line;

  getline(ground_truth, line);
  std::string id;
  std::string sequence;

  int count = 0;
  while(getline(ground_truth, line, '\t'))
  {
    id = line;
    getline(ground_truth, line, '\n');
    sequence = line;

    if (exp_m.contains(id) && lib_m.contains(sequence))
    {
      float tmp_score = calculateTanimotoScore(exp_m[id]->getBitset(), exp_m[id]->getBitCount(), lib_m[sequence]->getBitset(), lib_m[sequence]->getBitCount());
      if (tmp_score < lowest_score_) lowest_score_ = tmp_score;
      ground_truth_output << tmp_score << "\n";

      mean_score_ += tmp_score;
      ++count; 
      
      exp_m[id]->setGroundTruth();
      lib_m[sequence]->setGroundTruth();
    }
  }
  mean_score_ /= count;

  std::cout << "Lowest Tanimoto Score: " << lowest_score_ << std::endl;
  std::cout << "Mean Tanimoto Score: " << mean_score_ << std::endl;

  ground_truth.close();
  ground_truth_output.close();
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