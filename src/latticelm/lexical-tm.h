#pragma once

#include <latticelm/ll-stats.h>
#include <latticelm/sentence.h>
#include <latticelm/data-lattice.h>
#include <fst/vector-fst.h>
#include <fst/float-weight.h>
#include <cmath>
#include <iostream>
#include <unordered_map>

namespace latticelm {

class LexicalTM {

public:

  LexicalTM( SymbolSet<std::string> f_vocab, SymbolSet<std::string> e_vocab,
      float alpha, float gamma, const vector<string> & phonemes) {
    f_vocab_size_ = f_vocab.size();
    e_vocab_size_ = e_vocab.size();
    f_vocab_ = f_vocab;
    e_vocab_ = e_vocab;
    log_alpha_ = LogWeight(-log(alpha));
    log_gamma_ = LogWeight(-log(gamma));

    // Zero the count vectors. Assign uniform log probabilities to the CPD
    for(int i=0; i < e_vocab_size_; i++) {
      vector<fst::LogWeight> cpd_accumulator_row;
      vector<fst::LogWeight> base_dist_row;
      vector<int> counts_row;
      for(int j=0; j < f_vocab_size_; j++) {
        cpd_accumulator_row.push_back(fst::LogWeight::Zero());
        base_dist_row.push_back(fst::LogWeight(-log(1.0/f_vocab_size_)));
        counts_row.push_back(0);
      }
      cpd_accumulator_.push_back(cpd_accumulator_row);
      base_dist_.push_back(base_dist_row);
      counts_.push_back(counts_row);
    }

    phonemes_ = phonemes;
    // Create the `empty' lexicon that allows for phonemes to pass through as-is.
    lexicon_ = CreateEmptyLexicon(phonemes);

  }

  void RemoveSample(const Alignment & align);
  void AddSample(const Alignment & align);
  Alignment CreateSample(const DataLattice & lat, LLStats & stats);
  void ResampleParameters();
  fst::VectorFst<fst::LogArc> CreateReducedTM(const DataLattice & lattice);
  fst::VectorFst<fst::LogArc> CreateReducedTM(const DataLattice & lattice, const vector<vector<fst::LogWeight>> & cpd);
  void FindBestPaths(const vector<DataLatticePtr> & lattices, string align_fn);
  void FindBestPaths(const vector<DataLatticePtr> & lattices, string align_fn, const vector<vector<fst::LogWeight>> cpd);
  void FindBestPlainLatticePaths(const vector<DataLatticePtr> & lattices, string out_fn);
  void Normalize(int epochs);
  LogWeight DirichletProb(int e, int f);

  // Test methods to be moved elsewhere later
  void TestLogWeightSampling();

  // Helpful methods
  void PrintParams(string path);
  void PrintAvgParams(string path);
  void PrintParams(vector<vector<fst::LogWeight>> cpd, string path);
  void PrintCounts();

  vector<vector<fst::LogWeight>> load_TM(const string filename);

  // Related to the phoneme-based extensions
  vector<string> GetPhonemes(const vector<DataLatticePtr> & lattices);
  VectorFst<LogArc> CreateEmptyLexicon(const vector<string> & phonemes);
  VectorFst<LogArc> CreateTM(const DataLattice & lattice);
  void AddWord(VectorFst<LogArc> & lexicon, vector<WordId> phonemes);
  std::string PhonemeWord(vector<WordId> phonemes);
  Alignment PhonemeWordAlignment(const Alignment & ph_alignment);
  LogWeight DirichletProbNew(WordId e, WordId f);

protected:

  // Assuming our vocab fits in an int.
  int f_vocab_size_;
  int e_vocab_size_;
  SymbolSet<std::string> f_vocab_;
  SymbolSet<std::string> e_vocab_;
  vector<std::string> phonemes_;
  LogWeight log_alpha_; //Concentration parameter for the Dirichlet process.
  LogWeight log_gamma_; //Exponent used for the spelling model.

  // A grid that stores the sampling of the CPD at each iteration and gets
  // normalized after all the sampling is complete.
  vector<vector<fst::LogWeight>> cpd_accumulator_;
  // A uniform base distribution that the Dirichlet process will use.
  vector<vector<fst::LogWeight>> base_dist_;
  // The number of times we've seen a Foreign WordId align to an English WordId.
  vector<vector<int>> counts_;

  // Keys are pairs of foreign and English WordIds and values are couns of how
  // often the foreign word is aligned to the English word
  std::unordered_map<std::pair<WordId,WordId>, int> count_map_;
  // A map that stores the number of times a foreign word occurs.
  std::unordered_map<WordId, int> foreign_count_map_;

  VectorFst<LogArc> lexicon_;

};

}
