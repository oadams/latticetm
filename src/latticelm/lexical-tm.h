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

  LexicalTM(
      SymbolSet<std::string> f_vocab,
      SymbolSet<std::string> e_vocab,
      float alpha, float gamma,
      const std::unordered_set<std::string> & phonemes,
      std::string prior,
      float lambda,
      std::vector<float> starters) {
    f_vocab_size_ = f_vocab.size();
    e_vocab_size_ = e_vocab.size();
    phonemes_ = phonemes;
    f_vocab_ = f_vocab;
    e_vocab_ = e_vocab;
    log_alpha_ = LogWeight(-log(alpha));
    log_gamma_ = LogWeight(-log(gamma));
    gamma_ = gamma;

    // Zero the count std::vectors. Assign uniform log probabilities to the CPD
    for(int i=0; i < e_vocab_size_; i++) {
      std::vector<fst::LogWeight> cpd_accumulator_row;
      std::vector<fst::LogWeight> base_dist_row;
      std::vector<int> counts_row;
      for(int j=0; j < f_vocab_size_; j++) {
        cpd_accumulator_row.push_back(fst::LogWeight::Zero());
        base_dist_row.push_back(fst::LogWeight(-log(1.0/f_vocab_size_)));
        counts_row.push_back(0);
      }
      cpd_accumulator_.push_back(cpd_accumulator_row);
      base_dist_.push_back(base_dist_row);
      counts_.push_back(counts_row);
    }

    // Create the `empty' lexicon that allows for phonemes to pass through as-is.
    if(prior == "geom") {
      lexicon_ = CreateEmptyLexicon(phonemes_);
    } else if(prior == "pmp") {
      lexicon_ = CreateEmptyPMPLexicon(phonemes_, starters);
    } else if(prior == "poisson") {
      lexicon_ = CreateEmptyPoissonLexicon(phonemes_, lambda);
    }

  }

  void RemoveSample(const Alignment & align);
  void AddSample(const Alignment & align);
  Alignment CreateSample(const DataLattice & lat, LLStats & stats);
  void ResampleParameters();
  fst::VectorFst<fst::LogArc> CreateReducedTM(const DataLattice & lattice);
  fst::VectorFst<fst::LogArc> CreateReducedTM(const DataLattice & lattice, const std::vector<std::vector<fst::LogWeight>> & cpd);
  void FindBestPaths(const std::vector<DataLatticePtr> & lattices, std::string align_fn);
  void FindBestPaths(const std::vector<DataLatticePtr> & lattices, std::string align_fn, const std::vector<std::vector<fst::LogWeight>> cpd);
  void FindBestPaths(
      const std::vector<DataLatticePtr> & lattices,
      const std::string out_fn,
      // Will remove this dict, since LexicalTM has f_vocab_
      SymbolSet<std::string> & dict);
  void Normalize(int epochs);
  LogWeight DirichletProb(int e, int f);

  void WriteSortedCounts(std::string fn);

  // Test methods to be moved elsewhere later
  void TestLogWeightSampling();

  // Helpful methods
  void PrintParams(std::string path);
  void PrintAvgParams(std::string path);
  void PrintParams(std::vector<std::vector<fst::LogWeight>> cpd, std::string path);
  void PrintCounts();

  std::vector<std::vector<fst::LogWeight>> load_TM(const std::string filename);

  // Related to the phoneme-based extensions
  std::vector<std::string> GetPhonemes(const std::vector<DataLatticePtr> & lattices);
  VectorFst<LogArc> CreateEmptyLexicon(const std::unordered_set<std::string> & phonemes);
  VectorFst<LogArc> CreateEmptyPMPLexicon(
    const std::unordered_set<std::string> & phonemes,
    const std::vector<float> & starters);
  VectorFst<LogArc> CreateEmptyPoissonLexicon(
    const std::unordered_set<std::string> & phonemes, float lambda);
  VectorFst<LogArc> CreateTM(const DataLattice & lattice);
  void AddWord(VectorFst<LogArc> & lexicon, std::vector<WordId> phonemes, std::string phoneme_word);
  std::string PhonemeWord(std::vector<WordId> phonemes);
  Alignment PhonemeWordAlignment(const Alignment & ph_alignment);
  LogWeight DirichletProbNew(WordId e, WordId f);
  Alignment AssignUnks(const Alignment & unk_alignment, const Sentence & translation);
  void PrintAlign(const Alignment & align);
  void WriteSymbolSets();

protected:

  // Assuming our vocab fits in an int.
  int f_vocab_size_;
  int e_vocab_size_;
  SymbolSet<std::string> f_vocab_; // Used for both foreign words and phonemes.
  SymbolSet<std::string> e_vocab_;
  std::unordered_set<std::string> phonemes_;
  LogWeight log_alpha_; //Concentration parameter for the Dirichlet process.
  LogWeight log_gamma_; //Exponent used for the spelling model.
  float gamma_;

  // A grid that stores the sampling of the CPD at each iteration and gets
  // normalized after all the sampling is complete.
  std::vector<std::vector<fst::LogWeight>> cpd_accumulator_;
  // A uniform base distribution that the Dirichlet process will use.
  std::vector<std::vector<fst::LogWeight>> base_dist_;
  // The number of times we've seen a Foreign WordId align to an English WordId.
  std::vector<std::vector<int>> counts_;

  // Keys are pairs of foreign and English WordIds and values are couns of how
  // often the foreign word is aligned to the English word
  std::unordered_map<std::pair<WordId,WordId>, int> align_count_;
  // A map that stores the number of times a foreign word occurs.
  std::unordered_map<WordId, int> f_count_;
  // A map that stores the number of times an English word occurs.
  std::unordered_map<WordId, int> e_count_;

  VectorFst<LogArc> lexicon_;

  // The states associated with a given word in the lexicon.
  std::unordered_map<WordId, std::vector<VectorFst<LogArc>::StateId>> lexicon_states_;

};

}
