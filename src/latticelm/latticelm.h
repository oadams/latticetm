#pragma once

#include <string>
#include <memory>
#include <unordered_map>
#include <latticelm/timer.h>
#include <latticelm/symbol-set.h>
#include <latticelm/sentence.h>
#include <latticelm/data-lattice.h>
#include <latticelm/lexical-tm.h>

namespace latticelm {


class LatticeLM {

public:
  LatticeLM()  { }

  int main(int argc, char** argv);

  template <class LM>
  void PerformTraining(const std::vector<DataLatticePtr> & lattices, LM & lm);
  void PerformTrainingLexTM(const std::vector<DataLatticePtr> & lattices, LexicalTM & lm, int train_len, int test_len);
  void Prototyping(const std::vector<DataLatticePtr> & lattices);
  
protected:

  std::string file_format_;
  std::string model_in_file_, model_out_file_;
  std::string model_type_;

  SymbolSet<std::string> cids_;
  SymbolSet<std::string> trans_ids_; // For the vocabulary of the translations.

  int epochs_, beam_;
  int char_n_, word_n_;
  float lattice_weight_;
  float alpha_; //A concentration parameter, if you want it.
  float gamma_; //A parameter for the prior spelling model (ie Geometric dist or poisson dist)
  std::string outfile_;

  Timer time_;

};

}
