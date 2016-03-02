#include <latticelm/lexical-tm.h>
#include <latticelm/sampgen.h>
#include <fst/compose.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fst/shortest-path.h>
#include <latticelm/timer.h>
#include <latticelm/data-lattice.h>

using namespace latticelm;
using namespace fst;

void LexicalTM::RemoveSample(const Alignment & align) {
  //Reduce the counts for the alignments.
  for(int i = 0; i < align.size(); i++) {
    counts_[align[i].second][align[i].first]--;
    assert(counts_[align[i].second][align[i].first] >= 0);
  }
}

void LexicalTM::AddSample(const Alignment & align) {
  //Reduce the counts for the alignments.
  for(int i = 0; i < align.size(); i++) {
    counts_[align[i].second][align[i].first]++;
    assert(counts_[align[i].second][align[i].first] > 0);
  }
}

void LexicalTM::PrintCounts() {
  cout << endl << "Alignment counts: " << endl;
  cout << "\t";
  for(int j = 0; j < f_vocab_size_; j++) {
    cout << f_vocab_.GetSym(j) << "\t";
  }
  cout << endl;
  for(int i = 0; i < e_vocab_size_; i++) {
    cout << e_vocab_.GetSym(i) << "\t";
    for(int j = 0; j < f_vocab_size_; j++) {
      cout << counts_[i][j] << "\t";
    }
    cout << endl;
  }
  cout << endl;
}

void LexicalTM::PrintParams() {
  PrintParams(cpd_);
}

void LexicalTM::PrintParams(vector<vector<fst::LogWeight>> cpd) {
  cout << std::fixed << std::setw( 1 ) << std::setprecision( 3 );
  //cout << endl << "CPD parameters: " << endl;
  cout << "\t";
  for(int j = 0; j < f_vocab_size_; j++) {
    cout << f_vocab_.GetSym(j) << "\t";
  }
  cout << endl;
  for(int i = 0; i < e_vocab_size_; i++) {
    cout << e_vocab_.GetSym(i) << "\t";
    for(int j = 0; j < f_vocab_size_; j++) {
      cout << exp(-1*cpd[i][j].Value()) << "\t";
    }
    cout << endl;
  }
  cout << endl;
}

void LexicalTM::Normalize(int epochs) {
  for(int i = 0; i < e_vocab_size_; i++) {
    for(int j = 0; j < f_vocab_size_; j++) {
      cpd_accumulator_[i][j] = fst::Divide(cpd_accumulator_[i][j],LogWeight(-log(epochs)));
    }
  }
  //cout << endl << "Avg. CPD parameters: " << endl;
  //PrintParams(cpd_accumulator_);
  /*
  cout << std::fixed << std::setw( 1 ) << std::setprecision( 3 );
  cout << endl << "Average CPD parameters: " << endl;
  cout << "\t";
  for(int j = 0; j < f_vocab_size_; j++) {
    cout << f_vocab_.GetSym(j) << "\t";
  }
  cout << endl;
  for(int i = 0; i < e_vocab_size_; i++) {
    cout << e_vocab_.GetSym(i) << "\t";
    for(int j = 0; j < f_vocab_size_; j++) {
      cout << exp(-1*cpd_accumulator_[i][j].Value()) << "\t";
    }
    cout << endl;
  }
  cout << endl;
  */
}

/** Gives the number of times a word_id occurred in a sentence.*/
int in(WordId word_id, Sentence sentence) {
  int ret = 0;
  for(int i = 0; i < sentence.size(); i++) {
    if (word_id == sentence[i]) {
      ret++;
    }
  }
  return ret;
}

VectorFst<LogArc> LexicalTM::CreateReducedTM(const DataLattice & lattice) {
  return LexicalTM::CreateReducedTM(lattice, cpd_);
}

/** Create a TM based on the parameters that is constrained by the lattice's translation **/
VectorFst<LogArc> LexicalTM::CreateReducedTM(const DataLattice & lattice, const vector<vector<fst::LogWeight>> & cpd) {
  VectorFst<LogArc> reduced_tm;
  VectorFst<LogArc>::StateId only_state = reduced_tm.AddState();
  reduced_tm.SetStart(only_state);
  reduced_tm.SetFinal(only_state, LogArc::Weight::One());

  Sentence translation = lattice.GetTranslation();

  LogWeight total = LogWeight::Zero();
  for(int f : lattice.GetFWordIds()) {
    total = fst::Plus(total, cpd[0][f]);
  }
  for(int f : lattice.GetFWordIds()) {
    reduced_tm.AddArc(only_state, LogArc(f, 0, fst::Divide(cpd[0][f], total), only_state));
  }
  for(int e = 1; e < e_vocab_size_; e++) {
    LogWeight total = LogWeight::Zero();
    int times_in = in(e, translation);
    for(int i = 0; i < times_in; i++) {
      for(int f : lattice.GetFWordIds()) {
        total = fst::Plus(total, cpd[e][f]);
      }
    }
    if(times_in > 0) {
      for(int f : lattice.GetFWordIds()) {
        LogWeight dupCoef = fst::LogWeight(-1*log(times_in)); // So that we can multiply the weight of the arc by the number of times we see the English word.
        reduced_tm.AddArc(only_state, LogArc(f, e, fst::Divide(fst::Times(dupCoef, cpd[e][f]), total), only_state));
      }
    }
  }
  //ArcSortFst<LogArc, ILabelCompare<LogArc>>(reduced_tm, ILabelCompare<LogArc>());
  ArcSort(&reduced_tm, ILabelCompare<LogArc>());
  return reduced_tm;
}

/** Create a TM based on the parameters that is constrained by the lattice's translation **/
/*
VectorFst<LogArc> LexicalTM::CreateReducedTM(const DataLattice & lattice, const vector<vector<fst::LogWeight>> & cpd) {
  Timer time = Timer();
  VectorFst<LogArc> reduced_tm;
  VectorFst<LogArc>::StateId only_state = reduced_tm.AddState();
  reduced_tm.SetStart(only_state);
  reduced_tm.SetFinal(only_state, LogArc::Weight::One());

  Sentence translation = lattice.GetTranslation();

  cerr << "preloop: " << time.Elapsed() << endl;

  // %TODO: Perhaps this should be optimized at some point.

  // Starting at 1 because 0 represents an epsilon transition and we don't
  // accept epsilon transitions on the foreign side in the translation model.
  // That would result in loops in the composition.
  for(int f_word_id = 1; f_word_id < f_vocab_size_; f_word_id++) {
    // Normalizing the probabilities. Two steps:
    // 1. Find the total probability mass of the p(f|e) for each of the English words that occur in
    //    the translation given the foreign word.
    LogWeight total = LogWeight::Zero();
    // First add the probability of an epsilon (ie. null token) on the English side.
    total = fst::Plus(total, cpd[0][f_word_id]);
    // Then check each of the English words to see if they are in the
    // translation, and add probability mass if they are
    for(int e_word_id = 1; e_word_id < e_vocab_size_; e_word_id++) {
      int times_in = in(e_word_id, translation);
      for(int i = 0; i < times_in; i++) {
        total = fst::Plus(total, cpd[e_word_id][f_word_id]);
      }
    }
    // 2. Divide the conditional probability of each of the English words by the
    //    aforementioned total when adding a corresponding arc to the WFST.
    reduced_tm.AddArc(only_state, LogArc(f_word_id, 0, fst::Divide(cpd[0][f_word_id], total), only_state));
    for(int e_word_id = 1; e_word_id < e_vocab_size_; e_word_id++) {
      int times_in = in(e_word_id, translation);
      for(int i = 0; i < times_in; i++) {
        reduced_tm.AddArc(only_state, LogArc(f_word_id, e_word_id, fst::Divide(cpd[e_word_id][f_word_id], total), only_state));
      }
    }
    cerr << "loop " << f_word_id << ": " << time.Elapsed() << endl;
  }

  return reduced_tm;
}
*/

Alignment LexicalTM::CreateSample(const DataLattice & lattice, LLStats & stats) {

  //TestLogWeightSampling();
  //exit(0);

  // Perform reduction on TM to make it conform to the lattice.translation_
  Timer time;
  time = Timer();
  VectorFst<LogArc> reduced_tm = CreateReducedTM(lattice);
  cerr << "Creating reduced tm took: " << time.Elapsed() << endl;
  reduced_tm.Write("reduced_tm.fst");

  //lattice.GetFst().Write("lattice.fst");

  // Compose the lattice with the reduced tm.
  ComposeFst<LogArc> composed_fst(lattice.GetFst(), reduced_tm);
  //VectorFst<LogArc> vecfst(composed_fst);
  //vecfst.Write("composed.fst");

  // Sample from the composed Fst.
  Timer time2;
  time2 = Timer();
  VectorFst<LogArc> sample_fst;
  /*stats.lik_ +=*/ SampGen(composed_fst, sample_fst);
  //sample_fst.Write("sample.fst");
  cerr << "Sampling took: " << time2.Elapsed() << endl;

  /*
  vector<int> counts = {0,0,0,0,0,0,0};
  for(int e = 0; e < 1000; e++) {
    VectorFst<LogArc> sample_fst;
    SampGen(composed_fst, sample_fst);
    Alignment align = FstToAlign(sample_fst);
    counts[align[0].first]++;
    counts[align[1].first]++;
  }
  cout << counts << endl;
  exit(0);
  */

  Timer time3 = Timer();
  Alignment align = FstToAlign(sample_fst);
  cerr << "FstToAlign took: " << time3.Elapsed() << endl;

  return align;

}

void LexicalTM::ResampleParameters() {
  // Specify hyperparameters of the Dirichlet Process.
  // We assume a uniform distribution, base_dist_, which has been initialized to uniform.
  for(int i = 0; i < e_vocab_size_; i++) {
    double row_total = 0;
    // TODO Perhaps superfluous? This row_total should just end up adding to 1.
    for(int j = 0; j < f_vocab_size_; j++) {
      row_total += counts_[i][j];
    }
    for(int j = 0; j < f_vocab_size_; j++) {
      LogWeight numerator = fst::Plus(fst::Times(log_alpha_,base_dist_[i][j]), LogWeight(-log(counts_[i][j])));
      LogWeight denominator = fst::Plus(log_alpha_,LogWeight(-log(row_total)));
      cpd_[i][j] = fst::Divide(numerator,denominator);
      cpd_accumulator_[i][j] = fst::Plus(cpd_accumulator_[i][j], cpd_[i][j]);
    }
  }
}


/** Samples the best path through the lattice using the translations and
* average translation model parameters to inform the sample. **/
void LexicalTM::FindBestPaths(const vector<DataLatticePtr> & lattices) {

  long i = 1;
  for (auto latticep : lattices) {
      DataLattice lattice = *latticep;
      // Perform reduction on TM to make it conform to the lattice.translation_
      // %TODO: Make this use the normalized TM params.
      VectorFst<LogArc> reduced_tm = CreateReducedTM(lattice, cpd_accumulator_);
      // Compose the lattice with the reduced tm.
      ComposeFst<LogArc> composed_fst(lattice.GetFst(), reduced_tm);
      VectorFst<LogArc> vecfst(composed_fst);
      vecfst.Write("composedforbestpath.fst");
      // Find the shortest path.
      VectorFst<LogArc> * shortest_path = new VectorFst<LogArc>;
      //VectorFst<StdArc> tropfst(vecfst);
      DataLattice::Dijkstra(vecfst, f_vocab_, e_vocab_);
      //shortest_path->Write("sample_" + to_string(i) + ".fst");
      i++;
  }

}
