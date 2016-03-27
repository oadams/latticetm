#include <latticelm/lexical-tm.h>
#include <latticelm/sampgen.h>
#include <fst/compose.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fst/shortest-path.h>
#include <latticelm/timer.h>
#include <latticelm/data-lattice.h>
#include <algorithm>
#include <boost/bind.hpp>
#include <boost/algorithm/string.hpp>

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

void LexicalTM::PrintParams(string path) {
  PrintParams(cpd_accumulator_, path);
}

void LexicalTM::PrintParams(vector<vector<fst::LogWeight>> cpd, string path) {
  ofstream tm_file;
  tm_file.open(path);

  // %TODO I want to sort the English words by frequency, then show the top foreign translations of each.
  for(int i = 0; i < e_vocab_size_; i++) {
    tm_file << i << ", " << e_vocab_.GetSym(i) << endl;
    vector<pair<int, float>> dist(f_vocab_size_);
    for (int j = 0; j < f_vocab_size_; j++) {
      dist[j] = {j, cpd[i][j].Value()};
    }
    //vector<fst::LogWeight> dist = vector<fst::LogWeight>(cpd[i]);
    std::sort(dist.begin(), dist.end(), boost::bind(&std::pair<int, float>::second, _1) < boost::bind(&std::pair<int, float>::second, _2));
    //for(int j = 0; j < (10 < f_vocab_size_) ? 10 : f_vocab_size_; j++) {
    for(int j = 0; j < 7; j++) {
      tm_file << "\t" << dist[j].first << ", " << f_vocab_.GetSym(dist[j].first) << ": " << dist[j].second << endl;
    }
  }

  tm_file.close();

  /*
  tm_file << std::fixed << std::setw( 1 ) << std::setprecision( 3 );
  //tm_file << endl << "CPD parameters: " << endl;
  tm_file << "\t";
  for(int j = 0; j < f_vocab_size_; j++) {
    tm_file << f_vocab_.GetSym(j) << "\t";
  }
  tm_file << endl;
  for(int i = 0; i < e_vocab_size_; i++) {
    tm_file << e_vocab_.GetSym(i) << "\t";
    for(int j = 0; j < f_vocab_size_; j++) {
      tm_file << exp(-1*cpd[i][j].Value()) << "\t";
    }
    tm_file << endl;
  }
  tm_file << endl;
  */
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

LogWeight LexicalTM::DirichletProb(int e, int f) {
    // %TODO test the correctness of this function

    // Get the total counts of e.
    int e_total = 0;
    for(int f_prime = 0; f_prime < f_vocab_size_; f_prime++) {
      e_total += counts_[e][f_prime];
    }
    LogWeight numerator = fst::Plus(fst::Times(log_alpha_,base_dist_[e][f]), LogWeight(-log(counts_[e][f])));
    LogWeight denominator = fst::Plus(log_alpha_,LogWeight(-log(e_total)));
    return fst::Divide(numerator,denominator);
}

VectorFst<LogArc> LexicalTM::CreateReducedTM(const DataLattice & lattice) {
  VectorFst<LogArc> reduced_tm;
  VectorFst<LogArc>::StateId only_state = reduced_tm.AddState();
  reduced_tm.SetStart(only_state);
  reduced_tm.SetFinal(only_state, LogArc::Weight::One());

  Sentence translation = lattice.GetTranslation();

  std::set<WordId> translation_set(translation.begin(), translation.end());

  for(int e : translation_set) {
    LogWeight total = LogWeight::Zero();
    for(int f : lattice.GetFWordIds()) {
      total = fst::Plus(total, DirichletProb(e,f));
    }
    for(int f : lattice.GetFWordIds()) {
      if(total == LogWeight::Zero()) {
        reduced_tm.AddArc(only_state, LogArc(f, e, LogWeight::Zero(), only_state));
      } else {
        reduced_tm.AddArc(only_state, LogArc(f, e, fst::Divide(DirichletProb(e,f), total), only_state));
      }
    }
  }
  //ArcSortFst<LogArc, ILabelCompare<LogArc>>(reduced_tm, ILabelCompare<LogArc>());
  ArcSort(&reduced_tm, ILabelCompare<LogArc>());
  return reduced_tm;
}

/** Create a TM based on the parameters that is constrained by the lattice's translation **/
VectorFst<LogArc> LexicalTM::CreateReducedTM(const DataLattice & lattice, const vector<vector<fst::LogWeight>> & cpd) {
  VectorFst<LogArc> reduced_tm;
  VectorFst<LogArc>::StateId only_state = reduced_tm.AddState();
  reduced_tm.SetStart(only_state);
  reduced_tm.SetFinal(only_state, LogArc::Weight::One());

  Sentence translation = lattice.GetTranslation();

  std::set<WordId> translation_set(translation.begin(), translation.end());

  for(int e : translation_set) {
    LogWeight total = LogWeight::Zero();
    for(int f : lattice.GetFWordIds()) {
      total = fst::Plus(total, cpd[e][f]);
    }
    for(int f : lattice.GetFWordIds()) {
      if(total == LogWeight::Zero()) {
        reduced_tm.AddArc(only_state, LogArc(f, e, LogWeight::Zero(), only_state));
      } else {
        reduced_tm.AddArc(only_state, LogArc(f, e, fst::Divide(cpd[e][f], total), only_state));
      }
    }
  }
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

  // Perform reduction on TM to make it conform to the lattice.translation_
  VectorFst<LogArc> reduced_tm = CreateReducedTM(lattice);
  //reduced_tm.Write("reduced_tm.fst");

  //lattice.GetFst().Write("lattice.fst");

  // Compose the lattice with the reduced tm.
  ComposeFst<LogArc> composed_fst(lattice.GetFst(), reduced_tm);
  //VectorFst<LogArc> vecfst(composed_fst);
  //vecfst.Write("composed.fst");

  // Sample from the composed Fst.
  VectorFst<LogArc> sample_fst;
  /*stats.lik_ +=*/ SampGen(composed_fst, sample_fst);
  //sample_fst.Write("sample.fst");

  Alignment align = FstToAlign(sample_fst);

  return align;

}

void LexicalTM::ResampleParameters() {
  // Specify hyperparameters of the Dirichlet Process.
  // We assume a uniform distribution, base_dist_, which has been initialized to uniform.
  for(int i = 0; i < e_vocab_size_; i++) {
    double row_total = 0;
    for(int j = 0; j < f_vocab_size_; j++) {
      row_total += counts_[i][j];
    }
    for(int j = 0; j < f_vocab_size_; j++) {
      LogWeight numerator = fst::Plus(fst::Times(log_alpha_,base_dist_[i][j]), LogWeight(-log(counts_[i][j])));
      LogWeight denominator = fst::Plus(log_alpha_,LogWeight(-log(row_total)));
      cpd_accumulator_[i][j] = fst::Plus(cpd_accumulator_[i][j], fst::Divide(numerator,denominator));
    }
  }
}

void LexicalTM::FindBestPaths(const vector<DataLatticePtr> & lattices, string align_fn) {
  FindBestPaths(lattices, align_fn, cpd_accumulator_);
}

/** Samples the best path through the lattice using the translations and
* average translation model parameters to inform the sample. **/
void LexicalTM::FindBestPaths(
    const vector<DataLatticePtr> & lattices,
    string align_fn,
    vector<vector<fst::LogWeight>> cpd) {

  ofstream && align_file = ofstream();
  align_file.open(align_fn);

  long i = 0;
  for (auto latticep : lattices) {
      DataLattice lattice = *latticep;

      // Perform reduction on TM to make it conform to the lattice.translation_
      VectorFst<LogArc> reduced_tm = CreateReducedTM(lattice, cpd);
      reduced_tm.Write("data/out/lattices/reduced_tm" + to_string(i) + ".fst");
      // Compose the lattice with the reduced tm.
      ComposeFst<LogArc> composed_fst(lattice.GetFst(), reduced_tm);
      VectorFst<LogArc> vecfst(composed_fst);
      lattice.GetFst().Write("data/out/lattices/plain" + to_string(i) + ".fst");
      vecfst.Write("data/out/lattices/composed" + to_string(i++) + ".fst");
      // Find the shortest path.
      VectorFst<LogArc> * shortest_path = new VectorFst<LogArc>;
      vector<int> && prev_state = vector<int>();
      vector<pair<int,int>> && prev_align = vector<pair<int,int>>();
      DataLattice::Dijkstra(vecfst, prev_state, prev_align, f_vocab_, e_vocab_, false);

      int final_state = DataLattice::GetFinal(composed_fst);
      DataLattice::StringFromBacktrace(final_state, prev_state, prev_align, f_vocab_, cout);
      DataLattice::AlignmentFromBacktrace(final_state, prev_state, prev_align, f_vocab_, e_vocab_, align_file);
  }
  align_file.close();
}

void LexicalTM::FindBestPlainLatticePaths(const vector<DataLatticePtr> & lattices, string out_fn) {
    ofstream && out_file = ofstream();
    out_file.open(out_fn);

    long i = 0;
    for (auto latticep : lattices) {
      DataLattice lattice = *latticep;
      vector<int> && prev_state = vector<int>();
      vector<pair<int,int>> && prev_align = vector<pair<int,int>>();
      if(i == -8) {
        DataLattice::Dijkstra(lattice.GetFst(), prev_state, prev_align, f_vocab_, e_vocab_, true);
      } else {
        DataLattice::Dijkstra(lattice.GetFst(), prev_state, prev_align, f_vocab_, e_vocab_);
      }
      DataLattice::StringFromBacktrace(lattice.GetFinal(), prev_state, prev_align, f_vocab_, out_file);
      i++;
    }
    out_file.close();
}

vector<vector<fst::LogWeight>> LexicalTM::load_TM(const string filename) {
  ifstream in(filename);
  if(!in) THROW_ERROR("Could not open " << filename);

  // Initialize translation model
  vector<vector<fst::LogWeight>> tm;
  for(int e = 0; e < e_vocab_size_; e++) {
    vector<fst::LogWeight> row;
    for(int f = 0; f < f_vocab_size_; f++) {
      row.push_back(LogWeight::Zero());
    }
    tm.push_back(row);
  }

  string line;
  while(getline(in, line)) {
    vector<string> line_tokens;
    boost::split(line_tokens, line, boost::is_any_of(" "), boost::token_compress_on);
    if(line_tokens.size() != 3) {
      cerr << "tm line: " << line_tokens << endl;
      THROW_ERROR("Translation model line must consist of 3 tokens space delimited.");
    }
    string f = line_tokens[0];
    string e = line_tokens[1];
    float real_prob = stof(line_tokens[2]);
    LogWeight prob = LogWeight(-1*log(real_prob));

    if(real_prob <= 0.f) { 
      THROW_ERROR("TM probs loaded must be greater than 0.0");
    }

//    cout << "before suspect." << endl;
 //   cout << "e" << endl;
    if(e_vocab_.KeyInMap(e) && f_vocab_.KeyInMap(f)) {
      tm[(e != "NULL") ? e_vocab_.GetId(e) : 0][(f != "NULL") ? f_vocab_.GetId(f) : 0] = prob;
    }
  //  cout << "after suspect." << endl;
  }

  //cout << "P(super|student): " << tm[e_vocab_.GetId("student")][f_vocab_.GetId("super")] << endl;
  //cout << fst::Plus(tm[e_vocab_.GetId("student")][f_vocab_.GetId("super")],fst::LogWeight(1.3404)) << endl;

  // Check to make sure each row is a valid distribution
  /*
  for(vector<fst::LogWeight> row : tm) {
    fst::LogWeight total = LogWeight::Zero();
    for(fst::LogWeight prob : row) {
      total = fst::Plus(total, prob);
    }
    cout << total << endl;
  }
  */

  return tm;
}
