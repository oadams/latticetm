#include <latticelm/lexical-tm.h>
#include <latticelm/sampgen.h>
#include <fst/compose.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fst/shortest-path.h>
#include <latticelm/timer.h>

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
  cout << endl << "Avg. CPD parameters: " << endl;
  PrintParams(cpd_accumulator_);
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

  for(int f = 1; f < f_vocab_size_; f++) {
    reduced_tm.AddArc(only_state, LogArc(f, 0, cpd[0][f], only_state));
    for(int e = 1; e < f_vocab_size_; e++) {
      int times_in = in(e, translation);
      for(int i = 0; i < times_in; i++) {
        reduced_tm.AddArc(only_state, LogArc(f, e, cpd[e][f], only_state));
      }
    }
  }
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

/** Uses Dijkstra's shortest path algorithm to find the shortest path through
 * the lattice. This will be used for finding the best source sentence and
 * alignment in the composed lattice.

 *I'm implementing this because the FST ShortestPath implementation:
 *http://www.openfst.org/twiki/bin/view/FST/ShortestPathDoc indicates that the
 *`path' property must hold for the weights. But this path property does not
 *hold for log weights it would seem.*/
void LexicalTM::Dijkstra(const Fst<LogArc> & lattice, MutableFst<LogArc> * shortest_path) {
  VectorFst<LogArc>::StateId initial_state = lattice.Start();
  assert(initial_state == 0);
  //VectorFst<LogArc>::StateId final_state = lattice.NumStates()-1;

  vector<float> min_distance;
  min_distance.push_back(0.0);
  vector<int> prev_state;
  prev_state.push_back(-1);
  vector<pair<int,int>> prev_align;
  prev_align.push_back({-1,-1});
  set<pair<float,VectorFst<LogArc>::StateId>> active_vertices;
  active_vertices.insert( {0.0, initial_state} );

  while(!active_vertices.empty()) {
    int cur = active_vertices.begin()->second;
    active_vertices.erase(active_vertices.begin());
    fst::ArcIterator<Fst<LogArc>> arc_iter(lattice, cur);
    while(true) {
      if(arc_iter.Done()) break;
      const LogArc& arc = arc_iter.Value();
      //cout << arc.weight << " " << arc.ilabel << " " << arc.olabel << endl;
      // Expand min_distance if we need to.
      while(arc.nextstate+1 > min_distance.size()) {
        min_distance.push_back(std::numeric_limits<float>::max());
        prev_state.push_back(-1);
        pair<int,int> nullpair = {-1,-1};
        prev_align.push_back(nullpair);
      }
      if(fst::Times(min_distance[cur],arc.weight).Value() < min_distance[arc.nextstate]) {
        active_vertices.erase( { min_distance[arc.nextstate], arc.nextstate } );
        min_distance[arc.nextstate] = fst::Times(min_distance[cur], arc.weight).Value();
        prev_state[arc.nextstate] = cur;
        prev_align[arc.nextstate] = {arc.ilabel, arc.olabel};
        active_vertices.insert( { min_distance[arc.nextstate], arc.nextstate } );
      }
      arc_iter.Next();
    }
  }

  //cout << prev_state << endl;
  //cout << prev_align << endl;
  StringFromBacktrace(prev_state, prev_align);
  //cout << "Len of shortest path: " << min_distance[min_distance.size()-1] << endl;
  //cout << "---------------" << endl;
}

void LexicalTM::StringFromBacktrace(const vector<int> & prev_state, const vector<pair<int,int>> & prev_align) {
  int id = prev_state.size()-1;
  vector<string> foreign_source;
  while(true) {
    int wordid = prev_align[id].first;
    if(wordid == -1) break;
    foreign_source.push_back(f_vocab_.GetSym(wordid));
    id = prev_state[id];
  }
  for(int i = foreign_source.size()-1; i >= 0; i--){
      cout << foreign_source[i] << " ";
  }
  cout << endl;
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
      Dijkstra(vecfst, shortest_path);
      //shortest_path->Write("sample_" + to_string(i) + ".fst");
      i++;
  }

}

void LexicalTM::TestLogWeightSampling() {

  // Define some probabilities
  LogWeight eighty = LogWeight(0.22314);
  LogWeight twenty = LogWeight(1.60944);
  LogWeight sixty = LogWeight(0.51083);
  LogWeight ninetysix = LogWeight(0.04082);
  LogWeight five = LogWeight(-1.6094);

  // Create a trivial FST.
  VectorFst<LogArc> ifst;
  ifst.AddState();
  ifst.AddState();
  ifst.SetStart(0);
  ifst.SetFinal(1, LogArc::Weight::One());
  //ifst.AddArc(0, LogArc(1,1, eighty, 1));
  //ifst.AddArc(0, LogArc(2,2, twenty, 1));
  //ifst.AddArc(0, LogArc(1,1, fst::Times(eighty,eighty), 1));
  //ifst.AddArc(0, LogArc(2,2, fst::Times(twenty,twenty), 1));
  //ifst.AddArc(0, LogArc(1,1, sixty, 1));
  ifst.AddArc(0, LogArc(1,1, sixty, 1));
  ifst.AddArc(0, LogArc(2,2, twenty,1));
  ifst.AddArc(0, LogArc(2,2, twenty,1));

  cout << fst::Times(twenty,twenty) << endl;
  cout << fst::Plus(twenty,twenty) << endl;

  std::vector<int> count {0,0,0};
  for(int epoch=0; epoch < 10000; epoch++){
    VectorFst<LogArc> sample;
    SampGen(ifst, sample);
    Alignment align = FstToAlign(sample);
    count[align[0].first]++;
  }
  cout << count << endl;;

}
