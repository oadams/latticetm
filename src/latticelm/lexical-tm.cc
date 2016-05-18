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
#include <random>

using namespace latticelm;
using namespace fst;

/** Creates an empty lexicon with a Geometric prior**/
VectorFst<LogArc> LexicalTM::CreateEmptyLexicon(const unordered_set<string> & phonemes) {
  VectorFst<LogArc> lexicon;
  VectorFst<LogArc>::StateId home = lexicon.AddState();
  lexicon.SetStart(home);
  lexicon.SetFinal(home, LogArc::Weight::One());

  VectorFst<LogArc>::StateId phoneme_home = lexicon.AddState();

  for(auto it = phonemes.begin(); it != phonemes.end(); ++it ) {
    lexicon.AddArc(home, LogArc(f_vocab_.GetId(*it), f_vocab_.GetId("<eps>"), 
        Divide(LogWeight::One(), LogWeight(-log(phonemes_.size()))), phoneme_home));
    lexicon.AddArc(phoneme_home, LogArc(f_vocab_.GetId(*it), f_vocab_.GetId("<eps>"),
        Divide(log_gamma_, LogWeight(-log(phonemes_.size()))) , phoneme_home));
  }

  lexicon.AddArc(phoneme_home, LogArc(f_vocab_.GetId("<eps>"), f_vocab_.GetId("<unk>"),
      LogWeight(-log(1-gamma_)), home));
  lexicon.Write("data/phoneme-prototyping/lexicons/empty.fst");

  return lexicon;
}

/** Creates an empty lexicon with a poor man's Poisson prior. starters is a
 * list of k probabilities for P(1)..P(k). Everything thereafter is geometric
 * according to gamma. */
VectorFst<LogArc> LexicalTM::CreateEmptyPMPLexicon(
    const unordered_set<string> & phonemes,
    const vector<float> & starters) {

  VectorFst<LogArc> lexicon;
  VectorFst<LogArc>::StateId home = lexicon.AddState();
  lexicon.SetStart(home);
  lexicon.SetFinal(home, LogArc::Weight::One());

  // Add a state for each specified probability.
  VectorFst<LogArc>::StateId cur = lexicon.AddState();
  for(auto it = phonemes.begin(); it != phonemes.end(); ++it ) {
    lexicon.AddArc(home, LogArc(f_vocab_.GetId(*it), f_vocab_.GetId("<eps>"),
        Divide(LogWeight::One(), LogWeight(-log(phonemes_.size()))), cur));
  }

  LogWeight remainder = LogWeight::One();
  for(float prob : starters) {
    LogWeight end_prob = Divide(Divide(prob,remainder),
                                LogWeight(-log(phonemes_.size())));
    lexicon.AddArc(
        cur,
        LogArc(
            f_vocab_.GetId("<eps>"),
            f_vocab_.GetId("<unk>"),
            end_prob,
            home));
    remainder = LogWeight(-log(1-exp(-end_prob.Value())));
    VectorFst<LogArc>::StateId next = lexicon.AddState();
    for(auto it = phonemes.begin(); it != phonemes.end(); ++it) {
      lexicon.AddArc(
          cur,
          LogArc(
              f_vocab_.GetId(*it),
              f_vocab_.GetId("<eps>"),
              remainder,
              next));
    }
    cur = next;
  }

  // Add the Geometric component.
  VectorFst<LogArc>::StateId phoneme_home = lexicon.AddState();

  for(auto it = phonemes.begin(); it != phonemes.end(); ++it ) {
    lexicon.AddArc(cur, LogArc(f_vocab_.GetId(*it), f_vocab_.GetId("<eps>"), 
        Divide(log_gamma_, LogWeight(-log(phonemes_.size()))), phoneme_home));
    lexicon.AddArc(phoneme_home, LogArc(f_vocab_.GetId(*it), f_vocab_.GetId("<eps>"),
        Divide(log_gamma_, LogWeight(-log(phonemes_.size()))) , phoneme_home));
  }

  lexicon.AddArc(phoneme_home, LogArc(f_vocab_.GetId("<eps>"), f_vocab_.GetId("<unk>"),
      LogWeight(-log(1-gamma_)), home));
  lexicon.Write("data/phoneme-prototyping/lexicons/empty.fst");

  return lexicon;
}

VectorFst<LogArc> LexicalTM::CreateTM(const DataLattice & lattice) {
  VectorFst<LogArc> tm;
  VectorFst<LogArc>::StateId home = tm.AddState();
  tm.SetStart(home);
  tm.SetFinal(home, LogArc::Weight::One());

  // Adding x:e arcs
  for(auto e : lattice.GetTranslation()) {
    // Get the total counts of e
    int e_total = e_count_[e];

    // Add the <unk>:e arc.
    // Determine the probability
    LogWeight prob = Divide(log_alpha_, Plus(log_alpha_, LogWeight(-log(e_total))));
    /* Weight shouldn't be 1, but related to the base dist.*/
    tm.AddArc(home, LogArc(f_vocab_.GetId("<unk>"), e, prob, home));

    // Add the f:e arcs.
    for(auto it = f_count_.begin(); it != f_count_.end(); it++) {
      if (it->second > 0) {
        prob = Divide(LogWeight(-log(align_count_[{it->first,e}])), Plus(log_alpha_, LogWeight(-log(e_total))));
        tm.AddArc(home, LogArc(it->first, e, prob, home));
      }
    }
  }

  return tm;
}

// Takes a vector of phonemes and returns a string representation.
std::string LexicalTM::PhonemeWord(vector<WordId> phonemes) {
  ostringstream word_stream;
  word_stream << "w(";
  if(phonemes.size() > 0) {
    word_stream << f_vocab_.GetSym(phonemes[0]);
  }
  for(int i = 1; i < phonemes.size(); i++) {
    word_stream << "+" << f_vocab_.GetSym(phonemes[i]);
  }
  word_stream << ")";
  return word_stream.str();
}

/* Takes a phoneme--word Alignment and converts it to a phoneme-word--word alignment */
/*
Alignment LexicalTM::PhonemeWordAlignment(const Alignment & ph_alignment) {
  Alignment word_alignment;
  vector<WordId> && ph_buf = vector<WordId>();
  for(auto arrow : ph_alignment) {
    if(arrow.first == f_vocab_.GetId("<eps>")) {
      word_alignment.push_back({f_vocab_.GetId(PhonemeWord(ph_buf)), arrow.second});
      ph_buf = vector<WordId>();
    } else {
      ph_buf.push_back(arrow.first);
    }
  }
  return word_alignment;
}
*/

// Add a newly sampled sequence to the lexicon
void LexicalTM::AddWord(VectorFst<LogArc> & lexicon, vector<WordId> phonemes, std::string phoneme_word) {

  // Adding the word to the lexicon
  VectorFst<LogArc>::StateId home = lexicon.Start();
  VectorFst<LogArc>::StateId cur = home;
  ostringstream word_stream;
  vector<VectorFst<LogArc>::StateId> state_buf; /* Stores the lexicon StateIds associated with a word */
  if(phonemes.size() > 0) {
    VectorFst<LogArc>::StateId next = lexicon.AddState();
    state_buf.push_back(next);
    lexicon.AddArc(home, LogArc(phonemes[0], f_vocab_.GetId("<eps>"), LogWeight::One(), next));
    cur = next;
  }
  for(int i = 1; i < phonemes.size(); i++) {
    VectorFst<LogArc>::StateId next = lexicon.AddState();
    state_buf.push_back(next);
    lexicon.AddArc(cur, LogArc(phonemes[i], f_vocab_.GetId("<eps>"), LogWeight::One(), next));
    cur = next;
  }
  lexicon.AddArc(cur, LogArc(f_vocab_.GetId("<eps>"), f_vocab_.GetId(phoneme_word), LogWeight::One(), home));
  WriteSymbolSets();
  lexicon.Write("data/phoneme-prototyping/lexicon.fst");

  lexicon_states_[f_vocab_.GetId(phoneme_word)] = state_buf;
}

void LexicalTM::RemoveSample(const Alignment & align) {
  // The general idea is to decrement f_count_, e_count_, align_count_. If f_count_ is then empty, we want to remove that word from the lexicon.
  // Add words from the foreign side of the alignment to the lexicon and counts

  vector<WordId> && ph_buf = vector<WordId>();
  for(auto arrow : align) {
    if(arrow.first == f_vocab_.GetId("<eps>")) {
      // This means the current phoneme-word is finished and we can update the cache and lexicon

      // Construct a token that concatentates the phonemes.
      string phoneme_word = PhonemeWord(ph_buf);
      WordId ph_word_id = f_vocab_.GetId(phoneme_word);

      // Decrement the lexicon count.
      if(f_count_[ph_word_id] > 0) {
        f_count_[ph_word_id]--;
      }

/*
      if(f_count_[ph_word_id] == 0) {
        // Then we need to remove the path from the lexicon too.
        cout << "Deleting " << phoneme_word << " states: " << lexicon_states_[ph_word_id] << endl;
        lexicon_.DeleteStates(lexicon_states_[ph_word_id]);
      }
*/

      // Update the English cache.
      e_count_[arrow.second]--;
      assert(e_count_[arrow.second] >= 0);

      // Update the cache of alignments
      std::pair<WordId, WordId> word_arrow = {ph_word_id, arrow.second};
      align_count_[word_arrow]--;
      assert(align_count_[word_arrow] >= 0);

      // Reset ph_buf and state_buf for new word
      ph_buf = vector<WordId>();

    } else {
      // The phoneme word isn't over, so keep filling the buffer.
      ph_buf.push_back(arrow.first);
    }


  }

}

/* Getting phoneme alignments, updates lexicon and cache.*/
void LexicalTM::AddSample(const Alignment & align) {

  // Add words from the foreign side of the alignment to the lexicon and counts
  vector<WordId> && ph_buf = vector<WordId>();
  for(auto arrow : align) {
    if(arrow.first == f_vocab_.GetId("<eps>")) {
      // This means the current phoneme-word is finished and we can update the cache and lexicon

      // Construct a token that concatentates the phonemes.
      string phoneme_word = PhonemeWord(ph_buf);
      WordId ph_word_id = f_vocab_.GetId(phoneme_word);

      // If phoneme_word is in the lexicon, increment the count. Otherwise, put
      // it in the lexicon
      if(f_count_.count(ph_word_id) == 1 && f_count_[ph_word_id] > 0) {
        f_count_[ph_word_id]++;
      } else {
        // Then add the word to the lexicon update the foreign cache
        AddWord(lexicon_, ph_buf, phoneme_word);
        f_count_[ph_word_id] = 1;
      }

      // Update the English cache.
      if(e_count_.count(arrow.second) == 1 && e_count_[arrow.second] > 0) {
        e_count_[arrow.second]++;
      } else {
        e_count_[arrow.second] = 1;
      }

      // Update the cache of alignments
      std::pair<WordId, WordId> word_arrow = {ph_word_id, arrow.second};
      if(align_count_.count(word_arrow) == 1 && align_count_[word_arrow] > 0) {
        align_count_[word_arrow]++;
      } else {
        cout << word_arrow << " " << f_vocab_.GetSym(ph_word_id) << endl;
        align_count_[word_arrow] = 1;
      }

      // Reset ph_buf and state_buf for new word
      ph_buf = vector<WordId>();

    } else {
      // The phoneme word isn't over, so keep filling the buffer.
      ph_buf.push_back(arrow.first);
    }
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

// TODO Rename this function once we start clearing out the old code.
/* Yields P(f|e) using a Dirichlet process over observed foreign tokens.*/
/*
LogWeight LexicalTM::DirichletProbNew(WordId e, WordId f) {
  // Get the total counts of e
  int e_total = 0;
  for(auto it = count_map_.begin(); it != count_map_.end(); it++) {
    if(it->first.first == e) {
      e_total += it->second;
    }
  }

  LogWeight numer = fst::Plus(fst::Times(log_alpha_,1.0/foreign_count_map_.size()), LogWeight(-log(count_map_[{f,e}])));
  LogWeight denom = fst::Plus(log_alpha_, LogWeight(-log(e_total)));
  return fst::Divide(numer, denom);
}
*/

LogWeight LexicalTM::DirichletProb(int e, int f) {
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

void LexicalTM::WriteSymbolSets() {
  f_vocab_.Write("data/phoneme-prototyping/f_symbols.txt");
  e_vocab_.Write("data/phoneme-prototyping/e_symbols.txt");
}

Alignment LexicalTM::CreateSample(const DataLattice & lattice, LLStats & stats) {

  WriteSymbolSets();

  lattice.GetFst().Write("data/phoneme-prototyping/lattice.fst");
  cout << "Wrote lattice." << endl;

  ArcSort(&lexicon_, ILabelCompare<LogArc>());

  // Create a translation model that constrains its options to the translation of the lattice.
  VectorFst<LogArc> tm = CreateTM(lattice);
  tm.Write("data/phoneme-prototyping/tm.fst");
  cout << "Created TM" << endl;

  WriteSymbolSets();

  // Compose the lattice with the lexicon.
  ComposeFst<LogArc> latlex(lattice.GetFst(), lexicon_);
  cout << "Composed lattice with lexicon..." << endl;
  VectorFst<LogArc> veclatlex(latlex);
  veclatlex.Write("data/phoneme-prototyping/latlex.fst");

  ArcSort(&tm, ILabelCompare<LogArc>());

  ComposeFst<LogArc> composed_fst(latlex, tm);
  cout << "Composed latlex with tm..." << endl;

  VectorFst<LogArc> vecfst(composed_fst);
  vecfst.Write("data/phoneme-prototyping/composed.fst");

  // Sample from the composed Fst.
  VectorFst<LogArc> sample_fst;
  SampGen(composed_fst, sample_fst);
  sample_fst.Write("data/phoneme-prototyping/sample.fst");
  Alignment alignment = FstToAlign(sample_fst);
  PrintAlign(alignment);

  /*
  Alignment unk_alignment = PhonemeWordAlignment(align);
  PrintAlign(unk_alignment);
  */

  return alignment;

}

void LexicalTM::PrintAlign(const Alignment & alignment) {
  cout << "[";
  for(auto arrow : alignment) {
    if(f_vocab_.GetSym(arrow.first) == "<eps>") {
      cout << "<eps";
    } else {
      cout << "<" << f_vocab_.GetSym(arrow.first);
    }
    if(e_vocab_.GetSym(arrow.second) == "<eps>") {
      cout << ",eps>,";
    } else if(e_vocab_.GetSym(arrow.second) == "<unk>") {
      cout << ",unk>,";
    } else {
      cout << "," << e_vocab_.GetSym(arrow.second) << ">,";
    }
  }
  cout << "]" << endl;
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

/** Finds the best path through the lattices using the lexicon and translation model in this LexicalTM object. Therefore there is no averaging over a number of samples, the best path is just found using final cache counts.**/
void LexicalTM::FindBestPaths(
    const vector<DataLatticePtr> & lattices,
    const string out_fn) {

  ofstream && of = ofstream();
  of.open(out_fn);

  for(auto latticep : lattices) {
    DataLattice lattice = *latticep;
    VectorFst<LogArc> tm = CreateTM(lattice);

    ArcSort(&lexicon_, ILabelCompare<LogArc>());
    ComposeFst<LogArc> complatlex(lattice.GetFst(), lexicon_);
    VectorFst<LogArc> latlex(complatlex);

    ArcSort(&tm, ILabelCompare<LogArc>());
    ComposeFst<LogArc> compfst(latlex, tm);
    VectorFst<LogArc> composed(compfst);

    // Find the shortest path.
    VectorFst<LogArc> * shortest_path = new VectorFst<LogArc>;
    vector<int> && prev_state = vector<int>();
    vector<pair<int,int>> && prev_align = vector<pair<int,int>>();
    DataLattice::Dijkstra(composed, prev_state, prev_align);
    int final_state = DataLattice::GetFinal(composed);
    DataLattice::StringFromBacktrace(
        final_state, prev_state, prev_align, f_vocab_, of);
  }

  of.close();
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
      DataLattice::Dijkstra(vecfst, prev_state, prev_align);

      int final_state = DataLattice::GetFinal(composed_fst);
      DataLattice::StringFromBacktrace(final_state, prev_state, prev_align, f_vocab_, cout);
      DataLattice::AlignmentFromBacktrace(final_state, prev_state, prev_align, f_vocab_, e_vocab_, align_file);
  }
  align_file.close();
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

void LexicalTM::WriteSortedCounts() {
  std::ofstream f;
  f.open("align_counts.txt");

  vector<pair<pair<WordId,WordId>, int>> items;
  for(auto it = align_count_.begin(); it != align_count_.end(); it++) {
    items.push_back(*it);
  }

  /*
  vector<pair<pair<WordId,WordId>, float>> cond_items;
  for(auto item : items) {
    cond_items.push_back({item.first, float(item.second)/e_count_[item.first.second]});
  }
  */

  // Sort items by second element
  std::sort(items.begin(), items.end(),
            boost::bind(&std::pair<pair<WordId,WordId>, int>::second, _1) <
            boost::bind(&std::pair<pair<WordId,WordId>, int>::second, _2));

  /*
  std::sort(cond_items.begin(), cond_items.end(),
            boost::bind(&std::pair<pair<WordId,WordId>, float>::second, _1) <
            boost::bind(&std::pair<pair<WordId,WordId>, float>::second, _2));
  */


  for(auto item : items) {
    if(item.first.second > 0) {
      f << f_vocab_.GetSym(item.first.first) << " " << e_vocab_.GetSym(item.first.second) << " " << item.second << std::endl;
    }
  }

  f.close();

}
