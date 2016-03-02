#pragma once

#include <vector>
#include <set>
#include <memory>
#include <sstream>
#include <iterator>
#include <fst/vector-fst.h>
#include <latticelm/symbol-set.h>
#include <latticelm/sentence.h>

using namespace fst;

namespace latticelm {

class DataLattice;
typedef std::shared_ptr<DataLattice> DataLatticePtr;

class DataLattice {

public:
  DataLattice() { }
  ~DataLattice() { }
 
  static std::vector<DataLatticePtr> ReadFromFile(const std::string & format, float weight, const std::string & filename, const std::string & trans_filename, SymbolSet<std::string> & dict, SymbolSet<std::string> & trans_dict);
  static std::vector<DataLatticePtr> ReadFromTextFile(const std::string & filename, float weight, SymbolSet<std::string> & dict);
  static std::vector<DataLatticePtr> ReadFromOpenFSTFile(const std::string & filename, float weight, SymbolSet<std::string> & dict);
  static void ReadTranslations(vector<DataLatticePtr> data_lattices, const string & trans_filename, SymbolSet<std::string> & trans_dict);

  const VectorFst<LogArc> & GetFst() const { return fst_; }

  const Sentence GetTranslation() const {
    return translation_;
  }
  void SetTranslation(Sentence translation) {
    translation_ = translation;
  }

  const set<WordId> GetFWordIds() const {
    return f_wordids_;
  }

  static void Dijkstra(const fst::Fst<fst::LogArc> & lattice, SymbolSet<string> & dict, SymbolSet<string> & trans_dict);
  static void StringFromBacktrace(const vector<int> & prev_state, const vector<pair<int,int>> & prev_align, SymbolSet<string> & dict);
  static void AlignmentFromBacktrace(const vector<int> & prev_state, const vector<pair<int,int>> & prev_align, SymbolSet<string> & dict, SymbolSet<string> & trans_dict);

protected:
  fst::VectorFst<LogArc> fst_;
  // A word-tokenized English translation for building translation models.
  Sentence translation_;
  // A set of WordIds on the foreign side that indicate foreign tokens that
  // occur in this lattice. This is used to optimize the creation of reduced
  // TMs for composition so that they don't have any superfluous arcs.
  set<WordId> f_wordids_;

};

}
