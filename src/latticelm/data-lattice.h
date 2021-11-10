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
 
  static std::vector<DataLatticePtr> ReadFromFile(const std::string & format, float weight, const std::string & filename, const std::string & trans_filename, SymbolSet<std::string> & dict, SymbolSet<std::string> & trans_dict, std::unordered_set<std::string> & phonemes);
  static std::vector<DataLatticePtr> ReadFromTextFile(const std::string & filename, float weight, SymbolSet<std::string> & dict);
  static std::vector<DataLatticePtr> ReadFromOpenFSTFile(const std::string & filename, float weight, SymbolSet<std::string> & dict, std::unordered_set<std::string> & phonemes);
  static void ReadTranslations(std::vector<DataLatticePtr> data_lattices, const std::string & trans_filename, SymbolSet<std::string> & trans_dict);

  const VectorFst<LogArc> & GetFst() const { return fst_; }

  const Sentence GetTranslation() const {
    return translation_;
  }
  void SetTranslation(Sentence translation) {
    translation_ = translation;
  }

  const std::set<WordId> GetFWordIds() const {
    return f_wordids_;
  }

  static void Dijkstra(
      const fst::Fst<fst::LogArc> & lattice,
      std::vector<int> & prev_state,
      std::vector<std::pair<int,int>> & prev_align);
  static void StringFromBacktrace(
      const int final_state_id,
      const std::vector<int> & prev_state,
      const std::vector<std::pair<int,int>> & prev_align,
      SymbolSet<std::string> & dict,
      std::ostream & out_stream);
  static void AlignmentFromBacktrace(
      const int final_state_id,
      const std::vector<int> & prev_state,
      const std::vector<std::pair<int,int>> & prev_align,
      SymbolSet<std::string> & dict,
      SymbolSet<std::string> & trans_dict,
      std::ofstream & align_file);
  static void FindBestPaths(
      const std::vector<DataLatticePtr> & lattices,
      const std::string out_fn,
      SymbolSet<std::string> & dict);

  static int GetFinal(const fst::Fst<fst::LogArc> & fst) {
    for (StateIterator<Fst<LogArc>> iter(fst);
        !iter.Done();
        iter.Next()) {
      int state_id = iter.Value();
      if(fst.Final(state_id) == LogArc::Weight::One()) {
        return state_id;
      }
    }
    THROW_ERROR("No final state.");
  }

  int GetFinal() {
    return GetFinal(fst_);
  }

protected:
  fst::VectorFst<LogArc> fst_;
  // A word-tokenized English translation for building translation models.
  Sentence translation_;
  // A std::set of WordIds on the foreign side that indicate foreign tokens that
  // occur in this lattice. This is used to optimize the creation of reduced
  // TMs for composition so that they don't have any superfluous arcs.
  std::set<WordId> f_wordids_;

};

}
