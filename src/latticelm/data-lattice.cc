#include <latticelm/data-lattice.h>
#include <latticelm/macros.h>
#include <latticelm/sentence.h>
#include <boost/algorithm/string.hpp>
#include <fst/script/compile-impl.h>
#include <iostream>
#include <fstream>

using namespace latticelm;
using namespace std;
using namespace fst;

vector<DataLatticePtr> DataLattice::ReadFromFile(const std::string & format, float weight, const std::string & filename, const std::string & trans_filename, SymbolSet<string> & dict, SymbolSet<string> & trans_dict) {
  vector<DataLatticePtr> data_lattices;
  if(format == "text") {
    data_lattices = ReadFromTextFile(filename, weight, dict);
  } else if (format == "openfst") {
    data_lattices = ReadFromOpenFSTFile(filename, weight, dict);
  } else {
    THROW_ERROR("Illegal file format: " << format);
  }
  if(!trans_filename.empty()) {
    DataLattice::ReadTranslations(data_lattices, trans_filename, trans_dict);
  }
  return data_lattices;
}

vector<DataLatticePtr> DataLattice::ReadFromTextFile(const std::string & filename, float weight, SymbolSet<string> & dict) {
  string line;
  ifstream in(filename);
  if(!in) THROW_ERROR("Could not open " << filename);
  vector<DataLatticePtr> ret;
  while(getline(in, line)) {
    Sentence sent = ParseSentence(line, dict);
    if(*sent.rbegin() != 2) sent.push_back(2); // All sentences must end with a sentence ending
    DataLatticePtr ptr(new DataLattice);
    VectorFst<LogArc>::StateId last_id = ptr->fst_.AddState(), next_id;
    ptr->fst_.SetStart(last_id);
    for(auto wid : sent) {
      next_id = ptr->fst_.AddState();
      ptr->fst_.AddArc(last_id, LogArc(wid, wid, 0.f, next_id));
      ptr->f_wordids_.insert(wid);
      last_id = next_id;
    }
    ptr->fst_.SetFinal(last_id, LogArc::Weight::One());
    ret.push_back(ptr);
  }
  return ret;
}

vector<DataLatticePtr> DataLattice::ReadFromOpenFSTFile(const std::string & filename, float weight, SymbolSet<string> & dict) {
  /** Be wary of the assumptions this method makes:
    *   - The input file includes some number of FSTs, each of which is separated by a blank line.
    *   - An FST description is comprised of a number of lines. Each line
    *   represents an arc and is delimited by tabs or spaces. The first two values are
    *   the ids from and to. The second to are the transduction input and
    *   output. The final value is the weight.
    *   - The first state in the first listed arc of each FST is the sole start state.
    *   - The final state in the final listed arc of each FST is the sole final state.
    *   - StateIds created by Add state start from 0 and increment.
    *   - I assume I can implicitly cast or convert an integer to a stateid (as evidenced by the call to stoi())
    */
  string line;
  ifstream in(filename);
  if(!in) THROW_ERROR("Could not open " << filename);
  vector<DataLatticePtr> ret;
  // Initialize lattice
  DataLatticePtr ptr(new DataLattice);
  VectorFst<LogArc>::StateId last_id = ptr->fst_.AddState();
  ptr->fst_.SetStart(last_id);
  VectorFst<LogArc>::StateId num_states = last_id + 1;
  VectorFst<LogArc>::StateId to_state;
  int lineid = 0;
  while(getline(in, line)) {
    lineid++;
    if(line == "") {
      // If there are no more lines after this, let's leave this loop.
      if(!getline(in, line)) {
        break;
      }
      lineid++;
      // Otherwise wrap up this lattice and initialize a new one.
      ptr->fst_.SetFinal(to_state, LogArc::Weight::One());
      //DataLattice::Dijkstra(ptr->fst_, dict);
      ret.push_back(ptr);
      ptr = DataLatticePtr(new DataLattice);
      VectorFst<LogArc>::StateId last_id = ptr->fst_.AddState();
      ptr->fst_.SetStart(last_id);
      num_states = last_id + 1;
    }
    // Read in tokens
    vector<string> line_tokens;
    boost::split(line_tokens, line, boost::is_any_of("\t "), boost::token_compress_on);
    if(line_tokens.size() != 5) {
      cerr << "line " << lineid << ": " << line_tokens << endl;
      THROW_ERROR("Ill-formed FST input. Each line must consist of 5 tokens tab or space delimited.")
    }
    VectorFst<LogArc>::StateId from_state = stoi(line_tokens[0]);
    to_state = stoi(line_tokens[1]);
    WordId in = dict.GetId(line_tokens[2]);
    WordId out = dict.GetId(line_tokens[3]);
    LogWeight arc_weight = LogWeight(stof(line_tokens[4])*weight);
    // Add any necessary states before we add the arc.
    while(num_states < from_state+1 || num_states < to_state+1) {
      ptr->fst_.AddState();
      num_states += 1;
    }
    ptr->fst_.AddArc(from_state, LogArc(in, out, arc_weight, to_state));
    ptr->f_wordids_.insert(in);
  }
  // Wrap up the last uncompleted lattice.
  ptr->fst_.SetFinal(to_state, LogArc::Weight::One());
  dict.Write("data/out/lattices/isymbols.txt");
  ret.push_back(ptr);
  return ret;
}

void DataLattice::ReadTranslations(vector<DataLatticePtr> data_lattices, const string & trans_filename, SymbolSet<string> & trans_dict) {
  // Assuming each line in the translation file corresponds to one lattice, we iterate
  // through the already loaded lattices and give them their corresponding translation.
  string line;
  ifstream in(trans_filename);
  if(!in) THROW_ERROR("Could not open " << trans_filename);
  int i = 0;
  while(getline(in, line)) {
    // Tokenize the string using whitespace.
    Sentence sent = ParseSentence(line, trans_dict);
    data_lattices[i++]->SetTranslation(sent);
  }
  if(i != data_lattices.size()) THROW_ERROR("Number of lattices and number of translations are not equal.");
  trans_dict.Write("data/out/lattices/osymbols.txt");
}

/** Uses Dijkstra's shortest path algorithm to find the shortest path through
 * the lattice. This will be used for finding the best source sentence and
 * alignment in the composed lattice.

 * I'm implementing this because the FST ShortestPath implementation:
 * http://www.openfst.org/twiki/bin/view/FST/ShortestPathDoc indicates that the
 * `path' property must hold for the weights. But this path property does not
 * hold for log weights it would seem.

 * Though I could convert Fst<LogArc>s to Fst<StdArc>s, I had pretty much
 * implemented this by the time I found out that would be equivalent.*/
void DataLattice::Dijkstra(const Fst<LogArc> & lattice,
                            vector<int> & prev_state, vector<pair<int,int>> & prev_align,
                            SymbolSet<string> & dict, SymbolSet<string> & trans_dict, bool debug) {
  VectorFst<LogArc>::StateId initial_state = lattice.Start();
  assert(initial_state == 0);
  //VectorFst<LogArc>::StateId final_state = lattice.NumStates()-1;

  vector<float> min_distance;
  min_distance.push_back(0.0);
  prev_state.push_back(-1);
  prev_align.push_back({-1,-1});
  set<pair<float,VectorFst<LogArc>::StateId>> active_vertices;
  active_vertices.insert( {0.0, initial_state} );

  std::ofstream debug_stream;
  if (debug) debug_stream.open("data/out/debug_stream.txt");

  if (debug) debug_stream << "active_vertices.begin()->first: " << active_vertices.begin()->first << std::endl;

  while(!active_vertices.empty()) {
    int cur = active_vertices.begin()->second;
    active_vertices.erase(active_vertices.begin());
    fst::ArcIterator<Fst<LogArc>> arc_iter(lattice, cur);
    if (debug) debug_stream << "Iterating over the arcs from state " << cur << endl;
    while(true) {
      if(arc_iter.Done()) break;
      const LogArc& arc = arc_iter.Value();
      //cout << arc.weight << " " << arc.ilabel << " " << arc.olabel << endl;
      // Expand min_distance if we need to.
      while(arc.nextstate+1 > min_distance.size()) {
        min_distance.push_back(LogWeight::Zero().Value());
        prev_state.push_back(-1);
        pair<int,int> nullpair = {-1,-1};
        prev_align.push_back(nullpair);
      }
      if (debug) debug_stream << "min_distance: " << min_distance << endl;
      if (debug) debug_stream << "\tDealing with arc: " << arc.ilabel << ":" << arc.olabel << "/" << arc.weight << " to " << arc.nextstate << endl;
      if (debug) debug_stream << "\t\tmin_distance[cur]: " << min_distance[cur] << endl;
      if (debug) debug_stream << "\t\tfst::Times(min_distance[cur],arc.weight): " << fst::Times(min_distance[cur],arc.weight).Value() << endl;
      if (debug) debug_stream << "\t\tmin_distance[arc.nextstate]: " << min_distance[arc.nextstate] << endl;
      if(fst::Times(min_distance[cur],arc.weight).Value() <= min_distance[arc.nextstate]) {
        active_vertices.erase( { min_distance[arc.nextstate], arc.nextstate } );
        min_distance[arc.nextstate] = fst::Times(min_distance[cur], arc.weight).Value();
        if (debug) debug_stream << "\t\tmin_distance[arc.nextstate]: " << min_distance[arc.nextstate] << endl;
        prev_state[arc.nextstate] = cur;
        prev_align[arc.nextstate] = {arc.ilabel, arc.olabel};
        active_vertices.insert( { min_distance[arc.nextstate], arc.nextstate } );
      }
      arc_iter.Next();
    }
  }
  if (debug) debug_stream << "prev_state: " << prev_state << endl;
  if (debug) debug_stream << "prev_align: " << prev_align << endl;
  debug_stream.close();
}

void DataLattice::StringFromBacktrace(const int final_state_id, const vector<int> & prev_state, const vector<pair<int,int>> & prev_align, SymbolSet<string> & dict, ostream & out_stream) {
  int id = final_state_id;
  vector<string> foreign_source;
  while(true) {
    int wordid = prev_align[id].first;
    if(wordid == -1) break;
    foreign_source.push_back(dict.GetSym(wordid));
    id = prev_state[id];
  }
  for(int i = foreign_source.size()-1; i >= 0; i--){
      out_stream << foreign_source[i] << " ";
  }
  out_stream << endl;
}

void DataLattice::AlignmentFromBacktrace(const VectorFst<LogArc>::StateId final_state_id, const vector<int> & prev_state, const vector<pair<int,int>> & prev_align, SymbolSet<string> & dict, SymbolSet<string> & trans_dict, ofstream & align_file) {
  int id = final_state_id;
  vector<pair<string,string>> alignments;
  while(true) {
    int f_wordid = prev_align[id].first;
    int e_wordid = prev_align[id].second;
    if(f_wordid == -1) break;
    alignments.push_back({dict.GetSym(f_wordid), trans_dict.GetSym(e_wordid)});
    id = prev_state[id];
  }
  for(int i = alignments.size()-1; i >=0; i--) {
    align_file << alignments[i] << " ";
  }
  align_file << endl;
}

/* Returns the phonemes present in the lattice */
vector<string> DataLattice::GetPhonemes(const vector<DataLatticePtr> & lattices) {
  // TODO Implement
  for(auto lattice_ptr : lattices) {
    DataLattice lattice = *lattice_ptr;
  }
  vector<string> phonemes = {"h", "aU", "s", "O", "f"};
  return phonemes;
}

