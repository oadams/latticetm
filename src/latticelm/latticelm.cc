#include <iostream>
#include <fstream>
#include <string>
#include <numeric>
#include <boost/program_options.hpp>
#include <latticelm/timer.h>
#include <latticelm/latticelm.h>
#include <latticelm/data-lattice.h>
#include <latticelm/hierarchical-lm.h>
#include <latticelm/lexical-tm.h>
#include <latticelm/ll-stats.h>
#include <latticelm/macros.h>

#include <fst/compose.h>
#include <latticelm/sampgen.h>
#include <sstream>
#include <unordered_map>

using namespace std;
namespace po = boost::program_options;

namespace latticelm {

// TODO Perhaps this method should really be in LexicalTM.
void LatticeLM::PerformTrainingLexTM(const vector<DataLatticePtr> & all_lattices, LexicalTM & tm, int train_len, int test_len) {

  assert(train_len > 0);
  assert(test_len > 0);

  vector<DataLatticePtr>::const_iterator train_start;
  train_start = all_lattices.begin() + (all_lattices.size() - train_len);
  vector<DataLatticePtr> train_lattices(train_start, all_lattices.end());
  assert(train_lattices.size() == train_len);

  vector<DataLatticePtr>::const_iterator test_start;
  test_start = all_lattices.begin() + (all_lattices.size() - test_len);
  vector<DataLatticePtr> test_lattices(test_start, all_lattices.end());
  assert(test_lattices.size() == test_len);

  // Perform training
  vector<int> order(train_lattices.size()); std::iota(order.begin(), order.end(), 0);
  vector<Alignment> alignments(train_lattices.size());
  for(int epoch = 1; epoch <= epochs_; epoch++) {
    std::shuffle(order.begin(), order.end(), *GlobalVars::rndeng);
    LLStats ep_stats;
    int align_count = 0;
    for(int align_id : order) {
      cerr << endl;
      cerr << "align " << align_count << ", align_id: " << align_id << endl;
      cerr << "time: " << time_.Elapsed() << endl;
      align_count++;
      if(epoch != 1)
        tm.RemoveSample(alignments[align_id]);
      alignments[align_id] = tm.CreateSample(*train_lattices[align_id], ep_stats);

      if(align_count % 100 == 0) {
        tm.WriteSortedCounts();
      }

      tm.AddSample(alignments[align_id]);
    }
    cerr << "Finished epoch " << epoch << ": char=" << ep_stats.words_ << ", ppl=" << ep_stats.CalcPPL() << " (s=" << time_.Elapsed() << ")" << endl;
    //tm.ResampleParameters();
    //tm.PrintParams("data/out/params/tm.sample" + to_string(epoch));
  }
  //tm.Normalize(epochs_);
  //tm.PrintParams("data/out/params/tm.avg");
  //tm.FindBestPaths(test_lattices, "data/out/alignments.txt");
}

template <class LM>
void LatticeLM::PerformTraining(const vector<DataLatticePtr> & lattices, LM & lm) {

  // Perform training
  vector<int> order(lattices.size()); std::iota(order.begin(), order.end(), 0);
  vector<Sentence> sentences(lattices.size());
  for(int epoch = 1; epoch <= epochs_; epoch++) {
    std::shuffle(order.begin(), order.end(), *GlobalVars::rndeng);
    LLStats ep_stats;
    for(int sid : order) {
      if(epoch != 1)
        lm.RemoveSample(sentences[sid]);
      sentences[sid] = lm.CreateSample(*lattices[sid], ep_stats);
      lm.AddSample(sentences[sid]);
    }
    cerr << "Finished epoch " << epoch << ": char=" << ep_stats.words_ << ", ppl=" << ep_stats.CalcPPL() << " (s=" << time_.Elapsed() << ")" << endl;
    lm.ResampleParameters();
  }
}

int LatticeLM::main(int argc, char** argv) {
  po::options_description desc("*** latticelm (by Graham Neubig) ***");
  desc.add_options()
      ("help", "Produce help message")
      ("train_file", po::value<string>()->default_value(""), "Training file")
      ("train_ref", po::value<string>()->default_value(""), "Training reference file containing true phoneme strings (optional)")
      ("trans_file", po::value<string>()->default_value(""), "File containing word-tokenized translations of the training lattices in plain text.")
      ("file_format", po::value<string>()->default_value("text"), "The format of the lattices in the input file")
      ("model_type", po::value<string>()->default_value("pylm"), "Model type (hierlm to do segmentation and LM learning, pylm to just do lm learning)")
      ("beam", po::value<int>()->default_value(0), "Beam size")
      ("epochs", po::value<int>()->default_value(100), "Epochs")
      ("word_n", po::value<int>()->default_value(2), "Length of word n-grams")
      ("char_n", po::value<int>()->default_value(2), "Length of character n-grams")
      ("model_in", po::value<string>()->default_value(""), "The file to read the model to")
      ("model_out", po::value<string>()->default_value(""), "The file to write the final model to")
      ("seed", po::value<int>()->default_value(0), "The random seed, or 0 to change every time")
      ("lattice_weight", po::value<float>()->default_value(1.f), "Amount of weight to give to the lattice probabilities")
      ("verbose", po::value<int>()->default_value(1), "Verbosity of messages to print")
      ("concentration", po::value<float>()->default_value(1.0), "The concentration parameter for the Dirichlet process of the translation model.")
      ("plain_best_paths", po::value<string>()->default_value(""), "Just output the 1-best path through the supplied lattice.")
      ("train_len", po::value<int>()->default_value(-1), "Number of training sents")
      ("test_len", po::value<int>()->default_value(-1), "Number of test sents")
      ("using_external_tm", po::value<string>()->default_value(""), "For using an external TM to perform decoding")
      ;
  boost::program_options::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);   
  if (vm.count("help")) {
      cout << desc << endl;
      return 1;
  }

  // Temporary buffers
  string line;

  // Save various settings
  epochs_ = vm["epochs"].as<int>();
  beam_ = vm["beam"].as<int>();
  char_n_ = vm["char_n"].as<int>();
  word_n_ = vm["word_n"].as<int>();
  lattice_weight_ = vm["lattice_weight"].as<float>();
  file_format_ = vm["file_format"].as<string>();
  model_type_ = vm["model_type"].as<string>();
  alpha_ = vm["concentration"].as<float>();

  GlobalVars::Init(vm["verbose"].as<int>(), vm["seed"].as<int>());

  // Initialize the vocabulary
  cids_.GetId("<eps>");
  cids_.GetId("<unk>");
  //cids_.GetId("<s>");
  //cids_.GetId("</s>");

  // Initialize the translation vocabulary
  trans_ids_.GetId("<eps>");
  //trans_ids_.GetId("<s>");
  //trans_ids_.GetId("</s>");

  unordered_set<std::string> phonemes;

  // Load data
  vector<DataLatticePtr> lattices = DataLattice::ReadFromFile(file_format_, lattice_weight_, vm["train_file"].as<string>(), vm["trans_file"].as<string>(), cids_, trans_ids_, phonemes);

  unordered_map<pair<WordId,WordId>, int> map;
  pair<WordId,WordId> x (2,7);
  map.insert({x, 5});

  cout << map[x] << endl;

  //Prototyping(lattices);
  float gamma = 0.9;

  if(!vm["plain_best_paths"].as<string>().empty()) {
    DataLattice::FindBestPaths(
        lattices,
        "data/out/" + vm["plain_best_paths"].as<string>(),
        cids_);
    return 0;
  }

  if(!vm["using_external_tm"].as<string>().empty()) {
    LexicalTM tm(cids_, trans_ids_, alpha_, gamma, phonemes);
    vector<vector<fst::LogWeight>> tm_params = tm.load_TM(vm["using_external_tm"].as<string>());
    tm.FindBestPaths(lattices, "data/out/external_tm_alignments.txt", tm_params);
    return 0;
  }

  // Create the timer
  time_ = Timer();
  cerr << "Started training! (s=" << time_.Elapsed() << ")" << endl;

  cids_.Write("data/phoneme-prototyping/f_vocab_.txt");

  // Create the hierarchical LM
  if(model_type_ == "pylm") {
    Pylm pylm(cids_.size(), char_n_);
    PerformTraining(lattices, pylm);
  } else if(model_type_ == "hierlm") {
    HierarchicalLM hlm(cids_.size(), char_n_, word_n_);
    PerformTraining(lattices, hlm);
  } else if(model_type_ == "lextm") {
    LexicalTM tm(cids_, trans_ids_, alpha_, gamma, phonemes);
    tm.WriteSortedCounts();
    PerformTrainingLexTM(lattices, tm, vm["train_len"].as<int>(), vm["test_len"].as<int>());
  }

  return 0;

}

}
