#include <iostream>
#include <fst/fstlib.h>
#include <latticelm/sampgen.h>

using namespace latticelm;
using namespace fst;

int main(int argc, char *argv[])
{
  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << "input-LogFst\n";
  }

  // Load the model
  VectorFst<LogArc> *model = VectorFst<LogArc>::Read(argv[1]);

  //srand(0);
  // Sample a path
  int count2 = 0;
  int count1 = 0;
  for (int i = 0; i < 100000; i++) {
    VectorFst<LogArc> sample_fst;
    SampGen(*model, sample_fst);
    Sentence sent = FstToSent(sample_fst);
    cout << sent << "\n";
    cout << sent.size() << "\n";
    //if (sent.size() == 2) {
    //  count2++;
    //} else if (sent.size() == 1) {
    //  count1++;
   // }
    //sample_fst.Write("output.fst");
  }
  //cout << "1: " << count1 << "\n";
  //cout << "2: " << count2 << "\n";
}
