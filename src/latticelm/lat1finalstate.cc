#include <iostream>
#include <fst/fstlib.h>

using namespace fst;

int main(int argc, char *argv[])
{
  if (argc != 3) {
    std::cout << "Usage: " << argv[0] << "input-LogFst output-LogFst\n";
  }

  VectorFst<LogArc> *model = VectorFst<LogArc>::Read(argv[1]);

  // Create a new final state.
  LogArc::StateId final_state_id = model->AddState();
  // Set One True Final state to Zero for now, will change to One later.
  model->SetFinal(final_state_id, LogWeight::Zero());

  // Change the weights for all current final states, and add free epsilon
  // transitions to the One True Final state.
  for (StateIterator<VectorFst<LogArc>> siter(*model); !siter.Done(); siter.Next()) {
    LogArc::StateId state_id = siter.Value();

    // If the state is final, then use the weight in an epsilon transition to
    // the final state, and change the current state's weight to zero.
    if (model->Final(state_id) != LogWeight::Zero()) {
      LogWeight weight = model->Final(state_id);
      // Add an arc to the new final state.
      model->AddArc(state_id, LogArc(0, 0, weight, final_state_id));
      // Make the state non-final.
      model->SetFinal(state_id, LogWeight::Zero());
    }
  }

  model->SetFinal(final_state_id, LogWeight::One());

  model->Write(argv[2]);
}
