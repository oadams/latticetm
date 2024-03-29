#pragma once

#include <latticelm/sentence.h>
#include <latticelm/macros.h>
#include <fst/fst.h>
#include <fst/mutable-fst.h>
#include <vector>
#include <stdexcept>

using namespace std;

namespace latticelm {

inline float neglogsumexp(const vector<float> & probs) {
    float m = probs[0];
    for(int i = 1; i < (int)probs.size(); i++)
        m = fmin(m, probs[i]);
    float ret = 0;
    for(int i = 0; i < (int)probs.size(); i++)
        ret += exp(m-probs[i]);
    return log(ret)-m;
}  

// sample a single value appropriately from a vector of weights
unsigned SampleWeights(vector<float> & ws, float anneal = 1);

template<class A>
float SampGen(const fst::Fst<A> & ifst, fst::MutableFst<A> & ofst, unsigned nbest = 1, float anneal = 1) { 
  typedef fst::Fst<A> F;
  typedef typename F::Weight W;
  typedef typename A::StateId S;

  // sanity check
  if(ifst.Start() == -1)
    THROW_ERROR("Attempting to sample an empty FST");
  if(ifst.Final(ifst.Start()) != std::numeric_limits<float>::infinity())
    THROW_ERROR("Sampling FSTs where start states are final is not supported yet");

  // the number of remaining incoming arcs, and total weights of each state
  std::vector< int > incomingArcs;
  std::vector< vector< A > > backArcs;
  std::vector< W > stateWeights;
  unsigned i, statesFinished = 0;
  
  // intialize the data values
  for (fst::StateIterator< fst::Fst<A> > siter(ifst); !siter.Done(); siter.Next()) {
    S s = siter.Value();
    // cout << "state: " << s << endl;
    for(fst::ArcIterator< F > aiter(ifst, s); !aiter.Done(); aiter.Next()) {
      const A& a = aiter.Value();
      // cout << " -> " << a.nextstate << endl;
      while((unsigned)a.nextstate >= incomingArcs.size()) {
        incomingArcs.push_back(0);
        stateWeights.push_back(W::Zero());
        backArcs.push_back(std::vector<A>());
      }
      incomingArcs[a.nextstate]++;
      backArcs[a.nextstate].push_back(A(a.ilabel, a.olabel, a.weight, s));
    }
  }
  stateWeights[ifst.Start()] = W::One();
  incomingArcs[ifst.Start()] = 0;

  // calculate the number of arcs incoming to each state
  vector< S > stateQueue(1,ifst.Start());
  while(stateQueue.size() > 0) {
    unsigned s = stateQueue[stateQueue.size()-1];
    stateQueue.pop_back();
    for(fst::ArcIterator< F > aiter(ifst, s); !aiter.Done(); aiter.Next()) {
      const A& a = aiter.Value();
      // cout << "stateWeights[" << a.nextstate << "]: (" << stateWeights[a.nextstate]<<"+("<<stateWeights[s]<<"*"<<a.weight<<"))"<<endl;
      stateWeights[a.nextstate] = Plus(stateWeights[a.nextstate],Times(stateWeights[s],a.weight));
      // cout << " -> " << stateWeights[a.nextstate] << endl;
      if(--incomingArcs[a.nextstate] == 0)
        stateQueue.push_back(a.nextstate);
    }
    statesFinished++;
  }
  if(statesFinished != incomingArcs.size())
    throw std::runtime_error("Sampling cannot be performed on cyclic FSTs");


  // sample the states backwards from the final state
  ofst.AddState();
  ofst.SetStart(0);

  // find the final states and sample a final state
  vector< float > stateCandWeights;
  vector< S > stateCandIds;
  for (fst::StateIterator< fst::Fst<A> > siter(ifst); !siter.Done(); siter.Next()) {
    S s = siter.Value();
    // cerr << "stateWeights[" << s << "] == " << stateWeights[s].Value() << endl;
    float w = Times(ifst.Final(s),stateWeights[s]).Value();
    if(w != std::numeric_limits<float>::infinity()) {
      // cout << "Final state "<<s<<","<<w<<endl;
      stateCandWeights.push_back( w );
      stateCandIds.push_back( s );
    }
  }
  float ret = neglogsumexp(stateCandWeights);

  for(unsigned n = 0; n < nbest; n++) {

    S currState = stateCandIds[SampleWeights(stateCandWeights, anneal)];

    // add the final state
    S outState = (ifst.Start() != currState?ofst.AddState():0);
    ofst.SetFinal(outState, ifst.Final(currState));

    // sample the values in order
    while(outState != 0) {
      const vector<A> & arcs = backArcs[currState];
      vector<float> arcWeights(arcs.size(), 0);
      for(i = 0; i < arcs.size(); i++) 
        arcWeights[i] = Times(arcs[i].weight,stateWeights[arcs[i].nextstate]).Value();
      const A & myArc = arcs[SampleWeights(arcWeights, anneal)];
      S nextOutState = (myArc.nextstate != ifst.Start()?ofst.AddState():0);
      // cout << "Adding arc " << nextOutState << "--"<<myArc.ilabel<<"/"<<myArc.olabel<<":"<<myArc.weight<<"-->"<<outState<<endl;
      ofst.AddArc(nextOutState, A(myArc.ilabel,myArc.olabel,myArc.weight,outState));
      outState = nextOutState;
      currState = myArc.nextstate;
    }
  
  }

  return ret;

}

template<class A>
Sentence FstToSent(const fst::Fst<A> & ifst) {
  typedef fst::Fst<A> F;
  auto sid = ifst.Start();
  Sentence ret;
  while(true) {
    typename fst::ArcIterator< F > aiter(ifst, sid);
    if(aiter.Done()) break;
    const A& a = aiter.Value();
    if(a.olabel != 0) ret.push_back(a.olabel);
    sid = a.nextstate;
  }
  return ret;
}

template<class A>
Alignment FstToAlign(const fst::Fst<A> & ifst) {
  typedef fst::Fst<A> F;
  auto state_id = ifst.Start();
  Alignment ret;
  while(true) {
    typename fst::ArcIterator<F> aiter(ifst, state_id);
    if(aiter.Done()) break;
    const A& arc = aiter.Value();
    if(arc.olabel != 0 || arc.ilabel != 0) {
      ret.push_back(std::pair<WordId,WordId>(arc.ilabel, arc.olabel));
    }
    state_id = arc.nextstate;
  }
  return ret;
}

}
