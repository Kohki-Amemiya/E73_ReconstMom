// ====================================================================
//   RunAction.cc
//
// ====================================================================
#include "G4Run.hh"
#include "RunAction.hh"
#include "G4Timer.hh"
#include "globals.hh"
#include "TROOT.h"
#include "TRint.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"

//////////////////////
RunAction::RunAction()
//////////////////////
{

  tree = new TTree("tree", "tree");
  //  tree->Branch("beta1",&beta1,"beta1/S");
}

///////////////////////
RunAction::~RunAction()
///////////////////////
{
}

///////////////////////////////////////////////////
void RunAction::BeginOfRunAction(const G4Run *aRun)
///////////////////////////////////////////////////
{
}

/////////////////////////////////////////////////
void RunAction::EndOfRunAction(const G4Run *aRun)
/////////////////////////////////////////////////
{
  TFile *file = new TFile("rootfiles/E73simulator_run031.root",
                          "RECREATE", "Geant4 ROOT analysis");

  tree->Write();
  file->Close();
  delete file;
}
