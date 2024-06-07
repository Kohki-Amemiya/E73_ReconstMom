// ====================================================================
//   EventAction.cc
//
// ====================================================================
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"

#include "EventAction.hh"
#include "RCSD.hh"

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"

#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4ThreeVector.hh"
#include "TRandom3.h"

#include "RunAction.hh"

int EventAction::nevt = 0;

//////////////////////////
EventAction::EventAction()
//////////////////////////
{
}

///////////////////////////
EventAction::~EventAction()
///////////////////////////
{
}

////////////////////////////////////////////////////////////
void EventAction::BeginOfEventAction(const G4Event *anEvent)
////////////////////////////////////////////////////////////
{
}

//////////////////////////////////////////////////////////
void EventAction::EndOfEventAction(const G4Event *anEvent)
//////////////////////////////////////////////////////////
{

  G4SDManager *SDManager = G4SDManager::GetSDMpointer();

  // get "Hit Collection of This Event"
  G4HCofThisEvent *HCTE = anEvent->GetHCofThisEvent();
  if (!HCTE)
    return;

  double beam_mass = anEvent->GetPrimaryVertex()->GetPrimary()->GetMass();
  double beam_mom = anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum().mag();
  G4double beam_posx = anEvent->GetPrimaryVertex()->GetPosition().x();
  G4double beam_posy = anEvent->GetPrimaryVertex()->GetPosition().y();
  G4double beam_posz = anEvent->GetPrimaryVertex()->GetPosition().z();
  // G4double beam_gtime = anEvent->GetPrimaryVertex()->GetPrimary()->GetT0(); //no GetTime() ??

  G4cout << "beam mass = " << beam_mass / CLHEP::MeV << G4endl;
  G4cout << "beam mom = " << beam_mom / CLHEP::MeV << G4endl;
  G4cout << "beam pos = (" << beam_posx << ", " << beam_posy << ", " << beam_posz << ")" << G4endl;
  G4double bg = 0;
  bg = beam_mom / beam_mass;
  G4double beta = bg / sqrt(1. + bg * bg);
  //  G4cout <<"beam beta = "<<bg / sqrt(1. + bg * bg)<<G4endl;
  G4cout << "EvID =  " << anEvent->GetEventID() << G4endl;

  // get a hit collection
  static G4int idrc = -1;
  if (idrc < 0)
    idrc = SDManager->GetCollectionID("rangecounter");
  RCHitsCollection *hcrc = (RCHitsCollection *)HCTE->GetHC(idrc);
  if (!hcrc)
  {
    G4cerr << "No Hit Collection" << G4endl;
    return; // no hit collection
  }
  // get hits
  G4int nhits = hcrc->entries();

  std::cout << "nhits=" << nhits << std::endl;

  G4int Maxch = 0;
  G4double dE = 0;
  G4double adc[18] = {-999};
  G4double mom[18] = {-999};
  G4double momx[18] = {-999};
  G4double momy[18] = {-999};
  G4double momz[18] = {-999};
  G4double gtime[18] = {-999};
  G4double resgtime[18] = {-999};
  G4double totE = 0.;
  bool hitflag[20]; //=false;
  for (int i = 0; i < 20; i++)
  {
    hitflag[i] = false;
  }

  for (G4int idx = 0; idx < nhits; idx++)
  {
    G4int ich = (*hcrc)[idx]->GetID();
    G4double edep = (*hcrc)[idx]->GetEdep();
    G4double mom_ = (*hcrc)[idx]->GetMom();
    G4double momx_ = (*hcrc)[idx]->GetMomx();
    G4double momy_ = (*hcrc)[idx]->GetMomy();
    G4double momz_ = (*hcrc)[idx]->GetMomz();
    G4double gtime_ = (*hcrc)[idx]->GetTime();
    G4double resgtime_ = gRandom->Gaus(gtime_, 0.1); // add 100ps resolution

    G4cout << "idx= " << idx << ", ich=" << ich << ", idrc:" << idrc << ", edep=" << edep << ", Mom=" << mom_ << ", resgtime=" << resgtime_ << G4endl;
    if (idx == 0)
    {
      dE = edep;
    }

    adc[ich] = edep;
    mom[ich] = mom_;
    momx[ich] = momx_;
    momy[ich] = momy_;
    momz[ich] = momz_;
    gtime[ich] = gtime_;
    resgtime[ich] = resgtime_;
    totE += edep;
  }

  TTree *tree_tmp = (TTree *)gROOT->FindObject("tree");
  if (nevt == 0)
  {
    tree_tmp->Branch("beam_mom", &beam_mom, "beam_mom/D");
    tree_tmp->Branch("beam_posx", &beam_posx, "beam_posx/D");
    tree_tmp->Branch("beam_posy", &beam_posy, "beam_posy/D");
    tree_tmp->Branch("beam_posz", &beam_posz, "beam_posz/D");
    tree_tmp->Branch("nhits", &nhits, "nhits/I");
    tree_tmp->Branch("adc", adc, "adc[18]/D");
    tree_tmp->Branch("mom", mom, "mom[18]/D");
    tree_tmp->Branch("momx", momx, "momx[18]/D");
    tree_tmp->Branch("momy", momy, "momy[18]/D");
    tree_tmp->Branch("momz", momz, "momz[18]/D");
    tree_tmp->Branch("totE", &totE, "totE/D");
    tree_tmp->Branch("gtime", gtime, "gtime[18]/D");
    tree_tmp->Branch("resgtime", &resgtime, "resgtime[18]/D");
  }

  tree_tmp->Fill();
  nevt++;
}
