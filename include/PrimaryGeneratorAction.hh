// ====================================================================
//    PrimaryGeneratorAction.hh
//
// ====================================================================
#ifndef PRIMARY_GENERATOR_ACTION_H
#define PRIMARY_GENERATOR_ACTION_H
 
#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;

// ====================================================================
//
// class definition
//
// ====================================================================

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
private:
  // use G4 particle gun
  G4ParticleGun* particleGun;

public:
  PrimaryGeneratorAction();
  ~PrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event* anEvent);

};

#endif
