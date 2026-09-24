#ifndef WCSimSteppingAction_h
#define WCSimSteppingAction_h 1


#include "G4Event.hh"
#include "G4UserSteppingAction.hh"
#include "G4ThreeVector.hh"

#include "WCSimRunAction.hh"

class G4HCofThisEvent;
class G4Event;

// First we have a structure to hold the values of the emitted photon.
typedef struct {

  G4ThreeVector direction; // The direction of the emitted photon.
  G4ThreeVector position; // The position where the emitted photon was created (same as where the incident photon was absorbed).
  double energy; // The energy value of the emitted photon.
  double wavelength; // The wavelength of the emitted photon.

} EmittedPhoton;


class WCSimSteppingAction : public G4UserSteppingAction
{
 private:
  WCSimRunAction* runAction;
  WCSimDetectorConstruction* det;

public:
  WCSimSteppingAction(WCSimRunAction*,WCSimDetectorConstruction*);

  ~WCSimSteppingAction()
  { };

  void UserSteppingAction(const G4Step*);

  G4int G4ThreeVectorToWireTime(G4ThreeVector *pos3d,
				G4ThreeVector lArPos,
				G4ThreeVector start,
				G4int i);
  
  void Distortion(G4double x,
		  G4double y);

  G4double FieldLines(G4double x,
		      G4double y,
		      G4int xy);

  static G4int n_photons_through_mPMTLV;
  static G4int n_photons_through_acrylic;
  static G4int n_photons_through_gel;
  static G4int n_photons_on_blacksheet;
  static G4int n_photons_on_smallPMT;
  static int primaryMuonTrackID; // ADDed by sahar
  static int primaryMuonEventID; // ADDed by sahar
  static int GetPrimaryMuonTrackID() { return primaryMuonTrackID; }  // ADDed by sahar

  WCSimRunAction* GetRunAction(){return runAction;}

private:

  G4double ret[2];

};
#include <unordered_map>
#include <string>
#include "G4Track.hh"

struct AllPhotonBirthInfo {
  int   evt    = -1;
  int   parent = -1;
  float ex_cm=0, ey_cm=0, ez_cm=0, et_ns=0;
  float edir_x=0, edir_y=0, edir_z=0;
  float lambda_nm = 0;
  std::string creator;
  float src_ke_step_MeV = 0;
  int   src_step = -1;        //  parent step number at emission
  int   hit_pmt_id = -1;
};

extern std::unordered_map<const G4Track*, AllPhotonBirthInfo> g_allPhotonBirth;

#endif
