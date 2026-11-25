#ifndef WCSimPionAnalysisMessenger_h
#define WCSimPionAnalysisMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithADouble;

// ================================================================
//  Runtime configuration for the pion-analysis trees.
//  Defined in WCSimPionAnalysisMessenger.cc, set from the macro via
//  the commands below.
// ================================================================

// Record every step of the primary and its charged daughters into the
// ionizing_steps tree.  OFF by default: it is a large tree (one row per
// step, not per vertex) and is only wanted for dedicated multiple-scattering
// / path-length studies.
//   /WCSim/PionAnalysis/RecordIonizingSteps true
extern bool g_recordIonizingSteps;

// Refractive index used for the Cherenkov threshold (beta_thr = 1/n).
// Default 1.33 (water at ~400 nm).
//   /WCSim/PionAnalysis/CherenkovIndex 1.33
extern double g_cherenkovIndex;

class WCSimPionAnalysisMessenger : public G4UImessenger
{
public:
  WCSimPionAnalysisMessenger();
  ~WCSimPionAnalysisMessenger();

  void SetNewValue(G4UIcommand* command, G4String newValue);

private:
  G4UIdirectory*      fDir;
  G4UIcmdWithABool*   fRecordIonizingSteps;
  G4UIcmdWithADouble* fCherenkovIndex;
};

#endif
