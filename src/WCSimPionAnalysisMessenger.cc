#include "WCSimPionAnalysisMessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4ios.hh"

// ---- runtime configuration (declared extern in the header) ----
bool   g_recordIonizingSteps = false;   // OFF by default: large tree
double g_cherenkovIndex      = 1.33;    // water at ~400 nm


WCSimPionAnalysisMessenger::WCSimPionAnalysisMessenger()
{
  fDir = new G4UIdirectory("/WCSim/PionAnalysis/");
  fDir->SetGuidance("Controls for the pion-analysis trees.");

  fRecordIonizingSteps =
      new G4UIcmdWithABool("/WCSim/PionAnalysis/RecordIonizingSteps", this);
  fRecordIonizingSteps->SetGuidance(
      "Record every step of the primary and its charged daughters into the");
  fRecordIonizingSteps->SetGuidance(
      "ionizing_steps tree. Default false -- this tree is one row per STEP");
  fRecordIonizingSteps->SetGuidance(
      "and makes output files much larger. Enable only for path-length /");
  fRecordIonizingSteps->SetGuidance("multiple-scattering studies.");
  fRecordIonizingSteps->SetParameterName("record", true);
  fRecordIonizingSteps->SetDefaultValue(false);

  fCherenkovIndex =
      new G4UIcmdWithADouble("/WCSim/PionAnalysis/CherenkovIndex", this);
  fCherenkovIndex->SetGuidance(
      "Refractive index used for the Cherenkov threshold (beta_thr = 1/n)");
  fCherenkovIndex->SetGuidance(
      "when filling the ckv_* branches of secondary_tracks. Default 1.33.");
  fCherenkovIndex->SetParameterName("n", true);
  fCherenkovIndex->SetDefaultValue(1.33);
}


WCSimPionAnalysisMessenger::~WCSimPionAnalysisMessenger()
{
  delete fRecordIonizingSteps;
  delete fCherenkovIndex;
  delete fDir;
}


void WCSimPionAnalysisMessenger::SetNewValue(G4UIcommand* command,
                                             G4String     newValue)
{
  if (command == fRecordIonizingSteps) {
    g_recordIonizingSteps = fRecordIonizingSteps->GetNewBoolValue(newValue);
    G4cout << "[PionAnalysis] RecordIonizingSteps = "
           << (g_recordIonizingSteps ? "true" : "false") << G4endl;
  }
  else if (command == fCherenkovIndex) {
    const double n = fCherenkovIndex->GetNewDoubleValue(newValue);
    if (n > 1.0) {
      g_cherenkovIndex = n;
      G4cout << "[PionAnalysis] CherenkovIndex = " << g_cherenkovIndex
             << "  (beta_thr = " << 1.0/g_cherenkovIndex << ")" << G4endl;
    } else {
      G4cout << "[PionAnalysis] ignoring CherenkovIndex " << n
             << " -- must be > 1" << G4endl;
    }
  }
}
