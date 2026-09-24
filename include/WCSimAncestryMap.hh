#ifndef WCSIMANCESTRYMAP_HH
#define WCSIMANCESTRYMAP_HH

#include <unordered_map>
#include <string>

// ================================================================
//  SourceCategory enum
//
//  Identifies the physics origin of a Cherenkov photon or secondary
//  track by walking the full Geant4 track parentage chain.
//  The integer value is stored in pion_photons and secondary_tracks
//  so you can select on it directly in ROOT/Python.
//
//  Used for BOTH:
//    AncestryMap_ClassifyPhoton  – classifies a Cherenkov photon by
//                                  walking up from its parent
//    AncestryMap_ClassifyTrack   – classifies a non-photon track by
//                                  examining itself + walking up
// ================================================================
enum SourceCategory {
  kPrimaryPion   = 0,  // primary π± (G4 parentID == 0)
  kSecPionPlus   = 1,  // secondary π+  (e.g. from hadronic inelastic)
  kSecPionMinus  = 2,  // secondary π-
  kPionZero      = 3,  // π⁰, or particle in the π⁰→γ→e± chain
  kMuon          = 4,  // μ± (from π decay or any source)
  kMichelElec    = 5,  // e± from μ decay ("Decay", parent is μ±)
  kDeltaRay      = 6,  // e± from ionisation (hIoni / eIoni / muIoni / …)
  kNuclear       = 7,  // nuclear-breakup fragment: p, n, α, heavy recoil
  kGamma         = 8,  // γ (from nCapture, bremsstrahlung, unmatched π⁰, …)
  kOther         = 9
};

// Human-readable name (useful for debug prints)
const char* SourceCategoryName(SourceCategory cat);

// ================================================================
//  AncestryInfo – stored once per track when it is first seen
// ================================================================
struct AncestryInfo {
  int   pdg            = 0;
  int   parentID       = -1;
  float ke0_MeV        = 0.f;   // kinetic energy at track creation
  float x0_cm          = 0.f;
  float y0_cm          = 0.f;
  float z0_cm          = 0.f;
  float t0_ns          = 0.f;
  std::string creatorProcess = "";
};

// ================================================================
//  Global ancestry map  (trackID → AncestryInfo)
//
//  Populated : WCSimTrackingAction::PreUserTrackingAction
//  Cleared   : WCSimEventAction::BeginOfEventAction
// ================================================================
extern std::unordered_map<int, AncestryInfo> g_trackAncestry;

// ================================================================
//  Global terminating-process map  (trackID → end process name)
//
//  For a track that comes to rest and dies, the post-step point's
//  GetProcessDefinedStep() reports the along-step process (e.g.
//  "Scintillation"), NOT the AtRest process that actually killed the
//  track.  WCSimSteppingAction resolves the true terminating process
//  (Decay / pi+CaptureAtRest / pi-CaptureAtRest / muMinusCaptureAtRest
//  / pi+Inelastic / Transportation …) for pion/muon tracks and stores
//  it here; WCSimTrackingAction reads it to fill secondary_tracks
//  end_process correctly.  Cleared at BeginOfEventAction.
//  Defined in WCSimAncestryMap.cc.
// ================================================================
extern std::unordered_map<int, std::string> g_trackEndProc;

// Global per-event index shared by all custom pion-analysis trees.
// Assigned in WCSimEventAction::BeginOfEventAction; read in
// WCSimTrackingAction / WCSimSteppingAction.  Defined in WCSimAncestryMap.cc.
extern int g_wcsim_evt;

// Register a new track (call once per track from PreUserTrackingAction)
void AncestryMap_Register(int         trackID,
                           int         pdg,
                           int         parentID,
                           float       ke0_MeV,
                           float       x0_cm,
                           float       y0_cm,
                           float       z0_cm,
                           float       t0_ns,
                           const std::string& creatorProcess);

// Wipe all entries – call at BeginOfEventAction
void AncestryMap_Clear();

// ================================================================
//  AncestryMap_ClassifyPhoton
//
//  Classifies a Cherenkov photon by walking up the ancestry chain
//  from its direct parent (the Cherenkov-emitting particle).
//
//  Output parameters
//  -----------------
//  outSrcTrk    : photon's direct parent track ID
//  outSrcPDG    : direct parent PDG code
//  outSrcKE0    : direct parent birth KE [MeV]
//  outSrcCreator: process that created the direct parent
//
//  outAncTrk    : track ID of the "interesting" ancestor
//                 (e.g. the π⁰ rather than the e+ it produced)
//  outAncPDG    : PDG of the interesting ancestor
//  outAncKE0    : birth KE of interesting ancestor [MeV]
//  outAncCreator: process that created the interesting ancestor
//
//  Returns the SourceCategory.
// ================================================================
SourceCategory AncestryMap_ClassifyPhoton(
    int          parentID,
    int&         outSrcTrk,
    int&         outSrcPDG,
    float&       outSrcKE0,
    std::string& outSrcCreator,
    int&         outAncTrk,
    int&         outAncPDG,
    float&       outAncKE0,
    std::string& outAncCreator);

// ================================================================
//  AncestryMap_ClassifyTrack
//
//  Classifies a non-photon secondary track (π, μ, e±, γ, nuclear
//  fragment) by examining the track itself and, where needed,
//  walking up the ancestry chain.
//
//  Call this from WCSimTrackingAction::PreUserTrackingAction when
//  filling SecTrackBuffer::src_cat.
//
//  Parameters
//  ----------
//  trackID      : Geant4 track ID of the track being classified
//  pdg          : PDG code of the track
//  parentID     : Geant4 parent track ID
//  creatorProcess: name of the G4 process that created the track
//
//  Output parameters
//  -----------------
//  outAncTrk    : track ID of the "interesting" ancestor
//                 (for a Michel e± this is the μ; for a delta-ray
//                  this is the δ-ray itself; for a π this is itself)
//  outAncPDG    : PDG of the interesting ancestor
//  outAncKE0    : birth KE of interesting ancestor [MeV]
//  outAncCreator: process that created the interesting ancestor
//
//  Returns the SourceCategory.
// ================================================================
SourceCategory AncestryMap_ClassifyTrack(
    int          trackID,
    int          pdg,
    int          parentID,
    const std::string& creatorProcess,
    int&         outSrcTrk,        // immediate parent track ID
    int&         outSrcPDG,        // immediate parent PDG
    float&       outSrcKE0,        // immediate parent birth KE [MeV]
    std::string& outSrcCreator,    // process that created the immediate parent
    int&         outAncTrk,        // interesting ancestor track ID
    int&         outAncPDG,        // interesting ancestor PDG
    float&       outAncKE0,        // interesting ancestor birth KE [MeV]
    std::string& outAncCreator);   // process that created interesting ancestor
#endif  // WCSIMANCESTRYMAP_HH
