#ifndef WCSIMSECONDARYTRACKTREE_HH
#define WCSIMSECONDARYTRACKTREE_HH

#include "TTree.h"
#include "TFile.h"
#include <string>
#include <unordered_map>

// ================================================================
//  secondary_tracks tree
//
//  One entry per "interesting" secondary particle:
//    • π± (primary or secondary)
//    • π⁰
//    • μ± from any decay
//    • e± created by "Decay"  (Michel electrons)
//    • e± created by "conv"   (from γ → e+e-, i.e. π⁰ chain)
//    • γ  created by "Decay"  (π⁰ → γγ prompt photons)
//    • nuclear fragments      (p, n, heavy ions from breakup)
//
//  Birth kinematics (x0, ke0, …) filled in PreUserTrackingAction.
//  End  kinematics (x1, ke1, …) filled in PostUserTrackingAction.
//
//  src_cat follows the SourceCategory enum in WCSimAncestryMap.hh,
//  but is assigned to *this particle itself*, not to photons it made.
// ================================================================

extern TTree* secondary_tracks_tree;

// ---- Branch variables ----
extern int         st_evt;
extern int         st_trk;
extern int         st_pdg;
extern int         st_parent_trk;
extern int         st_parent_pdg;
extern std::string st_creator;
extern int         st_src_cat;

// Birth
extern float st_x0_cm, st_y0_cm, st_z0_cm, st_t0_ns, st_ke0_MeV;

// End (filled at PostUserTrackingAction)
extern float st_x1_cm, st_y1_cm, st_z1_cm, st_t1_ns, st_ke1_MeV;
extern std::string st_end_process;

// Direct parent info (from ClassifyTrack)
extern int         st_src_trk;
extern int         st_src_pdg;
extern float       st_src_ke0_MeV;
extern std::string st_src_creator;

// Interesting ancestor info (from ClassifyTrack)
extern int         st_anc_trk;
extern int         st_anc_pdg;
extern float       st_anc_ke0_MeV;
extern std::string st_anc_creator;
extern float st_px0, st_py0, st_pz0;
extern float st_px1, st_py1, st_pz1;
// ================================================================
//  SecTrackBuffer – pending-fill helper
//
//  We need two Geant4 hooks (Pre + Post TrackingAction) to collect
//  both the birth and end kinematics before calling Fill().
//  This struct holds birth info until the end of the track.
// ================================================================
struct SecTrackBuffer {
  int   evt        = -1;
  int   trk        = -1;
  int   pdg        = 0;
  int   parent_trk = -1;
  int   parent_pdg = 0;
  int         src_trk     = -1;
  int         src_pdg     = 0;
  float       src_ke0_MeV = 0.f;
  std::string src_creator;

  int         anc_trk     = -1;
  int         anc_pdg     = 0;
  float       anc_ke0_MeV = 0.f;
  std::string anc_creator;
  std::string creator;
  int   src_cat    = 9;    // kOther

  float x0_cm=0, y0_cm=0, z0_cm=0, t0_ns=0, ke0_MeV=0;
  // end fields — left at zero until PostUserTrackingAction fills them
  float x1_cm=0, y1_cm=0, z1_cm=0, t1_ns=0, ke1_MeV=0;
  std::string end_process;
  float px0, py0, pz0;   // birth momentum direction (unit vector)
  float px1, py1, pz1;   // end momentum direction (unit vector)
};

// Map: Geant4 trackID → pending buffer
extern std::unordered_map<int, SecTrackBuffer> g_secTrackBuffer;


// ----------------------------------------------------------------
//  Interface
// ----------------------------------------------------------------

// Decide whether a track should appear in this tree.
// pdg           : G4 PDG code of the track
// creatorProcess: name of the G4 process that created it
bool SecTrack_IsInteresting(int pdg, const std::string& creatorProcess);

// Book tree+branches (call from WCSimRunAction::BeginOfRunAction)
void SecondaryTracksTree_Book(TFile* fout);

// Store birth info; call from PreUserTrackingAction
void SecondaryTracksTree_CacheBirth(int trackID, const SecTrackBuffer& buf);

// Store end info and fill one row; call from PostUserTrackingAction
void SecondaryTracksTree_FillEnd(int         trackID,
                                  float       x1_cm,
                                  float       y1_cm,
                                  float       z1_cm,
                                  float       t1_ns,
                                  float       ke1_MeV,
                                  float       px1,
                                  float       py1,
                                  float       pz1,
                                  const std::string& end_process);

// Write tree to file (call from WCSimRunAction::EndOfRunAction)
void SecondaryTracksTree_Write();

// Clear pending buffer (call from WCSimEventAction::BeginOfEventAction)
void SecondaryTracksTree_ClearBuffer();

#endif  // WCSIMSECONDARYTRACKTREE_HH
