#ifndef WCSimIonizingStepsTree_h
#define WCSimIonizingStepsTree_h 1

#include <string>
#include "TTree.h"
class TFile;

// ================================================================
//  ionizing_steps  --  ONE ROW PER STEP (not per vertex)
//
//  Every step of the primary and of its charged daughters (the same
//  "interesting" set as secondary_tracks: pi, mu, Michel/conv e,
//  p/n/nuclear).  Intended for multiple-scattering, path-length and
//  dE/dx studies.
//
//  DISABLED BY DEFAULT -- this tree is far larger than the per-vertex
//  trees.  Enable in the macro with:
//      /WCSim/PionAnalysis/RecordIonizingSteps true
// ================================================================

extern TTree* ionizing_steps_tree;

extern int         is_evt;
extern int         is_trk;
extern int         is_pdg;
extern int         is_parent_trk;
extern int         is_src_cat;
extern int         is_step;          // Geant4 step number within the track

extern float       is_x_cm;          // post-step point
extern float       is_y_cm;
extern float       is_z_cm;
extern float       is_t_ns;

extern float       is_ke_pre_MeV;
extern float       is_ke_post_MeV;
extern float       is_edep_MeV;      // total energy deposited in this step

extern float       is_step_len_cm;   // length of THIS step
extern float       is_track_len_cm;  // cumulative track length after this step

extern float       is_beta_pre;      // v/c before the step
extern float       is_beta_post;     // v/c after the step
extern int         is_above_ckv;     // 1 if beta_post > 1/n (still radiating)

extern int         is_n_cherenkov;   // Cherenkov photons created in this step
extern std::string is_process;       // step-limiting process

void IonizingStepsTree_Book(TFile* fout);
void IonizingStepsTree_Fill();
void IonizingStepsTree_Write();

#endif
