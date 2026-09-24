#ifndef WCSIMPIONSTEPSTREE_HH
#define WCSIMPIONSTEPSTREE_HH

#include "TTree.h"
#include "TFile.h"
#include <string>

// ================================================================
//  pion_steps tree
//
//  One row per hadronic interaction or decay step of any pion track
//  (primary π±, secondary π±, π⁰).  Records the vertex position,
//  KE before and after, process name, and number of secondaries
//  produced.  Elastic scatters appear here as "hadElastic" with
//  ps_n_secondaries == 0 and a KE change of a few percent.
//
//  Source category (ps_src_cat) uses the SourceCategory enum from
//  WCSimAncestryMap.hh:
//    0 = PrimaryPion   1 = SecPionPlus   2 = SecPionMinus
//    3 = PionZero      (others not expected for pion tracks)
//
//  Process names you will see:
//    "hadElastic"        elastic scatter off nucleus
//    "pi+Inelastic"      inelastic: absorption, charge-exchange, scatter
//    "pi-Inelastic"      same for π-
//    "pi0Inelastic"      same for π⁰
//    "Decay"             pion decay (in flight or at rest)
//    "pi-CaptureAtRest"  π- nuclear capture at rest
// ================================================================

extern TTree* pion_steps_tree;

extern int         ps_evt;
extern int         ps_trk;          // Geant4 track ID of the pion
extern int         ps_pdg;          // PDG code of the pion
extern int         ps_src_cat;      // SourceCategory enum (int)

// Interaction vertex (post-step point)
extern float       ps_x_cm;
extern float       ps_y_cm;
extern float       ps_z_cm;
extern float       ps_t_ns;

// Pion kinetic energy before and after the step
extern float       ps_ke_pre_MeV;
extern float       ps_ke_post_MeV;

// Process that defined this step
extern std::string ps_process;

// Number of new secondary tracks created at this step
// (0 for elastic, >0 for inelastic/decay)
extern int         ps_n_secondaries;
extern float ps_dir_pre_x,  ps_dir_pre_y,  ps_dir_pre_z;
extern float ps_dir_post_x, ps_dir_post_y, ps_dir_post_z;
// Lifecycle
void PionStepsTree_Book(TFile* fout);
void PionStepsTree_Fill();
void PionStepsTree_Write();

#endif  // WCSIMPIONSTEPSTREE_HH
