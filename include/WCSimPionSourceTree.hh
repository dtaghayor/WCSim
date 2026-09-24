#ifndef WCSIMPIONTREE_HH
#define WCSIMPIONTREE_HH

#include "TTree.h"
#include "TFile.h"
#include <string>

// ================================================================
//  pion_photons tree
//
//  One row per Cherenkov photon (from any particle in the event)
//  that terminates — either at a PMT face or by absorption.
//
//  The src_cat integer lets you select on origin without string
//  matching; see SourceCategory in WCSimAncestryMap.hh.
//
//  KE columns
//  ----------
//  pp_src_ke_MeV  = KE of the Cherenkov emitter AT THE EMISSION STEP.
//  pp_src_ke0_MeV = KE of the Cherenkov emitter AT ITS BIRTH.
//  pp_anc_ke0_MeV = KE of the "interesting ancestor" at its birth.
//
//  PMT / PE fate
//  -------------
//  pp_hit_pmt      = WCSim tube ID where the photon arrived (else -1).
//  pp_made_pe      = 1 if photon survived QE/collection efficiency.
//  pp_pe_pmt       = tube ID where the PE was recorded (else -1).
//  pp_pe_time_ns   = recorded PE arrival time [ns].
//  pp_pe_lambda_nm = wavelength at PE acceptance [nm].
//
//  pp_src_step     = step number of the parent at the moment of
//                    photon emission.
// ================================================================

extern TTree* pion_photons_tree;

// Event / track identity
extern int   pp_evt;
extern int   pp_trk;

// Photon birth (emission vertex)
extern float pp_ex_cm, pp_ey_cm, pp_ez_cm;
extern float pp_et_ns;
extern float pp_edir_x, pp_edir_y, pp_edir_z;
extern float pp_lambda_nm;
extern std::string pp_creator;

// Direct Cherenkov emitter (photon's immediate G4 parent)
extern int         pp_src_trk;
extern int         pp_src_pdg;
extern float       pp_src_ke_MeV;
extern float       pp_src_ke0_MeV;
extern int         pp_src_step;       // NEW
extern std::string pp_src_creator;

// Interesting ancestor (root of the physics chain)
extern int         pp_anc_trk;
extern int         pp_anc_pdg;
extern float       pp_anc_ke0_MeV;
extern std::string pp_anc_creator;

// Source category (SourceCategory enum cast to int)
extern int pp_src_cat;

// Photon fate
extern int         pp_hit_pmt;
extern std::string pp_end_boundary;
extern std::string pp_end_process;

// NEW: PE-bridge information
extern int         pp_made_pe;
extern int         pp_pe_pmt;
extern float       pp_pe_time_ns;
extern float       pp_pe_lambda_nm;

// Lifecycle helpers
void PionPhotonsTree_Book(TFile* fout);
void PionPhotonsTree_Fill();
void PionPhotonsTree_Write();

// ================================================================
//  Per-event PE-acceptance bridge
//
//  WCSimWCSD::ProcessHits calls _NotePE after the QE roll succeeds.
//  WCSimSteppingAction calls _GetPE when the photon track terminates.
//  The map is reset per event.
// ================================================================
void PionPhotonsTree_ResetEvent(int evt);
void PionPhotonsTree_NotePE(int evt, int photon_trk,
                             int pe_pmt, float pe_time_ns, float pe_lambda_nm);
bool PionPhotonsTree_GetPE(int evt, int photon_trk,
                            int& pe_pmt, float& pe_time_ns, float& pe_lambda_nm);

#endif  // WCSIMPIONTREE_HH
