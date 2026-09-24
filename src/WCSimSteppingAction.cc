#include "WCSimSteppingAction.hh"
#include "G4SystemOfUnits.hh"
#include <stdlib.h>
#include <stdio.h>
#include <WCSimRootEvent.hh>
//#include <G4SIunits.hh>
#include <G4OpticalPhoton.hh>
//#include "G4OpticalPhoton.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4VParticleChange.hh"
#include "G4SteppingVerbose.hh"
#include "G4SteppingManager.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "TTree.h"
#include "TFile.h"
#include "TDirectory.h"
#include "G4RunManager.hh"
#include "G4PhysicalConstants.hh"
//JR EDIT
#include <unordered_map>
#include "G4SystemOfUnits.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"

// --- NEW: pion analysis ---
#include "WCSimAncestryMap.hh"
#include "WCSimPionSourceTree.hh"
#include "WCSimPionStepsTree.hh"
#include <set>

std::unordered_map<const G4Track*, AllPhotonBirthInfo> g_allPhotonBirth;  // definition
int WCSimSteppingAction::primaryMuonTrackID = -1;
int WCSimSteppingAction::primaryMuonEventID = -1;

// g_wcsim_evt is declared in WCSimAncestryMap.hh, defined in WCSimAncestryMap.cc



// ================================================================
// PION ANALYSIS: unified photon birth cache
// One entry per Cherenkov photon, created during the PARENT's step
// so we can capture the exact emitter KE at emission.
// Keyed on const G4Track* (photon track pointer).
// ================================================================



//JR EDIT END






G4int WCSimSteppingAction::n_photons_through_mPMTLV = 0;
G4int WCSimSteppingAction::n_photons_through_acrylic = 0;
G4int WCSimSteppingAction::n_photons_through_gel = 0;
G4int WCSimSteppingAction::n_photons_on_blacksheet = 0;
G4int WCSimSteppingAction::n_photons_on_smallPMT = 0;

///////////////////////////////////////////////
///// BEGINNING OF WCSIM STEPPING ACTION //////
///////////////////////////////////////////////


WCSimSteppingAction::WCSimSteppingAction(WCSimRunAction *myRun, WCSimDetectorConstruction *myDet) : runAction(myRun), det(myDet) {

}


void WCSimSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  
  	//Begin JR Edit

// --- Debug: print when the mu- actually dies (decay vs capture vs exit) ---
/*
{
  auto* trk = aStep->GetTrack();
  if (trk && trk->GetDefinition()->GetPDGEncoding() == 13) { // mu-
    // Detect "about to die": after this step, track won't be alive
    const auto status_after = trk->GetTrackStatus(); // status *after* processes applied to this step
    if (status_after != fAlive) {

      auto* post = aStep->GetPostStepPoint();
      const G4VProcess* proc = post ? post->GetProcessDefinedStep() : nullptr;

      const G4String endProc = proc ? proc->GetProcessName() : "NONE";

      const auto* pv = trk->GetVolume();
      const G4String volName = pv ? pv->GetName() : "NONE";

      const auto* mat = (pv && pv->GetLogicalVolume()) ? pv->GetLogicalVolume()->GetMaterial() : nullptr;
      const G4String matName = mat ? mat->GetName() : "NONE";

      const auto pos = trk->GetPosition();

      G4cout << "[MU-END] evt="
             << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()
             << " trk=" << trk->GetTrackID()
             << " endProc=" << endProc
             << " vol=" << volName
             << " mat=" << matName
             << " x(mm)=" << (pos.x()/mm)
             << " y(mm)=" << (pos.y()/mm)
             << " z(mm)=" << (pos.z()/mm)
             << " t(ns)=" << (trk->GetGlobalTime()/ns)
             << G4endl;
    }
  }
}	
	
	
*/	
	


//Debug code for mu- decay

	
auto* trk_ph = aStep->GetTrack();

// WCSimSteppingAction.cc

if (trk_ph->GetDefinition() == G4OpticalPhoton::Definition()) {

  //const double tcut = 500.0*ns;   // choose based on your DAQ window
  const double tcut = 20000.0*ns;
  const int    ncut = 5000;        // tighten until runtime is acceptable

  if (trk_ph->GetGlobalTime() > tcut || trk_ph->GetCurrentStepNumber() > ncut) {
    trk_ph->SetTrackStatus(fStopAndKill);
  }
}


if (trk_ph->GetDefinition() == G4OpticalPhoton::Definition()) {

  auto* post = aStep->GetPostStepPoint();
  auto* proc = post ? post->GetProcessDefinedStep() : nullptr;

  if (trk_ph->GetKineticEnergy() < 0) {
    static int n = 0;
    if (n++ < 20) {
      G4cout << "NEG OPTICAL TRACK AFTER STEP: "
             << trk_ph->GetKineticEnergy()/MeV << " MeV"
             << " proc=" << (proc ? proc->GetProcessName() : "NONE")
             << " vol=" << (trk_ph->GetVolume() ? trk_ph->GetVolume()->GetName() : "NULL")
             << G4endl;
    }
    trk_ph->SetTrackStatus(fStopAndKill);
    trk_ph->SetKineticEnergy(0.0);
  }
}




  const G4Event* evt = G4RunManager::GetRunManager()->GetCurrentEvent();
const int evid = evt ? evt->GetEventID() : -1;

// reset once per event
if (evid != primaryMuonEventID) {
  primaryMuonEventID = evid;
  primaryMuonTrackID = -1;
}


// Clear per-event to avoid stale track pointers from previous event
static int g_allPhoton_evt_latch = -1;
if (evid != g_allPhoton_evt_latch) {
  g_allPhoton_evt_latch = evid;
  g_allPhotonBirth.clear();
}
//Decide if you want to look for only the primary
bool do_primary_block = true;
// Get current track
const G4Track* stepTrack = aStep->GetTrack();
// Alias reused by the primary/pion blocks below (same pointer as stepTrack)
const G4Track* trk = stepTrack;

if (!stepTrack) do_primary_block = false;




if (do_primary_block){


// ------------------------------------------------------------
// Latch the actual physical primary track ID for this event
// ------------------------------------------------------------
if (primaryMuonTrackID < 0) {

  // reuse the event pointer 'evt' declared at the top of this function
  if (evt) {
    const G4PrimaryVertex* pvtx = evt->GetPrimaryVertex(0);
    if (pvtx) {
      const G4PrimaryParticle* pprim = pvtx->GetPrimary();
      if (pprim) {
        const int primaryPDG = pprim->GetPDGcode();

        // Latch the Geant4 track corresponding to the generated primary
        // Require:
        //   - no parent
        //   - same PDG as event primary
        //   - first step of that track (extra safety)
        if (stepTrack->GetParentID() == 0 &&
            stepTrack->GetDefinition()->GetPDGEncoding() == primaryPDG &&
            stepTrack->GetCurrentStepNumber() == 1) {
          primaryMuonTrackID = stepTrack->GetTrackID();
        }
      }
    }
  }
}

// If we still have not found the event primary, do nothing
if (primaryMuonTrackID < 0) do_primary_block = false;

// ------------------------------------------------------------
// Only consider steps of the latched physical primary
// ------------------------------------------------------------
//const G4Track* trk = aStep->GetTrack();
if (!trk) do_primary_block = false;

if (trk->GetTrackID() != primaryMuonTrackID) do_primary_block = false;
if (trk->GetParentID() != 0) do_primary_block = false;   // extra safety

// Optional extra check: confirm this matches the event primary PDG
//const G4Event* evt = G4RunManager::GetRunManager()->GetCurrentEvent();
if (!evt) do_primary_block = false;

const G4PrimaryVertex* pvtx = evt->GetPrimaryVertex(0);
if (!pvtx) do_primary_block = false;

const G4PrimaryParticle* pprim = pvtx->GetPrimary();
if (!pprim) do_primary_block = false;

const int primaryPDG = pprim->GetPDGcode();
const int pdg = trk->GetDefinition()->GetPDGEncoding();
if (pdg != primaryPDG) do_primary_block = false;
}
//End of guards for primary block

const auto* secs = aStep->GetSecondaryInCurrentStep();




if (secs) {
  for (auto* sec : *secs) {
    if (!sec) continue;
    if (sec->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) continue;

    const G4VProcess* cp = sec->GetCreatorProcess();
    if (!cp) continue;

    // Cache both Cerenkov and Scintillation; remove the second condition
    // if you want to restrict to prompt Cerenkov light only.
    const std::string& procName = cp->GetProcessName();
    if (procName != "Cerenkov" && procName != "Scintillation") continue;

    AllPhotonBirthInfo info;
    info.evt    = evid;
    info.parent = stepTrack->GetTrackID();

    const G4ThreeVector pos = sec->GetPosition();
    info.ex_cm = (float)(pos.x()/cm);
    info.ey_cm = (float)(pos.y()/cm);
    info.ez_cm = (float)(pos.z()/cm);
    info.et_ns = (float)(sec->GetGlobalTime()/ns);

    const G4ThreeVector dir = sec->GetMomentumDirection();
    info.edir_x = (float)dir.x();
    info.edir_y = (float)dir.y();
    info.edir_z = (float)dir.z();

    const G4double Eph = sec->GetTotalEnergy();
    info.lambda_nm = (float)((Eph > 0) ? (h_Planck * c_light / Eph) / nm : 0.0);
    info.creator   = procName;

    // Exact emitter KE at this step — only possible while the parent is stepping
    const G4StepPoint* pre = aStep->GetPreStepPoint();
    info.src_ke_step_MeV = (float)(pre ? pre->GetKineticEnergy()/MeV
                                       : stepTrack->GetKineticEnergy()/MeV);
    info.src_step = stepTrack->GetCurrentStepNumber();
    g_allPhotonBirth[sec] = info;
  }
}

// (Removed dead "Fill AllPhotonsTree" scaffolding that wrote to the
//  never-declared g_primaryPhotonBirth map; superseded by pion_photons.)

// In WCSimSteppingAction.cc, inside UserSteppingAction
const int piPDG = stepTrack->GetDefinition()->GetPDGEncoding();

if (std::abs(piPDG) == 211 || piPDG == 111 || std::abs(piPDG) == 13) {

  const bool isPion = (std::abs(piPDG) == 211 || piPDG == 111);
  const bool isMuon = (std::abs(piPDG) == 13);

  const G4VProcess* postProc =
      aStep->GetPostStepPoint()->GetProcessDefinedStep();

  // A track that dies AT REST has fStopAndKill AND ~zero KE going INTO the
  // fatal step (kePre ~ 0): it has already ranged out, and the AtRest
  // process fires in a separate zero-length step.  An in-flight interaction
  // that also ends the track (pi+Inelastic, decay-in-flight) has
  // fStopAndKill too, and its kePost is also ~0 — so kePre, NOT kePost, is
  // what distinguishes at-rest death from an in-flight interaction.
  const double kePre    = aStep->GetPreStepPoint()->GetKineticEnergy()  / MeV;
  const double kePost   = aStep->GetPostStepPoint()->GetKineticEnergy() / MeV;
  const bool   stopping = (stepTrack->GetTrackStatus() == fStopAndKill);
  const bool   atRest   = stopping && (kePre < 1.0e-3);   // MeV; died from rest

  // ----------------------------------------------------------------
  // Resolve the process for this step.
  //   In flight : GetProcessDefinedStep() is correct (hadElastic,
  //               pi+Inelastic, Decay, …) — use it as-is.
  //   At rest   : GetProcessDefinedStep() reports the along-step process
  //               (e.g. "Scintillation"), NOT the AtRest process, so read
  //               the creator process of the (non-optical) secondaries the
  //               step produced — that carries the real AtRest process name
  //               (Decay / pi+CaptureAtRest / pi-CaptureAtRest /
  //                muMinusCaptureAtRest).  This is the A1/A2 fix.
  //   The secondary-creator lookup is applied ONLY at rest: in-flight steps
  //   also spawn Cherenkov photons (creator "Cerenkov"), so applying it
  //   there would corrupt the process name of a real in-flight interaction.
  // ----------------------------------------------------------------
  std::string procName = postProc ? postProc->GetProcessName() : "";
  if (atRest) {
    std::string atRestProc;
    const auto* stopSecs = aStep->GetSecondaryInCurrentStep();
    if (stopSecs) {
      for (auto* stopSec : *stopSecs) {
        if (!stopSec) continue;
        if (stopSec->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
          continue;   // ignore Cerenkov / scintillation photons
        const G4VProcess* scp = stopSec->GetCreatorProcess();
        if (scp) { atRestProc = scp->GetProcessName(); break; }
      }
    }
    // If nothing resolvable (e.g. a capture with no tracked secondaries),
    // flag it explicitly rather than keeping the misleading along-step name.
    procName = atRestProc.empty() ? "StopAndKill_AtRest" : atRestProc;
  }

  // Record this track's terminating process for secondary_tracks (A2).
  // In-flight world-exit -> "Transportation" (correct); in-flight
  // interaction -> the discrete process; at-rest death -> resolved AtRest
  // process (or the "StopAndKill_AtRest" sentinel).
  if (stopping && !procName.empty())
    g_trackEndProc[stepTrack->GetTrackID()] = procName;

  // Discrete in-flight interactions/decays worth a pion_steps row.
  // NOTE: these are Geant4 *runtime* process-name strings and are
  // physics-list-dependent. If your list names muon single-scattering
  // "CoulombScattering" (or folds it into "msc"), adjust the muon branch.
  const bool isInteresting =
      (isPion &&
       (procName == "hadElastic"       ||
        procName == "pi+Inelastic"     ||
        procName == "pi-Inelastic"     ||
        procName == "pi0Inelastic"     ||
        procName == "Decay"            ||
        procName == "pi-CaptureAtRest" ||
        procName == "pi+CaptureAtRest"))
      ||
      (isMuon &&
       (procName == "Decay"                ||  // Michel decay vertex (mu+/-)
        procName == "muMinusCaptureAtRest" ||  // mu- capture at rest
        procName == "CoulombScat"));           // single hard Coulomb deflection

  // A1 FIX: always record the at-rest death vertex, even when the process
  // name could not be resolved to a known terminator.  In-flight world-exits
  // keep their KE (atRest is false) so they are not recorded.
  const bool record = isInteresting || atRest;

  if (record) {
    ps_evt        = evid;
    ps_trk        = stepTrack->GetTrackID();
    ps_pdg        = piPDG;

    const G4ThreeVector pos = aStep->GetPostStepPoint()->GetPosition();
    ps_x_cm  = pos.x() / cm;
    ps_y_cm  = pos.y() / cm;
    ps_z_cm  = pos.z() / cm;
    ps_t_ns  = aStep->GetPostStepPoint()->GetGlobalTime() / ns;

    ps_ke_pre_MeV  = kePre;
    ps_ke_post_MeV = kePost;
    // Flag the rare ambiguous case (stopped, no resolvable process) so it
    // is visible in the output rather than silently blank.
    ps_process     = procName.empty() ? "StopAndKill_AtRest" : procName;

    const G4ThreeVector dir_pre  = aStep->GetPreStepPoint()->GetMomentumDirection();
    const G4ThreeVector dir_post = aStep->GetPostStepPoint()->GetMomentumDirection();

    ps_dir_pre_x  = (float)dir_pre.x();
    ps_dir_pre_y  = (float)dir_pre.y();
    ps_dir_pre_z  = (float)dir_pre.z();
    ps_dir_post_x = (float)dir_post.x();
    ps_dir_post_y = (float)dir_post.y();
    ps_dir_post_z = (float)dir_post.z();
    // Count physics daughters only.  GetSecondaryInCurrentStep() also holds
    // the Cherenkov/scintillation photons produced along the step, which made
    // this branch scale with the emitter's velocity instead of counting the
    // daughters of the interaction (e.g. 73 at a CoulombScat vertex, which
    // has no daughters at all).
    {
      int nPhysSec = 0;
      const auto* allSecs = aStep->GetSecondaryInCurrentStep();
      if (allSecs) {
        for (auto* aSec : *allSecs) {
          if (!aSec) continue;
          if (aSec->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
            continue;
          ++nPhysSec;
        }
      }
      ps_n_secondaries = nPhysSec;
    }
    {
      auto it = g_trackAncestry.find(stepTrack->GetTrackID());
      if (it != g_trackAncestry.end() && it->second.parentID == 0) {
        // Primary particle: categorise by PDG so a primary muon isn't
        // mislabeled kPrimaryPion.
        if (std::abs(piPDG) == 211 || piPDG == 111) ps_src_cat = (int)kPrimaryPion;
        else if (std::abs(piPDG) == 13)             ps_src_cat = (int)kMuon;
        else                                        ps_src_cat = (int)kOther;
      }
      else {
        int dummySrcTrk=-1, dummySrcPDG=0, dummyAncTrk=-1, dummyAncPDG=0;
        float dummySrcKE=0, dummyAncKE=0;
        std::string dummySrcC, dummyAncC;
        ps_src_cat = (int)AncestryMap_ClassifyTrack(
            stepTrack->GetTrackID(),
            piPDG,
            stepTrack->GetParentID(),
            stepTrack->GetCreatorProcess()
                ? stepTrack->GetCreatorProcess()->GetProcessName() : "",
            dummySrcTrk, dummySrcPDG, dummySrcKE, dummySrcC,
            dummyAncTrk, dummyAncPDG, dummyAncKE, dummyAncC);
      }
    }
    PionStepsTree_Fill();
  }
}
/*
if (do_primary_block){

// ------------------------------------------------------------
// Save Cherenkov optical photons created in this primary's step
// ------------------------------------------------------------
//const auto* secs = aStep->GetSecondaryInCurrentStep();
if (secs) {
  for (auto* sec : *secs) {
    if (sec->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) continue;

    auto* cp = sec->GetCreatorProcess();
    if (!(cp && cp->GetProcessName() == "Cerenkov")) continue;

    a_evt = evt->GetEventID();
    a_trk = sec->GetTrackID();
    a_par = trk->GetTrackID();
    a_muAncestorID = trk->GetTrackID();   // name unchanged, but now means "primary ancestor ID"
    a_pdg = sec->GetDefinition()->GetPDGEncoding(); // optical photon often 0 in WCSim

    auto pos = sec->GetPosition();
    a_ex_cm = pos.x()/cm;
    a_ey_cm = pos.y()/cm;
    a_ez_cm = pos.z()/cm;
    a_et_ns = sec->GetGlobalTime()/ns;

    auto dir = sec->GetMomentumDirection();
    a_edir_x = dir.x();
    a_edir_y = dir.y();
    a_edir_z = dir.z();

    a_Estart_MeV = trk->GetVertexKineticEnergy()/MeV;

    // KE of the primary at this step (pre-step KE)
    const auto* pre = aStep->GetPreStepPoint();
    a_muKE_MeV = pre ? (pre->GetKineticEnergy()/MeV)
                     : (trk->GetKineticEnergy()/MeV);

    // Step number of the primary
    a_muStep = trk->GetCurrentStepNumber();

    const G4double Eph = sec->GetTotalEnergy();
    a_lambda_nm = (h_Planck * c_light / Eph) / nm;

    a_creator = cp->GetProcessName();

    AllPhotonsTree_Fill();
  }
}
}
*/

/* OLD
// latch the primary muon trackID (first time we see it)
const G4Track* stepTrack = aStep->GetTrack();
if (stepTrack && primaryMuonTrackID < 0) {
  const int pdg = stepTrack->GetDefinition()->GetPDGEncoding();
  if (stepTrack->GetParentID() == 0 && (pdg == 13 || pdg == -13)) {
    primaryMuonTrackID = stepTrack->GetTrackID();
  }
}
//auto* trk = aStep->GetTrack();
//For AllPhotonsTree:
//
// ... your existing primary-muon selection ...
if (!trk) return;

// Only consider steps of the primary muon
if (trk->GetTrackID() != primaryMuonTrackID) return;
const int pdg = trk->GetDefinition()->GetPDGEncoding();
if (!(pdg == 13 || pdg == -13)) return;
if (trk->GetParentID() != 0) return;   // extra safety


const auto* secs = aStep->GetSecondaryInCurrentStep();
if (secs) {
  for (auto* s : *secs) {
    if (s->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) continue;
    auto* cp = s->GetCreatorProcess();
    if (!(cp && cp->GetProcessName() == "Cerenkov")) continue;

    a_evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    a_trk = s->GetTrackID();
    a_par = trk->GetTrackID();
    a_muAncestorID = trk->GetTrackID();
    a_pdg = s->GetDefinition()->GetPDGEncoding(); // optical photon is usually 0

    auto pos = s->GetPosition();
    a_ex_cm = pos.x()/cm; a_ey_cm = pos.y()/cm; a_ez_cm = pos.z()/cm;
    a_et_ns = s->GetGlobalTime()/ns;

    auto dir = s->GetMomentumDirection();
    a_edir_x = dir.x(); a_edir_y = dir.y(); a_edir_z = dir.z();

    a_Estart_MeV = trk->GetVertexKineticEnergy()/MeV;
    
    // NEW: KE of the muon at this step (choose pre-step as "energy when step begins")
    const auto* pre = aStep->GetPreStepPoint();
    a_muKE_MeV = pre ? (pre->GetKineticEnergy()/MeV) : (trk->GetKineticEnergy()/MeV);

    // NEW: muon step number (1,2,3,...)
    a_muStep = trk->GetCurrentStepNumber();

    const G4double Eph = s->GetTotalEnergy();
    a_lambda_nm = (h_Planck*c_light/Eph)/nm;

    a_creator = cp->GetProcessName();

    AllPhotonsTree_Fill();
  }
}
*/

//Now fill  ck tree
  static TTree* ck_tree = nullptr;
/*
  static int   b_event   = -1;
  static int   b_trackID = -1;
  static int   b_step    = -1;
  static int   b_pdg     = 0;
  static int   b_nck     = 0;

  static float b_x = 0.f, b_y = 0.f, b_z = 0.f;  // cm
  static float b_t = 0.f;                        // ns
  static float b_ke = 0.f;                       // MeV
  static float b_step_len = 0.f;  // cm
*/

// --- Per-step quantities ---
static int   b_event    = -1;
static int   b_trackID  = -1;
static int   b_step     = -1;
static int   b_pdg      = 0;
static int   b_n_ck     = 0;

static float b_x = 0.f, b_y = 0.f, b_z = 0.f, b_t = 0.f;
static float b_ke = 0.f;
static float b_ke_pre  = 0.f;
static float b_ke_post = 0.f;
static float b_edep    = 0.f;
static float b_step_len = 0.f;

// ROOT requires std::string to persist
static std::string b_mat;
static std::string b_step_proc;

   // Create the tree the first time we enter this function
  if (!ck_tree) {
    ck_tree = new TTree("ck_step", "Cherenkov photons per step for primary mu-");

    ck_tree->Branch("event",    &b_event,   "event/I");
  ck_tree->Branch("trackID",  &b_trackID, "trackID/I");
  ck_tree->Branch("step",     &b_step,    "step/I");
  ck_tree->Branch("pdg",      &b_pdg,     "pdg/I");
  ck_tree->Branch("n_ck",     &b_n_ck,    "n_ck/I");

  ck_tree->Branch("x", &b_x, "x/F");
  ck_tree->Branch("y", &b_y, "y/F");
  ck_tree->Branch("z", &b_z, "z/F");
  ck_tree->Branch("t", &b_t, "t/F");

  ck_tree->Branch("ke",       &b_ke,       "ke/F");
  ck_tree->Branch("ke_pre",   &b_ke_pre,   "ke_pre/F");
  ck_tree->Branch("ke_post",  &b_ke_post,  "ke_post/F");
  ck_tree->Branch("edep",     &b_edep,     "edep/F");
  ck_tree->Branch("step_len", &b_step_len, "step_len/F");

  ck_tree->Branch("mat",        &b_mat);
  ck_tree->Branch("step_proc", &b_step_proc);
    /*
    ck_tree->Branch("event",   &b_event,   "event/I");
    ck_tree->Branch("trackID", &b_trackID, "trackID/I");
    ck_tree->Branch("step",    &b_step,    "step/I");
    ck_tree->Branch("pdg",     &b_pdg,     "pdg/I");
    ck_tree->Branch("n_ck",    &b_nck,     "n_ck/I");

    ck_tree->Branch("x",       &b_x,       "x/F");
    ck_tree->Branch("y",       &b_y,       "y/F");
    ck_tree->Branch("z",       &b_z,       "z/F");
    ck_tree->Branch("t",       &b_t,       "t/F");
    ck_tree->Branch("ke",      &b_ke,      "ke/F");
    ck_tree->Branch("step_len", &b_step_len, "step_len/F");
*/
  }

  //auto* trk = aStep->GetTrack();
  //
 if (trk && trk->GetTrackID() == 1 && std::abs(trk->GetDefinition()->GetPDGEncoding()) == 13) {
//if (!trk) return;

// Primary muon only
//if (trk->GetTrackID() != 1) return;
//if (trk->GetDefinition()->GetPDGEncoding() != 13) return;

// Step info
auto* pre  = aStep->GetPreStepPoint();
auto* post = aStep->GetPostStepPoint();

// Count Cherenkov photons created in this step
b_n_ck = 0;
//const auto* secs = aStep->GetSecondaryInCurrentStep();
if (secs) {
  for (auto* sec : *secs) {
    if (sec->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) continue;
    auto* cp = sec->GetCreatorProcess();
    if (cp && cp->GetProcessName() == "Cerenkov") b_n_ck++;
  }
}

// Fill scalars
b_event   = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
b_trackID = trk->GetTrackID();
b_step    = trk->GetCurrentStepNumber();
b_pdg     = trk->GetDefinition()->GetPDGEncoding();

b_x = pre->GetPosition().x() / cm;
b_y = pre->GetPosition().y() / cm;
b_z = pre->GetPosition().z() / cm;
b_t = pre->GetGlobalTime()  / ns;

b_ke      = trk->GetKineticEnergy() / MeV;
b_ke_pre  = pre->GetKineticEnergy()  / MeV;
b_ke_post = post->GetKineticEnergy() / MeV;

b_edep     = aStep->GetTotalEnergyDeposit() / MeV;
b_step_len = aStep->GetStepLength() / cm;

// Material + process
auto* mat = pre->GetMaterial();
b_mat = mat ? mat->GetName() : "NONE";

auto* proc = post->GetProcessDefinedStep();
b_step_proc = proc ? proc->GetProcessName() : "NONE";

// Fill tree (DO NOT return early — keep all steps)
ck_tree->Fill();

}
  
  /*
  // --- Select the primary mu- (your run: trackID==1, pdg==13) ---
  auto* trk = aStep->GetTrack();
  if (!trk) return;

  if (trk->GetTrackID() != 1) return;
  if (trk->GetDefinition()->GetPDGEncoding() != 13) return;

  // --- Count Cherenkov optical photons produced in THIS step ---
  int n_ck = 0;
  const auto* secs = aStep->GetSecondaryInCurrentStep();
  if (secs) {
    for (auto* s : *secs) {
      if (s->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) continue;
      auto* cp = s->GetCreatorProcess();
      if (cp && cp->GetProcessName() == "Cerenkov") n_ck++;
    }
  }

  //if (n_ck <= 0) return;  // only write steps that produced Cherenkov photons

  // --- Fill branches ---
  auto* pre = aStep->GetPreStepPoint();
  if (!pre) return;

  b_event   = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  b_trackID = trk->GetTrackID();
  b_step    = trk->GetCurrentStepNumber();
  b_pdg     = trk->GetDefinition()->GetPDGEncoding();
  b_nck     = n_ck;

  const auto pos = pre->GetPosition();
  b_x  = pos.x() / cm;
  b_y  = pos.y() / cm;
  b_z  = pos.z() / cm;
  b_t  = pre->GetGlobalTime() / ns;
  b_ke = trk->GetKineticEnergy() / MeV;
  b_step_len = aStep->GetStepLength() / cm;

  ck_tree->Fill();
*/


/*    
    auto* trk = aStep->GetTrack();
  if (!trk) return;
  
  if (trk->GetCurrentStepNumber() == 1 && trk->GetParentID() == 0) {
  int pdg = trk->GetDefinition()->GetPDGEncoding();
  double ke = trk->GetKineticEnergy()/MeV;
  double vz = trk->GetVertexPosition().z()/cm;

  auto* vol = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  G4String vname = vol ? vol->GetName() : "NONE";

  G4cout << "PRIMARY TRACK STEP1: "
         << "trackID=" << trk->GetTrackID()
         << " pdg=" << pdg
         << " KE=" << ke << " MeV"
         << " vertex_z=" << vz << " cm"
         << " vol=" << vname
         << G4endl;
}


  if (trk->GetTrackID() != 1) return;
if (trk->GetDefinition()->GetPDGEncoding() != 13) return; // mu-

int n_ck = 0;
const auto* secs = aStep->GetSecondaryInCurrentStep();
if (secs) {
  for (auto* s : *secs) {
    if (s->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) continue;
    auto* cp = s->GetCreatorProcess();
    if (cp && cp->GetProcessName() == "Cerenkov") n_ck++;
  }
}

if (n_ck > 0) {
  auto* pre = aStep->GetPreStepPoint();
  G4cout << "MU CK STEP: step=" << trk->GetCurrentStepNumber()
         << " z=" << pre->GetPosition().z()/cm << " cm"
         << " KE=" << trk->GetKineticEnergy()/MeV << " MeV"
         << " n_ck=" << n_ck
         << G4endl;
}


  // Reset per event
  int evtID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  if (evtID != g_lastEventID) {
    g_lastEventID = evtID;
    g_realMuonTrackID = -1;
  }

  // --- Step 1: identify the real GPS mu- (once per event) ---
  if (g_realMuonTrackID < 0) {
    if (trk->GetCurrentStepNumber() == 1 &&
        trk->GetDefinition()->GetPDGEncoding() == 13 &&  // mu-
        trk->GetParentID() == 0)                         // primary
    {
      double ke = trk->GetKineticEnergy()/MeV;
      double vz = trk->GetVertexPosition().z()/cm;

      // Reject the dummy: KE=0 and absurd vertex (your dummy is at -100000 cm)
      // Keep "real" muons: KE ~ O(100 MeV) and vertex in detector-scale region
      if (ke > 1.0 && std::abs(vz) < 10000.0) { // 10 km is a generous sanity bound
        g_realMuonTrackID = trk->GetTrackID();
        G4cout << "REAL GPS MUON FOUND: event=" << evtID
               << " trackID=" << g_realMuonTrackID
               << " vertex_z=" << vz << " cm"
               << " KE=" << ke << " MeV"
               << G4endl;
      } else {
        // Debug print so you can see the dummy candidates too
        G4cout << "REJECT MUON CANDIDATE: event=" << evtID
               << " trackID=" << trk->GetTrackID()
               << " vertex_z=" << vz << " cm"
               << " KE=" << ke << " MeV"
               << G4endl;
      }
    }
  }

  // If we still haven't identified it, nothing else to do yet
  if (g_realMuonTrackID < 0) return;

  // --- Step 2: only process steps of the real muon ---
  if (trk->GetTrackID() != g_realMuonTrackID) return;
  if (trk->GetDefinition()->GetPDGEncoding() != 13) return;

  // Count Cherenkov optical photons produced this step
  int n_ck = 0;
  const auto* secs = aStep->GetSecondaryInCurrentStep();
  if (secs) {
    for (auto* s : *secs) {
      if (s->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) continue;
      auto* cp = s->GetCreatorProcess();
      if (cp && cp->GetProcessName() == "Cerenkov") n_ck++;
    }
  }

  if (n_ck > 0) {
    auto* pre = aStep->GetPreStepPoint();
    auto pos = pre->GetPosition();
    G4cout << "Muon step " << trk->GetCurrentStepNumber()
           << " z=" << pos.z()/cm << " cm"
           << " KE=" << trk->GetKineticEnergy()/MeV << " MeV"
           << " n_ck=" << n_ck
           << G4endl;
  }
*/
    //End JR edit
	
    const G4Event *event = G4EventManager::GetEventManager()->GetConstCurrentEvent();
    if(event->IsAborted() || event->GetEventID() < 0)
      return;
  //DISTORTION must be used ONLY if INNERTUBE or INNERTUBEBIG has been defined in BidoneDetectorConstruction.cc
  
  const G4Track* track       = aStep->GetTrack();
 

  //JR EDIT BEGIN


  // Not used:
  //const G4Event* evt = G4RunManager::GetRunManager()->GetCurrentEvent();
  //G4VPhysicalVolume* volume  = track->GetVolume();
  //G4SDManager* SDman   = G4SDManager::GetSDMpointer();
  //G4HCofThisEvent* HCE = evt->GetHCofThisEvent();

  
  // Debug for photon tracking
  G4StepPoint* thePrePoint = aStep->GetPreStepPoint();
  G4VPhysicalVolume* thePrePV = thePrePoint->GetPhysicalVolume();

  G4StepPoint* thePostPoint = aStep->GetPostStepPoint();
  G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();

  //G4OpBoundaryProcessStatus boundaryStatus=Undefined;
  //static G4ThreadLocal G4OpBoundaryProcess* boundary=NULL;  //doesn't work and needs #include tls.hh from Geant4.9.6 and beyond
  G4OpBoundaryProcess* boundary=NULL;
  
  //find the boundary process only once
  if(!boundary){
    G4ProcessManager* pm
      = aStep->GetTrack()->GetDefinition()->GetProcessManager();
    G4int nprocesses = pm->GetProcessListLength();
    G4ProcessVector* pv = pm->GetProcessList();
    G4int i;
    for( i=0;i<nprocesses;i++){
      if((*pv)[i]->GetProcessName()=="OpBoundary"){
	boundary = (G4OpBoundaryProcess*)(*pv)[i];
	break;
      }
    }
  }

  G4ParticleDefinition *particleType = track->GetDefinition();
  if(particleType == G4OpticalPhoton::OpticalPhotonDefinition() && thePrePV && thePostPV){
    
    //JR EDIT BEGIN
  


// ================================================================
// PION ANALYSIS: fill pion_photons tree when any photon terminates
// ================================================================
if (track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {

  const bool dying = (track->GetTrackStatus() == fStopAndKill);
  auto itp = g_allPhotonBirth.find(track);

  if (dying && itp != g_allPhotonBirth.end()) {
    const AllPhotonBirthInfo& info = itp->second;

    // ---- Look up PE acceptance — only fill the row if this photon made a PE ----
    int   pe_pmt;
    float pe_time;
    float pe_lambda;
    const bool madePE = PionPhotonsTree_GetPE(info.evt, track->GetTrackID(),
                                                pe_pmt, pe_time, pe_lambda);

    if (madePE) {
      // ---- Ancestry classification ----
      int         srcTrk=-1, srcPDG=0, ancTrk=-1, ancPDG=0;
      float       srcKE0=0,  ancKE0=0;
      std::string srcCreator, ancCreator;
      SourceCategory cat = AncestryMap_ClassifyPhoton(
          info.parent,
          srcTrk, srcPDG, srcKE0, srcCreator,
          ancTrk, ancPDG, ancKE0, ancCreator);

      // ---- Fill pion_photons branches ----
      pp_evt = info.evt;
      pp_trk = track->GetTrackID();

      pp_ex_cm = info.ex_cm;  pp_ey_cm = info.ey_cm;  pp_ez_cm = info.ez_cm;
      pp_et_ns = info.et_ns;
      pp_edir_x = info.edir_x; pp_edir_y = info.edir_y; pp_edir_z = info.edir_z;
      pp_lambda_nm = info.lambda_nm;
      pp_creator   = info.creator;

      pp_src_trk     = srcTrk;
      pp_src_pdg     = srcPDG;
      pp_src_ke_MeV  = info.src_ke_step_MeV;
      pp_src_ke0_MeV = srcKE0;
      pp_src_step    = info.src_step;          // NEW
      pp_src_creator = srcCreator;

      pp_anc_trk     = ancTrk;
      pp_anc_pdg     = ancPDG;
      pp_anc_ke0_MeV = ancKE0;
      pp_anc_creator = ancCreator;

      pp_src_cat = (int)cat;

      pp_hit_pmt      = info.hit_pmt_id;
      pp_end_boundary = thePostPV ? std::string(thePostPV->GetName()) : "NONE";
      const G4VProcess* endP = aStep->GetPostStepPoint() ?
                               aStep->GetPostStepPoint()->GetProcessDefinedStep() : nullptr;
      pp_end_process  = endP ? std::string(endP->GetProcessName()) : "NONE";

      // NEW: PE-bridge information
      pp_made_pe      = 1;
      pp_pe_pmt       = pe_pmt;
      pp_pe_time_ns   = pe_time;
      pp_pe_lambda_nm = pe_lambda;

      PionPhotonsTree_Fill();
    }

    g_allPhotonBirth.erase(itp);
  }
  else if (dying) {
    g_allPhotonBirth.erase(track);
  }
}

/*
    // ============================================================
// JR: Michel photon finalization (PMT hit or photon termination)
// ============================================================
{
  // Only care about photons we cached at creation time
  const int trkID = track->GetTrackID();
  if (g_michelPhotonFilled.find(track) != g_michelPhotonFilled.end()) return;

  auto it = g_michelPhotonBirth.find(track);
  if (it != g_michelPhotonBirth.end()) {

  // -------- Detect a PMT "hit" (best-effort from geometry) --------
int hitPMT = -1;

// Require valid PVs
if (thePrePV && thePostPV) {
  const G4String preName  = thePrePV->GetName();
  const G4String postName = thePostPV->GetName();

  // Only treat as a PMT hit when the photon ENTERS a "pmt" volume
  // (avoids pmt->pmt spam and internal steps)
  const bool enteredPMT = (postName == "pmt") && (preName != "pmt");

  if (enteredPMT) {
    int pmt_in_module = -1;
    int module_copy   = -1;

    const auto* postPoint = aStep->GetPostStepPoint();
    const auto* touch = postPoint ? postPoint->GetTouchable() : nullptr;

    if (touch) {
      // Walk the touchable history to find both the inner PMT copy and the module copy
      for (int d = 0; d <= touch->GetHistoryDepth(); ++d) {
        const auto* vol = touch->GetVolume(d);
        if (!vol) continue;

        const G4String nm = vol->GetName();

        if (nm == "pmt") {
          pmt_in_module = touch->GetCopyNumber(d);
        } else if (nm == "WCMultiPMT") {
          module_copy = touch->GetCopyNumber(d);
        }
      }

      // TEMP unique ID: module*100 + pmtIndex (good enough for now)
      if (pmt_in_module >= 0) {
        hitPMT = (module_copy >= 0) ? (module_copy * 100 + pmt_in_module)
                                    : pmt_in_module;
      }
	
      // Small debug dump for first few PMT entries
      static int dumped = 0;
      if (dumped < 5) {
        G4cout << "[MICHEL-PMT] evt="
               << (G4RunManager::GetRunManager()->GetCurrentEvent()
                     ? G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()
                     : -1)
               << " photonTrk=" << track->GetTrackID()
               << " pre=" << preName
               << " post=" << postName
               << " moduleCopy=" << module_copy
               << " pmtCopy=" << pmt_in_module
               << " hitPMT=" << hitPMT
               << G4endl;

        G4cout << "[TOUCH] depth=" << touch->GetHistoryDepth()
               << " trk=" << track->GetTrackID() << G4endl;
        for (int d = 0; d <= touch->GetHistoryDepth(); ++d) {
          const auto* v = touch->GetVolume(d);
          if (!v) continue;
          G4cout << "  d=" << d
                 << " name=" << v->GetName()
                 << " copy=" << touch->GetCopyNumber(d)
                 << G4endl;
        }
        dumped++;
      }
      
    }
  }
}

// At this point hitPMT is either -1 (no PMT entry) or our derived ID.
// Don't return here — just leave hitPMT=-1 if it wasn't a PMT entry.
// Your finalization code below should store hitPMT and end_boundary/process accordingly.


    // -------- Detect a PMT "hit" (best-effort from geometry) --------
    // In WCSim, PV names often contain "WCPMT" / "PMT" depending on geometry
    int hitPMT = -1;
    if (thePostPV) {
      const std::string postName = thePostPV->GetName();
      const int copyNo = thePostPV->GetCopyNo();

     // accept either "pmt" or WCMultiPMT (you’re seeing both)
bool isPMT = (postName == "pmt") || (postName == "WCMultiPMT");
if (!isPMT) return;

int hitPMT = copyNo;   // for your geometry this is what you’ve been printing as "copy="

const auto* post = aStep->GetPostStepPoint();
const auto* touch = post->GetTouchable();
static int dumped = 0;
if (dumped < 5) {
  G4cout << "[TOUCH] depth=" << touch->GetHistoryDepth()
         << " trk=" << track->GetTrackID() << G4endl;
  for (int d = 0; d <= touch->GetHistoryDepth(); ++d) {
    G4cout << "  d=" << d
           << " name=" << touch->GetVolume(d)->GetName()
           << " copy=" << touch->GetCopyNumber(d)
           << G4endl;
  }
  dumped++;
} 
   }

    

//G4cout << "[DBG-FIND]  key(ptr)=" << track << " trkID=" << track->GetTrackID() << G4endl;
    // -------- Determine if photon is ending this step --------
    // We fill if it:
    //   (a) enters a PMT volume (hitPMT>=0), OR
    //   (b) is being killed (Absorption, Boundary kill, etc.)
    const bool killedNow = (track->GetTrackStatus() == fStopAndKill);

    if (hitPMT >= 0 || killedNow) {

      const auto& info = it->second;


       if (hitPMT >= 0 && info.evt < 10) {
  G4cout << "[MICHEL-PMT] evt=" << info.evt
         << " photonTrk=" << trkID
         << " postPV=" << thePostPV->GetName()
         << " copy=" << hitPMT
         << G4endl;
}

      // Only proceed if it's a michel photon we cached

      // ---- Fill your michel_photons branches (mcp_*) ----
      mcp_evt        = info.evt;
      //mcp_trk        = info.trk;
      mcp_trk = track->GetTrackID();   // real, nonzero
      mcp_parent     = info.parent;
      mcp_pdg        = info.pdg;
      mcp_michel_trk = info.michel_trk;

      mcp_sx_cm = info.sx_cm;
      mcp_sy_cm = info.sy_cm;
      mcp_sz_cm = info.sz_cm;
      mcp_st_ns = info.st_ns;

      mcp_dir_x = info.dir_x;
      mcp_dir_y = info.dir_y;
      mcp_dir_z = info.dir_z;

      mcp_lambda_nm = info.lambda_nm;

      mcp_hit_pmt = hitPMT;

      // Where it ended / what boundary it was at
      mcp_end_boundary = (thePostPV ? thePostPV->GetName() : "NONE");

      // Process that limited this step (Absorption, OpBoundary, etc.)
      const G4VProcess* endP = thePostPoint ? thePostPoint->GetProcessDefinedStep() : nullptr;
      mcp_end_process = endP ? endP->GetProcessName() : "NONE";

      // Creator stored from birth cache (should be "Cerenkov")
      mcp_creator = info.creator;

      MichelPhotonsTree_Fill();

      // Remove so we do not double-fill this photon
      g_michelPhotonFilled.insert(track);
      g_michelPhotonBirth.erase(it);
    }
  }
}
*/
//JR EDIT END

    if( (thePrePV->GetName().find("MultiPMT") != std::string::npos) &&
	(thePostPV->GetName().find("vessel") != std::string::npos) )
      n_photons_through_mPMTLV++;
    
    if( (thePostPV->GetName().find("container") != std::string::npos) &&
	(thePrePV->GetName().find("vessel") != std::string::npos))
      n_photons_through_acrylic++;
    
    if( (thePrePV->GetName().find("container") != std::string::npos) )
      n_photons_through_gel++;

    if( (thePrePV->GetName().find("container") != std::string::npos) &&
	(thePostPV->GetName().find("inner") != std::string::npos) )
      n_photons_on_blacksheet++;

    if( (thePrePV->GetName().find("container") != std::string::npos) &&
	(thePostPV->GetName().find("pmt") != std::string::npos) )
      n_photons_on_smallPMT++;
	

    /*
    if( (thePrePV->GetName().find("pmt") != std::string::npos)){
      G4cout << "Photon between " << thePrePV->GetName() <<
	" and " << thePostPV->GetName() << " because " << 
	thePostPoint->GetProcessDefinedStep()->GetProcessName() << 
	" and boundary status: " <<  boundary->GetStatus() << " with track status " << track->GetTrackStatus() << G4endl;
      
	}*/



    if(track->GetTrackStatus() == fStopAndKill){
      if(boundary->GetStatus() == NoRINDEX){
	G4cout << "Optical photon is killed because of missing refractive index in either " << thePrePoint->GetMaterial()->GetName() << " or " << thePostPoint->GetMaterial()->GetName() <<
	  " (transition from " << thePrePV->GetName() << " to " << thePostPV->GetName() << ")" <<
	  " : could also be caused by Overlaps with volumes with logicalBoundaries." << G4endl;
	
      }
      /* Debug :  
      if( (thePrePV->GetName().find("PMT") != std::string::npos) ||
	  (thePrePV->GetName().find("pmt") != std::string::npos)){
	
	if(boundary->GetStatus() != StepTooSmall){
	  //	if(thePostPoint->GetProcessDefinedStep()->GetProcessName() != "Transportation")
	  G4cout << "Killed photon between " << thePrePV->GetName() <<
	    " and " << thePostPV->GetName() << " because " << 
	    thePostPoint->GetProcessDefinedStep()->GetProcessName() << 
	    " and boundary status: " <<  boundary->GetStatus() <<
	    G4endl;
	}
	}	*/
      
    }
  }



}


G4int WCSimSteppingAction::G4ThreeVectorToWireTime(G4ThreeVector *pos3d,
						    G4ThreeVector lArPos,
						    G4ThreeVector start,
						    G4int i)
{
  G4double x0 = start[0]-lArPos[0];
  G4double y0 = start[1]-lArPos[1];
//   G4double y0 = 2121.3;//mm
  G4double z0 = start[2]-lArPos[2];

  G4double dt=0.8;//mm
//   G4double midt = 2651.625;
  G4double pitch = 3;//mm
//   G4double midwir = 1207.10;
  G4double c45 = 0.707106781;
  G4double s45 = 0.707106781;

  G4double w1;
  G4double w2;
  G4double t;

//   G4double xField(0.);
//   G4double yField(0.);
//   G4double zField(0.);

//   if(detector->getElectricFieldDistortion())
//     {
//       Distortion(pos3d->getX(),
// 		 pos3d->getY());
      
//       xField = ret[0];
//       yField = ret[1];
//       zField = pos3d->getZ();
      
//       w1 = (int)(((zField+z0)*c45 + (x0-xField)*s45)/pitch); 
//       w2 = (int)(((zField+z0)*c45 + (x0+xField)*s45)/pitch); 
//       t = (int)(yField+1);

//       //G4cout<<" x orig "<<pos3d->getX()<<" y orig "<<(pos3d->getY()+y0)/dt<<G4endl;
//       //G4cout<<" x new "<<xField<<" y new "<<yField<<G4endl;
//     }
//   else 
//     {
      
  w1 = (int) (((pos3d->getZ()+z0)*c45 + (x0-pos3d->getX())*s45)/pitch); 
  w2 = (int)(((pos3d->getZ()+z0)*c45 + (x0+pos3d->getX())*s45)/pitch); 
  t  = (int)((pos3d->getY()+y0)/dt +1);
//     }

  if (i==0)
    return (int)w1;
  else if (i==1)
    return (int)w2;
  else if (i==2)
    return (int)t;
  else return 0;
} 


void WCSimSteppingAction::Distortion(G4double /*x*/,G4double /*y*/)
{
 
//   G4double theta,steps,yy,y0,EvGx,EvGy,EField,velocity,tSample,dt;
//   y0=2121.3;//mm
//   steps=0;//1 mm steps
//   tSample=0.4; //micros
//   dt=0.8;//mm
//   LiquidArgonMedium medium;  
//   yy=y;
//   while(y<y0 && y>-y0 )
//     {
//       EvGx=FieldLines(x,y,1);
//       EvGy=FieldLines(x,y,2);
//       theta=atan(EvGx/EvGy);
//       if(EvGy>0)
// 	{
// 	  x+=sin(theta);
// 	  y+=cos(theta);
// 	}
//       else
// 	{
// 	  y-=cos(theta);
// 	  x-=sin(theta);
// 	}
//       EField=sqrt(EvGx*EvGx+EvGy*EvGy);//kV/mm
//       velocity=medium.DriftVelocity(EField*10);// mm/microsec
//       steps+=1/(tSample*velocity);

//       //G4cout<<" step "<<steps<<" x "<<x<<" y "<<y<<" theta "<<theta<<" Gx "<<eventaction->Gx->Eval(x,y)<<" Gy "<<eventaction->Gy->Eval(x,y)<<" EField "<<EField<<" velocity "<<velocity<<G4endl;
//     }

//   //numbers
//   //EvGx=FieldLines(0,1000,1);
//   //EvGy=FieldLines(0,1000,2);
//   //EField=sqrt(EvGx*EvGx+EvGy*EvGy);//kV/mm
//   //velocity=medium.DriftVelocity(EField*10);// mm/microsec
//   //G4double quenching;
//   //quenching=medium.QuenchingFactor(2.1,EField*10);
//   //G4cout<<" Gx "<<EvGx<<" Gy "<<EvGy<<" EField "<<EField<<" velocity "<<velocity<<" quenching "<<quenching<<G4endl;


  //ret[0]=5;
//   if(yy>0)
//     ret[1]=2*y0/dt -steps;
//   else
//     ret[1]=steps; 
}


double WCSimSteppingAction::FieldLines(G4double /*x*/,G4double /*y*/,G4int /*coord*/)
{ //0.1 kV/mm = field
  //G4double Radius=302;//mm
//   G4double Radius=602;//mm
//   if(coord==1) //x coordinate
//     return  (0.1*(2*Radius*Radius*x*abs(y)/((x*x+y*y)*(x*x+y*y))));
//   else //y coordinate
//     return 0.1*((abs(y)/y)*(1-Radius*Radius/((x*x+y*y)*(x*x+y*y))) + abs(y)*(2*Radius*Radius*y/((x*x+y*y)*(x*x+y*y))));
  return 0;
}
