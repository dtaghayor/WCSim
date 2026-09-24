#include "WCSimTrackingAction.hh"
#include "WCSimTrajectory.hh"
#include "G4ParticleTypes.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "G4RunManager.hh"
#include "WCSimTrackInformation.hh"
#include "WCSimTrackingMessenger.hh"
#include "WCSimPrimaryGeneratorAction.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"

// --- NEW: pion analysis ---
#include <cstdlib>   // std::abs(int)
#include "WCSimAncestryMap.hh"
#include "WCSimSecondaryTrackTree.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

WCSimTrackingAction::WCSimTrackingAction()
{

  ProcessList.insert("Decay") ;                         // Michel e- from pi+ and mu+
  ProcessList.insert("conv") ;                         // Products of gamma conversion

  //ProcessList.insert("muMinusCaptureAtRest") ;          // Includes Muon decay from K-shell: for Michel e- from mu0. This dominates/replaces the mu- decay (noticed when switching off this process in PhysicsList)                                                   // TF: IMPORTANT: ONLY USE FROM G4.9.6 onwards because buggy/double counting before.
  ////////// ToDo: switch ON the above when NuPRISM uses G4 >= 4.9.6
  ProcessList.insert("nCapture");

//   ProcessList.insert("conv");

  // F. Nova One can check here if the photon comes from WLS
  ProcessList.insert("OpWLS");

  ParticleList.insert(111); // pi0
  ParticleList.insert(211); // pion+
  ParticleList.insert(-211);
  ParticleList.insert(321);
  ParticleList.insert(-321); // kaon-
  ParticleList.insert(311); // kaon0
  ParticleList.insert(-311); // kaon0 bar
  //ParticleList.insert(22); // I add photons (B.Q)
  ParticleList.insert(11); // e-
  ParticleList.insert(-11); // e+
  ParticleList.insert(13); // mu-
  ParticleList.insert(-13); // mu+
  // don't put gammas there or there'll be too many

  //TF: add protons and neutrons
  ParticleList.insert(2212);
  ParticleList.insert(2112);

  percentageOfCherenkovPhotonsToDraw = 0.0;
#ifndef WCSIM_SAVE_PHOTON_HISTORY
  SAVE_PHOTON_HISTORY = false;
#else
  SAVE_PHOTON_HISTORY = true;
#endif

  messenger = new WCSimTrackingMessenger(this);

  // Max time for radioactive decay:
  fMaxTime    = 1. * CLHEP::second;
  fTime_birth = 0.;
}

WCSimTrackingAction::~WCSimTrackingAction(){;}

void WCSimTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{

// ================================================================
// PION ANALYSIS — Step 1: register every track in the ancestry map.
// This must happen before any of this track's secondaries are stepped.
// ================================================================
if (aTrack) {
  const G4ThreeVector vtx = aTrack->GetVertexPosition();
  const G4VProcess*   cp  = aTrack->GetCreatorProcess();

  AncestryMap_Register(
      aTrack->GetTrackID(),
      aTrack->GetDefinition()->GetPDGEncoding(),
      aTrack->GetParentID(),
      (float)(aTrack->GetKineticEnergy() / MeV),
      (float)(vtx.x() / cm),
      (float)(vtx.y() / cm),
      (float)(vtx.z() / cm),
      (float)(aTrack->GetGlobalTime() / ns),
      cp ? std::string(cp->GetProcessName()) : "");
}

// ================================================================
// PION ANALYSIS — Step 2: cache birth info for any "interesting"
// secondary (π, μ, Michel e, conv e, nuclear fragment, …).
// The end kinematics are added in PostUserTrackingAction.
// ================================================================
if (aTrack) {
  const int         pdg     = aTrack->GetDefinition()->GetPDGEncoding();
  const G4VProcess* cp      = aTrack->GetCreatorProcess();
  const std::string creator = cp ? std::string(cp->GetProcessName()) : "";

  // B1 FIX: a track is suspended and resumed many times, and Geant4 calls
  // PreUserTrackingAction on EVERY resume.  Cache birth info only the first
  // time this track is seen, so x0/y0/z0/t0/ke0 stay the true track origin
  // instead of being overwritten with the resume point on each segment.
  if (SecTrack_IsInteresting(pdg, creator) &&
      !g_secTrackBuffer.count(aTrack->GetTrackID())) {

    // Look up parent PDG from the ancestry map
    int parentPDG = 0;
    {
      auto it = g_trackAncestry.find(aTrack->GetParentID());
      if (it != g_trackAncestry.end()) parentPDG = it->second.pdg;
    }

    SourceCategory cat = kOther;
    int         srcTrk=-1, srcPDG=0, ancTrk=-1, ancPDG=0;
    float       srcKE0=0.f, ancKE0=0.f;
    std::string srcC, ancC;

    if (aTrack->GetParentID() == 0) {
      // Primary particle: categorise by PDG so a non-pion primary (e.g. a
      // muon) isn't mislabeled kPrimaryPion.
      if (std::abs(pdg) == 211 || pdg == 111) cat = kPrimaryPion;
      else if (std::abs(pdg) == 13)           cat = kMuon;
      else                                    cat = kOther;
      srcTrk = aTrack->GetTrackID();
      srcPDG = pdg;
      ancTrk = aTrack->GetTrackID();
      ancPDG = pdg;
      srcKE0 = ancKE0 = (float)(aTrack->GetKineticEnergy() / MeV);
      srcC = ancC = "";
    } else {
      cat = AncestryMap_ClassifyTrack(
          aTrack->GetTrackID(),
          pdg,
          aTrack->GetParentID(),
          creator,
          srcTrk, srcPDG, srcKE0, srcC,
          ancTrk, ancPDG, ancKE0, ancC);
    }

    const G4ThreeVector vtx = aTrack->GetVertexPosition();
    SecTrackBuffer buf;
    buf.evt        = g_wcsim_evt;
    buf.trk        = aTrack->GetTrackID();
    buf.pdg        = pdg;
    buf.parent_trk = aTrack->GetParentID();
    buf.parent_pdg = parentPDG;
    buf.creator    = creator;
    buf.src_cat    = (int)cat;
    buf.x0_cm      = (float)(vtx.x() / cm);
    buf.y0_cm      = (float)(vtx.y() / cm);
    buf.z0_cm      = (float)(vtx.z() / cm);
    buf.t0_ns      = (float)(aTrack->GetGlobalTime() / ns);
    buf.ke0_MeV    = (float)(aTrack->GetKineticEnergy() / MeV);
    buf.src_trk     = srcTrk;
    buf.src_pdg     = srcPDG;
    buf.src_ke0_MeV = srcKE0;
    buf.src_creator = srcC;
    buf.anc_trk     = ancTrk;
    buf.anc_pdg     = ancPDG;
    buf.anc_ke0_MeV = ancKE0;
    buf.anc_creator = ancC;
    const G4ThreeVector dir0 = aTrack->GetMomentumDirection();
    buf.px0 = (float)dir0.x();
    buf.py0 = (float)dir0.y();
    buf.pz0 = (float)dir0.z();
    SecondaryTracksTree_CacheBirth(aTrack->GetTrackID(), buf);
  }
}



  if ( aTrack->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()
       || G4UniformRand() < percentageOfCherenkovPhotonsToDraw/100. )
    {
      WCSimTrajectory* thisTrajectory = new WCSimTrajectory(aTrack);
      thisTrajectory->SetSavePhotonTrack(true); // mark to save photon track in ROOT
      fpTrackingManager->SetTrajectory(thisTrajectory);
      fpTrackingManager->SetStoreTrajectory(true);
    }
  else if (SAVE_PHOTON_HISTORY)
    {
      // Keep the trajectory but not saving in ROOT
      WCSimTrajectory* thisTrajectory = new WCSimTrajectory(aTrack);
      thisTrajectory->SetSavePhotonTrack(false);
      fpTrackingManager->SetTrajectory(thisTrajectory);
      fpTrackingManager->SetStoreTrajectory(true);
    }
  else 
    fpTrackingManager->SetStoreTrajectory(false);
	
  WCSimPrimaryGeneratorAction *primaryGenerator = (WCSimPrimaryGeneratorAction *) (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  if(!primaryGenerator->IsConversionFound()) {
    if(aTrack->GetParentID()==0){
      primaryID = aTrack->GetTrackID();
    }
    else if(aTrack->GetParentID() == primaryID) {
      if (aTrack->GetCreatorProcess()->GetProcessName() == "conv") {
          primaryGenerator->FoundConversion();
      }
      //JR EDIT - UNCOMMENT FOR ORIGINAL CODE BELOW
      //G4EventManager::GetEventManager()->AbortCurrentEvent();
      //G4EventManager::GetEventManager()->GetNonconstCurrentEvent()->SetEventAborted();
    }
  }

  // Kill nucleus generated after TrackID 1
  G4ParticleDefinition* particle = aTrack->GetDefinition();
  G4String name   = particle->GetParticleName();
  G4double fCharge = particle->GetPDGCharge();

  G4Track* tr = (G4Track*) aTrack;
  if ( aTrack->GetTrackID() == 1 ) {
  	// Re-initialize time
  	fTime_birth = 0;
  	// Ask G4 to kill the track when all secondary are done (will exclude other decays)
  	if ( fCharge > 2. )
  		tr->SetTrackStatus(fStopButAlive);
  }

  if ( aTrack->GetTrackID() == 2 ) {
  	// First track of the decay save time
  	fTime_birth = aTrack->GetGlobalTime(); 
  }



  // ---------------------------------------------------------
//JR EDIT END
}

void WCSimTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
// ================================================================
// PION ANALYSIS — fill secondary_tracks end info for this track
// ================================================================
// B1 FIX: PostUserTrackingAction also fires on every SUSPEND, not just at
// the true end of the track.  Filling there produced one row per track
// segment (486 contiguous rows for a single 235 MeV muon).  Emit a row only
// once the track has actually finished, so secondary_tracks holds exactly
// one row per track: x0/ke0 = birth, x1/ke1 = death.
const G4TrackStatus secTrkStatus = aTrack ? aTrack->GetTrackStatus() : fAlive;
const bool secTrkFinished = (secTrkStatus == fStopAndKill ||
                             secTrkStatus == fKillTrackAndSecondaries);

if (aTrack && secTrkFinished && g_secTrackBuffer.count(aTrack->GetTrackID())) {

  const G4ThreeVector endpos = aTrack->GetPosition();

  const G4VProcess* endProc = nullptr;
  if (aTrack->GetStep() && aTrack->GetStep()->GetPostStepPoint())
    endProc = aTrack->GetStep()->GetPostStepPoint()->GetProcessDefinedStep();
  std::string endProcName = endProc ? std::string(endProc->GetProcessName()) : "NONE";

  // A2 FIX: GetProcessDefinedStep() on the post-step point reports the
  // along-step process (e.g. "Scintillation") for a track that dies at
  // rest, not the AtRest process that actually terminated it. For pion/
  // muon tracks, WCSimSteppingAction has already resolved the true
  // terminating process into g_trackEndProc — prefer it when present.
  {
    auto ep = g_trackEndProc.find(aTrack->GetTrackID());
    if (ep != g_trackEndProc.end() && !ep->second.empty())
      endProcName = ep->second;
  }
  const G4ThreeVector dir1 = aTrack->GetMomentumDirection();

  SecondaryTracksTree_FillEnd(
      aTrack->GetTrackID(),
      (float)(endpos.x() / cm),
      (float)(endpos.y() / cm),
      (float)(endpos.z() / cm),
      (float)(aTrack->GetGlobalTime() / ns),
      (float)(aTrack->GetKineticEnergy() / MeV),
      (float)dir1.x(), (float)dir1.y(), (float)dir1.z(),  
      endProcName);
}

  // added by M Fechner
  const G4VProcess* creatorProcess = aTrack->GetCreatorProcess();

  WCSimTrackInformation* anInfo;
  if (aTrack->GetUserInformation())
    anInfo = (WCSimTrackInformation*)(aTrack->GetUserInformation());   //eg. propagated to all secondaries blelow.
  else anInfo = new WCSimTrackInformation();

    // Simplified tracking: track all primaries, particles from the chosen list of particle types or creation processes,
    // Optionally, all particles producing cherekov hits and their parents, grandparents, etc. will be added later.

    // is it a primary ?
    // is the process in the set ? 
    // is the particle in the set ?
    if( aTrack->GetParentID()==0
	|| ((creatorProcess!=0) && ProcessList.count(creatorProcess->GetProcessName()))
	|| (ParticleList.count(aTrack->GetDefinition()->GetPDGEncoding()))
      )
    {
      // if so the track is worth saving
      anInfo->WillBeSaved(true);
    }
  else {
    anInfo->WillBeSaved(false);
  }


  G4Track* theTrack = (G4Track*)aTrack;
  theTrack->SetUserInformation(anInfo);

  // pass parent trajectory to children
  G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
  WCSimTrajectory *currentTrajectory = (WCSimTrajectory*)fpTrackingManager->GimmeTrajectory();
  if(currentTrajectory && !anInfo->GetMyTrajectory())
    anInfo->SetMyTrajectory(currentTrajectory);
  if(secondaries)
  {
    size_t nSeco = secondaries->size();
    if(nSeco>0)
    {
      for(size_t i=0;i<nSeco;i++)
      {
        WCSimTrackInformation* infoSec = new WCSimTrackInformation();
        infoSec->SetParentTrajectory(anInfo->GetMyTrajectory());
        (*secondaries)[i]->SetUserInformation(infoSec);
      }
    } 
  }

  // If this track produces a hit, traverse back through parent trajectories to flag that the parents produce a hit
  if (anInfo->GetProducesHit() && saveHitProducingTracks){
      WCSimTrajectory* parentTrajectory = anInfo->GetParentTrajectory();
      while(parentTrajectory != 0 && !parentTrajectory->GetProducesHit()){
          if (parentTrajectory->GetPDGEncoding()==0 && !parentTrajectory->GetSavePhotonTrack()) break; // do not save unwanted photon tracks
          parentTrajectory->SetProducesHit(true);
          parentTrajectory = parentTrajectory->GetParentTrajectory();
      }
  }

  if (currentTrajectory)
    //   if (aTrack->GetDefinition()->GetPDGCharge() == 0) 
  {

    G4ThreeVector currentPosition      = aTrack->GetPosition();
    G4VPhysicalVolume* currentVolume   = aTrack->GetVolume();

    currentTrajectory->SetStoppingPoint(currentPosition);
    currentTrajectory->SetStoppingVolume(currentVolume);

    if (aTrack->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
    {
      currentTrajectory->SetSaveFlag(anInfo->isSaved());// mark it for WCSimEventAction ;
      // don't save the optical photon tracks themselves simply when they produce a hit, since that info is already saved in WCSimRootCherenkovHitTime
      // optical photon tracks can still be saved if wanted by explicitly adding appropriate entries to the ParticleList or ProcessList via mac file commands
      currentTrajectory->SetProducesHit(anInfo->GetProducesHit());
    }
    else if (currentTrajectory->GetSavePhotonTrack()) // only save the wanted photon tracks
      currentTrajectory->SetSaveFlag(anInfo->isSaved()); 
    else 
      currentTrajectory->SetSaveFlag(false);
  }
	
  WCSimPrimaryGeneratorAction *primaryGenerator = (WCSimPrimaryGeneratorAction *) (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  if(!primaryGenerator->IsConversionFound() && 
     aTrack->GetTrackID() == primaryID &&
     aTrack->GetStep()->GetPostStepPoint()->GetProcessDefinedStep() &&
     aTrack->GetStep()->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() == "conv"){
      for(int i=0; i<2; i++){
          primaryGenerator->SetConversionProductParticle(i, secondaries->at(i)->GetParticleDefinition());
          primaryGenerator->SetConversionProductMomentum(i, secondaries->at(i)->GetMomentum());
      }
  }
}




