#include "WCSimWCSD.hh"
#include "WCSimSteppingAction.hh"
#include "WCSimDetectorConstruction.hh"
#include "WCSimTrackInformation.hh"
#include "G4ParticleTypes.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "Randomize.hh"
#include "G4ios.hh"

#include "WCSimPionSourceTree.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <unordered_set>

#include <sstream>

#include "TTree.h"
#include <string>


WCSimWCSD::WCSimWCSD(G4String CollectionName,
                     G4String name,
                     WCSimDetectorConstruction* myDet,
                     G4String inDetectorElement)
     :G4VSensitiveDetector(name), detectorElement(inDetectorElement)
{
  // Place the name of this collection on the list.  We can have more than one
  // in principle.  CollectionName is a vector.

  // Note there is some sort of problem here.  If I use the name
  // Which has a "/" in it, I can find this collection later using 
  // GetCollectionID()

  collectionName.insert(CollectionName);
  
  fdet = myDet;
  
  HCID = -1;
}

WCSimWCSD::~WCSimWCSD() {}

void WCSimWCSD::Initialize(G4HCofThisEvent* HCE)
{
  // Make a new hits collection. With the name we set in the constructor
  hitsCollection = new WCSimWCHitsCollection
    (SensitiveDetectorName,collectionName[0]);

  // This is a trick.  We only want to do this once.  When the program
  // starts HCID will equal -1.  Then it will be set to the pointer to
  // this collection.

  
  // Get the Id of the "0th" collection
  if (HCID<0){
    HCID =  GetCollectionID(0); 
  }  
  // Add it to the Hit collection of this event.
  HCE->AddHitsCollection( HCID, hitsCollection );  

  // Initialize the Hit map to all tubes not hit.
  PMTHitMap.clear();
  // Trick to access the static maxPE variable.  This will go away with the 
  // variable.

  WCSimWCHit* newHit = new WCSimWCHit();
  newHit->SetMaxPe(0);
  delete newHit;
}

G4bool WCSimWCSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{ 
  G4StepPoint*       preStepPoint = aStep->GetPreStepPoint();
  G4TouchableHandle  theTouchable = preStepPoint->GetTouchableHandle();
  G4VPhysicalVolume* thePhysical  = theTouchable->GetVolume();





  // Triggered by photons from InteriorWCPMT hitting photocathode 
  if (thePhysical->GetName()=="InteriorWCPMT") return ProcessHits_boundary(aStep,0);

  //XQ 3/30/11 try to get the local position try to add the position and direction
  G4ThreeVector worldPosition  = preStepPoint->GetPosition();
  G4ThreeVector localPosition  = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPosition);
  G4ThreeVector worldDirection = preStepPoint->GetMomentumDirection();
  G4ThreeVector localDirection = theTouchable->GetHistory()->GetTopTransform().TransformAxis(worldDirection);

  WCSimTrackInformation* trackinfo 
    = (WCSimTrackInformation*)(aStep->GetTrack()->GetUserInformation());

  G4int parentSavedTrackID = -1;
  G4float       photonStartTime;
  G4float       photonStartEnergy;
  G4float       photonEndEnergy;
  G4ThreeVector photonStartPos;
  G4ThreeVector photonStartDir;
  
  parentSavedTrackID   = aStep->GetTrack()->GetParentID();
  photonStartTime      = aStep->GetTrack()->GetGlobalTime() - aStep->GetTrack()->GetLocalTime(); // current time minus elapsed time of track
  photonStartEnergy    = aStep->GetTrack()->GetVertexKineticEnergy();
  photonEndEnergy      = aStep->GetTrack()->GetKineticEnergy();
  photonStartPos       = aStep->GetTrack()->GetVertexPosition();
  photonStartDir       = aStep->GetTrack()->GetVertexMomentumDirection();
 
  // Need to create a NONE string in case the Hit has no creatorProcess, such a Dark Noise Hit.
  const G4VProcess* process = aStep->GetTrack()->GetCreatorProcess();
  ProcessType_t photonCreatorProcess(kUnknownProcess);
  if (process) {
    photonCreatorProcess = WCSimEnumerations::ProcessTypeStringToEnum(process->GetProcessName());
  }

  G4int    trackID           = aStep->GetTrack()->GetTrackID();
  G4String volumeName        = aStep->GetTrack()->GetVolume()->GetName();
  
  G4double energyDeposition  = aStep->GetTotalEnergyDeposit();
  G4double hitTime           = aStep->GetPreStepPoint()->GetGlobalTime();

  G4ParticleDefinition *particleDefinition = 
    aStep->GetTrack()->GetDefinition();

  if ( particleDefinition != G4OpticalPhoton::OpticalPhotonDefinition() 
       && energyDeposition == 0.0) 
    return false;
  // MF : I don't see why other particles should register hits
  // they don't in skdetsim. 
  if ( particleDefinition != G4OpticalPhoton::OpticalPhotonDefinition())
    return false;

  G4String WCCollectionName;
  if(detectorElement=="tank") WCCollectionName = fdet->GetIDCollectionName();
  else if (detectorElement=="tankPMT2") WCCollectionName = fdet->GetIDCollectionName2();
  else if (detectorElement=="OD") WCCollectionName = fdet->GetODCollectionName();

  //= fdet->GetIDCollectionName();
  //WCSimPMTObject *PMT = fdet->GetPMTPointer(WCIDCollectionName);//B.Q
  //G4cout << "PMT associated to collection = " << PMT->GetPMTName() << G4endl;
  //if(volumeName == ) WCIDCollectionName = fdet->GetIDCollectionName();
  //G4String WCIDCollectionName2 = fdet->GetIDCollectionName2();

  // M Fechner : too verbose
  //  if (aStep->GetTrack()->GetTrackStatus() == fAlive)G4cout << "status is fAlive\n";
  if ((aStep->GetTrack()->GetTrackStatus() == fAlive )
      &&(particleDefinition == G4OpticalPhoton::OpticalPhotonDefinition()))
    return false;
  
  // TF: Problem: photons can go through the sensitive detector (glass)
  // and be killed in the next volume by Absorption, eg. in the BlackSheet, but will then
  // still be counted as hits here. So, either recode as extended/optical/LXe cathode implementation
  // or easier: make sure they cross the cathode, as the logical Boundary "cathode" can't be a
  // sensitive detector (I think).
  
  G4StepPoint        *postStepPoint = aStep->GetPostStepPoint();
  G4VPhysicalVolume  *postVol = postStepPoint->GetPhysicalVolume();
  //if (thePhysical)  G4cout << " thePrePV:  " << thePhysical->GetName()  << G4endl;
  //if (postVol) G4cout << " thePostPV: " << postVol->GetName() << G4endl;
  
  //Optical Photon must pass through glass into PMT interior!
  // What about the other way around? TF: current interior won't keep photons alive like in reality
  // Not an issue yet, because then interior needs to be a sensitive detector, when postStepPoint is the glass.
  if(postVol->GetName() != "InteriorWCPMT")
    return false;
  

  //  if ( particleDefinition ==  G4OpticalPhoton::OpticalPhotonDefinition() ) 
  // G4cout << volumeName << " hit by optical Photon! " << G4endl;
    
  // Make the tubeTag string based on the replica numbers
  // See WCSimDetectorConstruction::DescribeAndRegisterPMT() for matching
  // tag construction.

  std::stringstream tubeTag;

  // Start tubeTag with mother to distinguish different PMT hierarchies
//  G4LogicalVolume *theMother = thePhysical->GetMotherLogical();
//  if (theMother != NULL)
//    tubeTag << theMother->GetName() << ":";

//  tubeTag << thePhysical->GetName(); 
  for (G4int i = theTouchable->GetHistoryDepth()-1 ; i >= 0; i--){
    tubeTag << ":" << theTouchable->GetVolume(i)->GetName();
    tubeTag << "-" << theTouchable->GetCopyNumber(i);
  }
  //  tubeTag << ":" << theTouchable->GetVolume(i)->GetCopyNo(); 

  // Debug:
  //G4cout << "================================================" << G4endl;
  //G4cout << tubeTag.str() << G4endl;
  //G4cout << "================================================" << G4endl;

  // Get the tube ID from the tubeTag
  G4int replicaNumber;
  if(detectorElement=="tank") replicaNumber = WCSimDetectorConstruction::GetTubeID(tubeTag.str());
  else if(detectorElement=="tankPMT2") replicaNumber = WCSimDetectorConstruction::GetTubeID2(tubeTag.str());
  else if(detectorElement=="OD") replicaNumber = WCSimDetectorConstruction::GetODTubeID(tubeTag.str());
  else G4cout << "detectorElement not defined..." << G4endl;
  // Mark this photon as having hit a PMT in the pion_photons birth cache
  {
    const G4Track* currentTrack = aStep->GetTrack();
    auto it_birth = g_allPhotonBirth.find(currentTrack);
    if (it_birth != g_allPhotonBirth.end()) {
      it_birth->second.hit_pmt_id = (int)replicaNumber;
    }
  }
// NOTE: PionPhotonsTree_NotePE is called below, INSIDE the accepted-hit block
// (past both the QE roll and the angular-efficiency roll), next to AddPe, so
// the per-PE truth recorded here matches the PEs that feed the digitizer.


  G4double theta_angle = 0.;
  G4double effectiveAngularEfficiency = 0.;

  //XQ Add the wavelength there
  G4double  wavelength = (2.0*M_PI*197.3)/( aStep->GetTrack()->GetTotalEnergy()/eV);
  G4double ratio = 1.;
  G4double maxQE = 0.;
  G4double photonQE = 0.;
  if (fdet->GetPMT_QE_Method()==1 || fdet->GetPMT_QE_Method() == 4){
    photonQE = 1.1;
  }else if (fdet->GetPMT_QE_Method()==2){
    // maxQE = fdet->GetPMTQE(WCIDCollectionName,wavelength,0,200,700,ratio);
    maxQE = fdet->GetPMTQE(WCCollectionName,wavelength,0,200,660,ratio);
    photonQE = fdet->GetPMTQE(volumeName, wavelength,1,200,660,ratio);
    photonQE = photonQE/maxQE;
  }else if (fdet->GetPMT_QE_Method() == 3){
    ratio = 1./(1.-0.25);
    photonQE = fdet->GetPMTQE(WCCollectionName, wavelength,1,200,660,ratio);
    //photonQE = 0.5; //JR try making QE uniform for a test
  }

// G4cout<<fdet->GetPMT_QE_Method()<<G4endl;
 // G4cout<<photonQE<<G4endl;  
  if (G4UniformRand() <= photonQE){
    
     G4double local_x = localPosition.x();
     G4double local_y = localPosition.y();
     G4double local_z = localPosition.z();
     theta_angle = acos(fabs(local_z)/sqrt(pow(local_x,2)+pow(local_y,2)+pow(local_z,2)))/3.1415926*180.;
     effectiveAngularEfficiency = fdet->GetPMTCollectionEfficiency(theta_angle, volumeName);
    // G4cout << "Theta "<< theta_angle
    // << "ang_eff " << effectiveAngularEfficiency <<G4endl;	     
/*
     static int qprint = 0;
if (qprint < 20) {
  G4cout << "JR QE debug: wl=" << wavelength
         << " photonQE=" << photonQE
         << " theta=" << theta_angle
         << " collEff=" << effectiveAngularEfficiency
         << " method=" << fdet->GetPMT_QE_Method()
         << G4endl;
  qprint++;
}
*/

     if (G4UniformRand() <= effectiveAngularEfficiency || fdet->UsePMT_Coll_Eff()==0){
       //Retrieve the pointer to the appropriate hit collection. 
       //Since volumeName is the same as the SD name, this works. 
       G4SDManager* SDman = G4SDManager::GetSDMpointer();
       G4RunManager* Runman = G4RunManager::GetRunManager();
       G4int collectionID = SDman->GetCollectionID(volumeName);
       const G4Event* currentEvent = Runman->GetCurrentEvent();
       G4HCofThisEvent* HCofEvent = currentEvent->GetHCofThisEvent();
       hitsCollection = (WCSimWCHitsCollection*)(HCofEvent->GetHC(collectionID));

       // mark the track as having produced a hit
       if(!trackinfo)
           trackinfo = new WCSimTrackInformation();
       trackinfo->SetProducesHit(true);

       // Record this accepted PE (per-PE truth) for the pion_photons PE-bridge.
       // Keyed on (eventID, photon trackID) -- the same key GetPE uses at photon
       // termination in WCSimSteppingAction.  Placed past both the QE roll and
       // the angular-efficiency roll so it matches the PEs that feed the digitizer.
       PionPhotonsTree_NotePE(currentEvent ? currentEvent->GetEventID() : -1,
                              trackID,
                              (int)replicaNumber,
                              (float)(hitTime / ns),
                              (float)wavelength);
      /*
       // ===== JR: compute cos_inc for this accepted PE (once) =====
float cos_inc_f = -999.f;
{
  auto* pre  = aStep->GetPreStepPoint();
  auto* post = aStep->GetPostStepPoint();
  const G4ThreeVector dir = pre->GetMomentumDirection(); // unit

  // pre-side normal
  auto preTouch = pre->GetTouchableHandle();
  auto preSolid = preTouch->GetSolid();
  auto preLocal = preTouch->GetHistory()->GetTopTransform().TransformPoint(pre->GetPosition());
  auto npre_w   = preTouch->GetHistory()->GetTopTransform()
                    .TransformAxis(preSolid->SurfaceNormal(preLocal)).unit();
  double cos_pre = -dir.dot(npre_w);

  // post-side normal
  auto postTouch = post->GetTouchableHandle();
  auto postSolid = postTouch->GetSolid();
  auto postLocal = postTouch->GetHistory()->GetTopTransform().TransformPoint(post->GetPosition());
  auto npost_w   = postTouch->GetHistory()->GetTopTransform()
                    .TransformAxis(postSolid->SurfaceNormal(postLocal)).unit();
  double cos_post = -dir.dot(npost_w);

  double cos_inc = std::max(cos_pre, cos_post);
  cos_inc = std::max(-1.0, std::min(1.0, cos_inc));
  cos_inc_f = (float)cos_inc;
}
//JR END EDIT
*/

       // If this tube hasn't been hit add it to the collection
       if (PMTHitMap[replicaNumber] == 0)
       //if (PMTHitMap.find(replicaNumber) == PMTHitMap.end())  TF attempt to fix
	 {
	   WCSimWCHit* newHit = new WCSimWCHit();
	   newHit->SetTubeID(replicaNumber);
	   //newHit->SetTubeType(volumeName);//B. Quilain
	   newHit->SetEdep(energyDeposition); 
	   newHit->SetLogicalVolume(thePhysical->GetLogicalVolume());
	   
	   G4AffineTransform aTrans = theTouchable->GetHistory()->GetTopTransform();
	   newHit->SetRot(aTrans.NetRotation());
	   
	   aTrans.Invert();
	   newHit->SetPos(aTrans.NetTranslation());
	   // Set the hitMap value to the collection hit number
	   PMTHitMap[replicaNumber] = hitsCollection->insert( newHit );

	   
//	   (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddCosInc(cos_inc_f);
	   
	   (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPe(hitTime);
     (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddTrackID(trackID);
	   (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddParentID(parentSavedTrackID);
	   (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonStartTime(photonStartTime);
	   (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonStartEnergy(photonStartEnergy);
	   (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonEndEnergy(photonEndEnergy);
	   (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonStartPos(photonStartPos);
	   (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonEndPos(worldPosition);
	   (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonStartDir(photonStartDir);
	   (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonEndDir(worldDirection);
     (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonCreatorProcess(photonCreatorProcess);
	   
	   //     if ( particleDefinition != G4OpticalPhoton::OpticalPhotonDefinition() )
	   //       newHit->Print();
	     
	 }
       else {

//	 (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddCosInc(cos_inc_f);

	 (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPe(hitTime);
   (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddTrackID(trackID);
	 (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddParentID(parentSavedTrackID);
	 (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonStartTime(photonStartTime);
	 (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonStartEnergy(photonStartEnergy);
	 (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonEndEnergy(photonEndEnergy);
	 (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonStartPos(photonStartPos);
	 (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonEndPos(worldPosition);
	 (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonStartDir(photonStartDir);
	 (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonEndDir(worldDirection);
   (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonCreatorProcess(photonCreatorProcess);
	 
       }
     }
  }

  return true;
}

G4bool WCSimWCSD::ProcessHits_boundary(G4Step* aStep, G4TouchableHistory*)
{ 
  
  //JR EDIT BEGIN
  

/* UNCOMMENT IF THINGS GO WRONG
  // ===== JR DEBUG: unconditional Fill for optical photons (no creator filter) =====
  if (pmt_photon_tree) {
    G4Track* trk = aStep->GetTrack();
    if (trk && trk->GetDefinition()->GetPDGEncoding() == 0) { // optical photon
     
      const G4Event* gevt = G4RunManager::GetRunManager()->GetCurrentEvent();
      b_evt  = gevt ? gevt->GetEventID() : -1;

        b_par = trk->GetParentID();
        b_trk  = trk->GetTrackID();
        b_pdg      = trk->GetDefinition()->GetPDGEncoding();

        // emission (vertex) info
        b_ex_cm = trk->GetVertexPosition().x()/cm;
        b_ey_cm = trk->GetVertexPosition().y()/cm;
        b_ez_cm = trk->GetVertexPosition().z()/cm;
        b_et_ns = (trk->GetGlobalTime() - trk->GetLocalTime())/ns;
	
	// Emission direction (as stored at track vertex). This can be (0,0,0) sometimes.
	//G4ThreeVector edir = trk->GetVertexMomentumDirection();

	// --- Photon direction at PMT boundary (robust) ---
G4StepPoint* post = aStep->GetPostStepPoint();
G4ThreeVector edir = post->GetMomentumDirection();

b_edir_x = (float)edir.x();
b_edir_y = (float)edir.y();
b_edir_z = (float)edir.z();


        // detection point = where it hits the PMT boundary
        G4StepPoint* pre = aStep->GetPreStepPoint();
        b_dx_cm = pre->GetPosition().x()/cm;
        b_dy_cm = pre->GetPosition().y()/cm;
        b_dz_cm = pre->GetPosition().z()/cm;
        b_dt_ns = pre->GetGlobalTime()/ns;

        // photon energies
        b_Estart_MeV = trk->GetVertexKineticEnergy()/MeV;
        b_Eend_MeV   = trk->GetKineticEnergy()/MeV;

      if (b_Eend_MeV > 0) {
        const double hc_MeV_nm = 1.239841984e-3; // MeV*nm
        b_lambda_nm = hc_MeV_nm / b_Eend_MeV;
      } else {
        b_lambda_nm = -1.0;
      }

      // leave PMT id dummy for now
      b_pmt = 0;

      // creator string (if you have it)
      const G4VProcess* cp = trk->GetCreatorProcess();
      if (cp) {
        // if your branch is std::string:
        b_creator = cp->GetProcessName();
        // if your branch is char array, replace with:
        // std::snprintf(b_creator, sizeof(b_creator), "%s", cp->GetProcessName().c_str());
      } else {
        b_creator = "NULL";
      }

      pmt_photon_tree->Fill();
    }
  }
*/

  /*
  static int jr_calls = 0;
if (jr_calls < 5) {
  G4cout << "JR: ProcessHits_boundary called vol="
         << aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName()
         << " pdg=" << aStep->GetTrack()->GetDefinition()->GetPDGEncoding()
         << " parentID=" << aStep->GetTrack()->GetParentID()
         << G4endl;
  jr_calls++;
}
    // ---- JR DEBUG: minimal fill test ----
  if (pmt_photon_tree) {
    G4Track* trk = aStep->GetTrack();
    if (trk && trk->GetDefinition()->GetPDGEncoding() == 0) { // optical photon
      const G4VProcess* cp = trk->GetCreatorProcess();
      if (cp && cp->GetProcessName() == "Cerenkov") {

        const G4Event* gevt = G4RunManager::GetRunManager()->GetCurrentEvent();
        b_evt  = gevt ? gevt->GetEventID() : -1;

        b_par = trk->GetParentID();
        b_trk  = trk->GetTrackID();
        b_pdg      = trk->GetDefinition()->GetPDGEncoding();

        // emission (vertex) info
        b_ex_cm = trk->GetVertexPosition().x()/cm;
        b_ey_cm = trk->GetVertexPosition().y()/cm;
        b_ez_cm = trk->GetVertexPosition().z()/cm;
        b_et_ns = (trk->GetGlobalTime() - trk->GetLocalTime())/ns;

        // detection point = where it hits the PMT boundary
        G4StepPoint* pre = aStep->GetPreStepPoint();
        b_dx_cm = pre->GetPosition().x()/cm;
        b_dy_cm = pre->GetPosition().y()/cm;
        b_dz_cm = pre->GetPosition().z()/cm;
        b_dt_ns = pre->GetGlobalTime()/ns;

        // photon energies
        b_Estart_MeV = trk->GetVertexKineticEnergy()/MeV;
        b_Eend_MeV   = trk->GetKineticEnergy()/MeV;

        // wavelength from end energy (E in MeV)
        if (b_Eend_MeV > 0) {
          const double hc_MeV_nm = 1.239841984e-3; // MeV*nm
          b_lambda_nm = hc_MeV_nm / b_Eend_MeV;
        } else {
          b_lambda_nm = -1.0;
        }

        // TEMP: set PMT id to replicaNumber later; for now just 0
        b_pmt = 0;

        pmt_photon_tree->Fill();
      }
    }
  }
  // ---- END JR DEBUG ----
  */

//JR EDIT END
	
  G4StepPoint*       preStepPoint = aStep->GetPreStepPoint();
  G4StepPoint*       postStepPoint = aStep->GetPostStepPoint();
  G4TouchableHandle  theTouchable = postStepPoint->GetTouchableHandle();
  G4VPhysicalVolume* thePhysical  = theTouchable->GetVolume();

  //XQ 3/30/11 try to get the local position try to add the position and direction
  G4ThreeVector worldPosition  = postStepPoint->GetPosition();
  G4ThreeVector localPosition  = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPosition);
  G4ThreeVector worldDirection = preStepPoint->GetMomentumDirection();

  WCSimTrackInformation* trackinfo 
    = (WCSimTrackInformation*)(aStep->GetTrack()->GetUserInformation());

  G4int parentSavedTrackID = -1;
  G4float       photonStartTime;
  G4float       photonStartEnergy;
  G4float       photonEndEnergy;
  G4ThreeVector photonStartPos;
  G4ThreeVector photonStartDir;
  
  parentSavedTrackID   = aStep->GetTrack()->GetParentID();
  photonStartTime      = aStep->GetTrack()->GetGlobalTime() - aStep->GetTrack()->GetLocalTime(); // current time minus elapsed time of track
  photonStartEnergy    = aStep->GetTrack()->GetVertexKineticEnergy();
  photonEndEnergy      = aStep->GetTrack()->GetKineticEnergy();
  photonStartPos       = aStep->GetTrack()->GetVertexPosition();
  photonStartDir       = aStep->GetTrack()->GetVertexMomentumDirection();
 
  // Need to create a NONE string in case the Hit has no creatorProcess, such a Dark Noise Hit.
  const G4VProcess* process = aStep->GetTrack()->GetCreatorProcess();
  ProcessType_t photonCreatorProcess(kUnknownProcess);
  if (process) {
    photonCreatorProcess = WCSimEnumerations::ProcessTypeStringToEnum(process->GetProcessName());
  }

  G4int    trackID           = aStep->GetTrack()->GetTrackID();
  G4String volumeName        = thePhysical->GetName();
  
  G4double energyDeposition  = aStep->GetTotalEnergyDeposit();
  G4double hitTime           = postStepPoint->GetGlobalTime();

  G4ParticleDefinition *particleDefinition = 
    aStep->GetTrack()->GetDefinition();

  if ( particleDefinition != G4OpticalPhoton::OpticalPhotonDefinition())
    return false;

  G4String WCCollectionName;
  if(detectorElement=="tank") WCCollectionName = fdet->GetIDCollectionName();
  else if (detectorElement=="tankPMT2") WCCollectionName = fdet->GetIDCollectionName2();
  else if (detectorElement=="OD") WCCollectionName = fdet->GetODCollectionName();

  std::stringstream tubeTag;

  for (G4int i = theTouchable->GetHistoryDepth()-1 ; i >= 0; i--){
    tubeTag << ":" << theTouchable->GetVolume(i)->GetName();
    tubeTag << "-" << theTouchable->GetCopyNumber(i);
  }

  // Debug:
  // G4cout << "================================================" << G4endl;
  // G4cout << tubeTag.str() << G4endl;
  // G4cout << "================================================" << G4endl;

  // Get the tube ID from the tubeTag
  G4int replicaNumber;
  if(detectorElement=="tank") replicaNumber = WCSimDetectorConstruction::GetTubeID(tubeTag.str());
  else if(detectorElement=="tankPMT2") replicaNumber = WCSimDetectorConstruction::GetTubeID2(tubeTag.str());
  else if(detectorElement=="OD") replicaNumber = WCSimDetectorConstruction::GetODTubeID(tubeTag.str());
  else G4cout << "detectorElement not defined..." << G4endl;
  // Mark this photon as having hit a PMT in the pion_photons birth cache
  {
    const G4Track* currentTrack = aStep->GetTrack();
    auto it_birth = g_allPhotonBirth.find(currentTrack);
    if (it_birth != g_allPhotonBirth.end()) {
      it_birth->second.hit_pmt_id = (int)replicaNumber;
      
    }
  }
  G4double theta_angle = 0.;
  G4double effectiveAngularEfficiency = 0.;

  //XQ Add the wavelength there
  G4double  wavelength = (2.0*M_PI*197.3)/( aStep->GetTrack()->GetTotalEnergy()/eV);
  G4double ratio = 1.;
  G4double maxQE = 0.;
  G4double photonQE = 0.;
  if (fdet->GetPMT_QE_Method()==1 || fdet->GetPMT_QE_Method() == 4){
    photonQE = 1.1;
  }else if (fdet->GetPMT_QE_Method()==2){
    maxQE = fdet->GetPMTQE(WCCollectionName,wavelength,0,200,660,ratio);
    photonQE = fdet->GetPMTQE(volumeName, wavelength,1,200,660,ratio);
    photonQE = photonQE/maxQE;
  }else if (fdet->GetPMT_QE_Method() == 3){
    ratio = 1./(1.-0.25);
    photonQE = fdet->GetPMTQE(WCCollectionName, wavelength,1,200,660,ratio);
  }
  
  if (G4UniformRand() <= photonQE)
  {
    G4double local_x = localPosition.x();
    G4double local_y = localPosition.y();
    G4double local_z = localPosition.z();
    theta_angle = acos(fabs(local_z)/sqrt(pow(local_x,2)+pow(local_y,2)+pow(local_z,2)))/3.1415926*180.;
    effectiveAngularEfficiency = fdet->GetPMTCollectionEfficiency(theta_angle, volumeName);

    if (G4UniformRand() <= effectiveAngularEfficiency || fdet->UsePMT_Coll_Eff()==0)
    {
      //Retrieve the pointer to the appropriate hit collection. 
      //Since volumeName is the same as the SD name, this works. 
      G4SDManager* SDman = G4SDManager::GetSDMpointer();
      G4RunManager* Runman = G4RunManager::GetRunManager();
      G4int collectionID = SDman->GetCollectionID(volumeName);
      const G4Event* currentEvent = Runman->GetCurrentEvent();
      G4HCofThisEvent* HCofEvent = currentEvent->GetHCofThisEvent();
      hitsCollection = (WCSimWCHitsCollection*)(HCofEvent->GetHC(collectionID));

      // mark the track as having produced a hit
      if(!trackinfo)
          trackinfo = new WCSimTrackInformation();
      trackinfo->SetProducesHit(true);

      // Record this accepted PE (per-PE truth) for the pion_photons PE-bridge.
      // Same key/placement rationale as in ProcessHits: keyed on (eventID,
      // photon trackID), past both the QE and angular-efficiency rolls.  This
      // covers the mPMT / InteriorWCPMT boundary path.
      PionPhotonsTree_NotePE(currentEvent ? currentEvent->GetEventID() : -1,
                             trackID,
                             (int)replicaNumber,
                             (float)(hitTime / ns),
                             (float)wavelength);

      // If this tube hasn't been hit add it to the collection
      if (PMTHitMap[replicaNumber] == 0)
      {
        WCSimWCHit* newHit = new WCSimWCHit();
        newHit->SetTubeID(replicaNumber);
        newHit->SetEdep(energyDeposition); 
        newHit->SetLogicalVolume(thePhysical->GetLogicalVolume());
        
        G4AffineTransform aTrans = theTouchable->GetHistory()->GetTopTransform();
        newHit->SetRot(aTrans.NetRotation());
        
        aTrans.Invert();
        newHit->SetPos(aTrans.NetTranslation());
        // Set the hitMap value to the collection hit number
        PMTHitMap[replicaNumber] = hitsCollection->insert( newHit );
        (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPe(hitTime);
        (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddTrackID(trackID);
        (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddParentID(parentSavedTrackID);
        (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonStartTime(photonStartTime);
        (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonStartEnergy(photonStartEnergy);
        (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonEndEnergy(photonEndEnergy);
        (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonStartPos(photonStartPos);
        (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonEndPos(worldPosition);
        (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonStartDir(photonStartDir);
        (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonEndDir(worldDirection);
        (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonCreatorProcess(photonCreatorProcess);
      }
      else 
      {
        (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPe(hitTime);
        (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddTrackID(trackID);
        (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddParentID(parentSavedTrackID);
        (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonStartTime(photonStartTime);
        (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonStartEnergy(photonStartEnergy);
        (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonEndEnergy(photonEndEnergy);
        (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonStartPos(photonStartPos);
        (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonEndPos(worldPosition);
        (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonStartDir(photonStartDir);
        (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonEndDir(worldDirection);
        (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPhotonCreatorProcess(photonCreatorProcess);
        
      }
    }
  }

  return true;
}

void WCSimWCSD::EndOfEvent(G4HCofThisEvent* HCE)
{

 
  if (verboseLevel>1) 
  { 
    //Need to specify which collection in case multiple geometries are built

    G4String WCIDCollectionName = fdet->GetIDCollectionName();
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    G4int collectionID = SDman->GetCollectionID(WCIDCollectionName);
    hitsCollection = (WCSimWCHitsCollection*)HCE->GetHC(collectionID);

    G4int numHits = hitsCollection->entries();

    G4cout << "There are " << numHits << " hits in the "<<detectorElement<<" : "<< G4endl;
    for (G4int i=0; i < numHits; i++) {
      G4cout<<"ihit ID = "<<i<<G4endl;
      (*hitsCollection)[i]->Print();
    }

    //Added by B. Quilain for the hybrid version
    G4cout<<"Tube hit list finalized"<<G4endl;
    G4cout<<"Geometry is hybrid = "<<fdet->GetHybridPMT()<<G4endl;
    if (fdet->GetHybridPMT()) { //TD fdet->2019/07/13
      G4String WCIDCollectionName2 = fdet->GetIDCollectionName2();
      G4int collectionID2 = SDman->GetCollectionID(WCIDCollectionName2);
      hitsCollection2 = (WCSimWCHitsCollection*)HCE->GetHC(collectionID2);
    
      G4int numHits2 = hitsCollection2->entries();

      G4cout << "There are " << numHits2 << " tubes hit in the WC: " << G4endl;
      for (G4int i=0; i < numHits2; i++){
	G4cout<<"ihit ID = "<<i<<G4endl;
	(*hitsCollection2)[i]->Print();
      }
    }
      /*
    {
      if(abs((*hitsCollection)[i]->GetTubeID() - 1584)  < 5){
	  G4cout << (*hitsCollection)[i]->GetTubeID() << G4endl;
	  (*hitsCollection)[i]->Print();
      }
      }*/

    /* Detailed debug:
    G4cout << "Through mPMTLV " << WCSimSteppingAction::n_photons_through_mPMTLV << G4endl;
    G4cout << "Through Acrylic " << WCSimSteppingAction::n_photons_through_acrylic << G4endl;
    G4cout << "Through Gel " << WCSimSteppingAction::n_photons_through_gel << G4endl;
    G4cout << "On Blacksheet " << WCSimSteppingAction::n_photons_on_blacksheet << G4endl;
    G4cout << "On small PMT " << WCSimSteppingAction::n_photons_on_smallPMT << G4endl;
    */
  } 
}

