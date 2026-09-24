#include "WCSimSecondaryTrackTree.hh"
#include "TFile.h"
#include "TTree.h"
#include <cmath>

// ---- global storage ----

TTree* secondary_tracks_tree = nullptr;

int   st_evt        = -1;
int   st_trk        = -1;
int   st_pdg        = 0;
int   st_parent_trk = -1;
int   st_parent_pdg = 0;
std::string st_creator     = "";
int   st_src_cat    = 9;

float st_x0_cm=0, st_y0_cm=0, st_z0_cm=0, st_t0_ns=0, st_ke0_MeV=0;
float st_x1_cm=0, st_y1_cm=0, st_z1_cm=0, st_t1_ns=0, st_ke1_MeV=0;
std::string st_end_process = "";
int         st_src_trk     = -1;
int         st_src_pdg     = 0;
float       st_src_ke0_MeV = 0.f;
std::string st_src_creator = "";
int         st_anc_trk     = -1;
int         st_anc_pdg     = 0;
float       st_anc_ke0_MeV = 0.f;
std::string st_anc_creator = "";
float st_px0=0, st_py0=0, st_pz0=0;
float st_px1=0, st_py1=0, st_pz1=0;
std::unordered_map<int, SecTrackBuffer> g_secTrackBuffer;


// ----------------------------------------------------------------
bool SecTrack_IsInteresting(int pdg, const std::string& creatorProcess)
{
  const int a = std::abs(pdg);

  // Charged pions (π±) and neutral pion (π⁰)
  if (a == 211 || pdg == 111) return true;

  // Muons
  if (a == 13) return true;

  // Electrons/positrons from muon decay (Michel) or gamma conversion
  if (a == 11 && (creatorProcess == "Decay" || creatorProcess == "conv"))
    return true;

  // Prompt gammas from π⁰ decay (π⁰ → γγ)
  if (pdg == 22 && creatorProcess == "Decay") return true;

  // Nuclear fragments: proton, neutron, heavy ions (Geant4 ion encoding)
  if (a == 2212 || a == 2112 || a > 1000000000) return true;

  return false;
}


// ----------------------------------------------------------------
void SecondaryTracksTree_Book(TFile* fout)
{
  if (secondary_tracks_tree) return;
  if (!fout)                  return;
  fout->cd();

  secondary_tracks_tree = new TTree("secondary_tracks",
      "Interesting secondaries: pi, mu, Michel e, pi0, nuclear fragments");

  secondary_tracks_tree->Branch("evt",        &st_evt,        "evt/I");
  secondary_tracks_tree->Branch("trk",        &st_trk,        "trk/I");
  secondary_tracks_tree->Branch("pdg",        &st_pdg,        "pdg/I");
  secondary_tracks_tree->Branch("parent_trk", &st_parent_trk, "parent_trk/I");
  secondary_tracks_tree->Branch("parent_pdg", &st_parent_pdg, "parent_pdg/I");
  secondary_tracks_tree->Branch("creator",    &st_creator);
  secondary_tracks_tree->Branch("src_cat",    &st_src_cat,    "src_cat/I");

  secondary_tracks_tree->Branch("x0_cm",   &st_x0_cm,   "x0_cm/F");
  secondary_tracks_tree->Branch("y0_cm",   &st_y0_cm,   "y0_cm/F");
  secondary_tracks_tree->Branch("z0_cm",   &st_z0_cm,   "z0_cm/F");
  secondary_tracks_tree->Branch("t0_ns",   &st_t0_ns,   "t0_ns/F");
  secondary_tracks_tree->Branch("ke0_MeV", &st_ke0_MeV, "ke0_MeV/F");

  secondary_tracks_tree->Branch("x1_cm",   &st_x1_cm,   "x1_cm/F");
  secondary_tracks_tree->Branch("y1_cm",   &st_y1_cm,   "y1_cm/F");
  secondary_tracks_tree->Branch("z1_cm",   &st_z1_cm,   "z1_cm/F");
  secondary_tracks_tree->Branch("t1_ns",   &st_t1_ns,   "t1_ns/F");
  secondary_tracks_tree->Branch("ke1_MeV", &st_ke1_MeV, "ke1_MeV/F");

  secondary_tracks_tree->Branch("end_process", &st_end_process);

  secondary_tracks_tree->Branch("src_trk",     &st_src_trk,     "src_trk/I");
  secondary_tracks_tree->Branch("src_pdg",     &st_src_pdg,     "src_pdg/I");
  secondary_tracks_tree->Branch("src_ke0_MeV", &st_src_ke0_MeV, "src_ke0_MeV/F");
  secondary_tracks_tree->Branch("src_creator", &st_src_creator);

  secondary_tracks_tree->Branch("anc_trk",     &st_anc_trk,     "anc_trk/I");
  secondary_tracks_tree->Branch("anc_pdg",     &st_anc_pdg,     "anc_pdg/I");
  secondary_tracks_tree->Branch("anc_ke0_MeV", &st_anc_ke0_MeV, "anc_ke0_MeV/F");
  secondary_tracks_tree->Branch("anc_creator", &st_anc_creator);
  secondary_tracks_tree->Branch("px0", &st_px0, "px0/F");
  secondary_tracks_tree->Branch("py0", &st_py0, "py0/F");
  secondary_tracks_tree->Branch("pz0", &st_pz0, "pz0/F");
  secondary_tracks_tree->Branch("px1", &st_px1, "px1/F");
  secondary_tracks_tree->Branch("py1", &st_py1, "py1/F");
  secondary_tracks_tree->Branch("pz1", &st_pz1, "pz1/F");
}


// ----------------------------------------------------------------
void SecondaryTracksTree_CacheBirth(int trackID, const SecTrackBuffer& buf)
{
  g_secTrackBuffer[trackID] = buf;
}


// ----------------------------------------------------------------
void SecondaryTracksTree_FillEnd(int         trackID,
                                  float       x1_cm,
                                  float       y1_cm,
                                  float       z1_cm,
                                  float       t1_ns,
                                  float       ke1_MeV,
                                  float px1, float py1, float pz1,
                                  const std::string& end_process)
{
  auto it = g_secTrackBuffer.find(trackID);
  if (it == g_secTrackBuffer.end()) return;

  const SecTrackBuffer& buf = it->second;

  st_evt        = buf.evt;
  st_trk        = buf.trk;
  st_pdg        = buf.pdg;
  st_parent_trk = buf.parent_trk;
  st_parent_pdg = buf.parent_pdg;
  st_creator    = buf.creator;
  st_src_cat    = buf.src_cat;

  st_x0_cm   = buf.x0_cm;
  st_y0_cm   = buf.y0_cm;
  st_z0_cm   = buf.z0_cm;
  st_t0_ns   = buf.t0_ns;
  st_ke0_MeV = buf.ke0_MeV;

  st_x1_cm    = x1_cm;
  st_y1_cm    = y1_cm;
  st_z1_cm    = z1_cm;
  st_t1_ns    = t1_ns;
  st_ke1_MeV  = ke1_MeV;
  st_end_process = end_process;
  st_src_trk     = buf.src_trk;
  st_src_pdg     = buf.src_pdg;
  st_src_ke0_MeV = buf.src_ke0_MeV;
  st_src_creator = buf.src_creator;
  st_anc_trk     = buf.anc_trk;
  st_anc_pdg     = buf.anc_pdg;
  st_anc_ke0_MeV = buf.anc_ke0_MeV;
  st_anc_creator = buf.anc_creator;
  st_px0 = buf.px0;
  st_py0 = buf.py0;
  st_pz0 = buf.pz0;
  st_px1 = px1;
  st_py1 = py1;
  st_pz1 = pz1;
  if (secondary_tracks_tree) secondary_tracks_tree->Fill();

  g_secTrackBuffer.erase(it);
}


// ----------------------------------------------------------------
void SecondaryTracksTree_Write()
{
  if (!secondary_tracks_tree) return;
  secondary_tracks_tree->Write("", TObject::kOverwrite);
}


// ----------------------------------------------------------------
void SecondaryTracksTree_ClearBuffer()
{
  g_secTrackBuffer.clear();
}
