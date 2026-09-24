#include "WCSimPionSourceTree.hh"
#include "TFile.h"
#include "TTree.h"
#include <unordered_map>

// ---- global storage ----

TTree* pion_photons_tree = nullptr;

int   pp_evt = -1;
int   pp_trk = -1;

float pp_ex_cm  = 0.f, pp_ey_cm  = 0.f, pp_ez_cm  = 0.f;
float pp_et_ns  = 0.f;
float pp_edir_x = 0.f, pp_edir_y = 0.f, pp_edir_z = 0.f;
float pp_lambda_nm = 0.f;
std::string pp_creator = "";

int         pp_src_trk      = -1;
int         pp_src_pdg      = 0;
float       pp_src_ke_MeV   = 0.f;
float       pp_src_ke0_MeV  = 0.f;
int         pp_src_step     = -1;
std::string pp_src_creator  = "";

int         pp_anc_trk      = -1;
int         pp_anc_pdg      = 0;
float       pp_anc_ke0_MeV  = 0.f;
std::string pp_anc_creator  = "";

int pp_src_cat = 9;

int         pp_hit_pmt      = -1;
std::string pp_end_boundary = "";
std::string pp_end_process  = "";

// NEW: PE-bridge globals
int   pp_made_pe       = 0;
int   pp_pe_pmt        = -1;
float pp_pe_time_ns    = 0.f;
float pp_pe_lambda_nm  = 0.f;


// ----------------------------------------------------------------
void PionPhotonsTree_Book(TFile* fout)
{
  if (pion_photons_tree) return;
  if (!fout)             return;

  fout->cd();

  pion_photons_tree = new TTree("pion_photons",
      "Cherenkov photons with full ancestry + PE acceptance (pion analysis)");

  pion_photons_tree->Branch("evt", &pp_evt, "evt/I");
  pion_photons_tree->Branch("trk", &pp_trk, "trk/I");

  pion_photons_tree->Branch("ex_cm",     &pp_ex_cm,     "ex_cm/F");
  pion_photons_tree->Branch("ey_cm",     &pp_ey_cm,     "ey_cm/F");
  pion_photons_tree->Branch("ez_cm",     &pp_ez_cm,     "ez_cm/F");
  pion_photons_tree->Branch("et_ns",     &pp_et_ns,     "et_ns/F");
  pion_photons_tree->Branch("edir_x",    &pp_edir_x,    "edir_x/F");
  pion_photons_tree->Branch("edir_y",    &pp_edir_y,    "edir_y/F");
  pion_photons_tree->Branch("edir_z",    &pp_edir_z,    "edir_z/F");
  pion_photons_tree->Branch("lambda_nm", &pp_lambda_nm, "lambda_nm/F");
  pion_photons_tree->Branch("creator",   &pp_creator);

  pion_photons_tree->Branch("src_trk",     &pp_src_trk,     "src_trk/I");
  pion_photons_tree->Branch("src_pdg",     &pp_src_pdg,     "src_pdg/I");
  pion_photons_tree->Branch("src_ke_MeV",  &pp_src_ke_MeV,  "src_ke_MeV/F");
  pion_photons_tree->Branch("src_ke0_MeV", &pp_src_ke0_MeV, "src_ke0_MeV/F");
  pion_photons_tree->Branch("src_step",    &pp_src_step,    "src_step/I");
  pion_photons_tree->Branch("src_creator", &pp_src_creator);

  pion_photons_tree->Branch("anc_trk",     &pp_anc_trk,     "anc_trk/I");
  pion_photons_tree->Branch("anc_pdg",     &pp_anc_pdg,     "anc_pdg/I");
  pion_photons_tree->Branch("anc_ke0_MeV", &pp_anc_ke0_MeV, "anc_ke0_MeV/F");
  pion_photons_tree->Branch("anc_creator", &pp_anc_creator);

  pion_photons_tree->Branch("src_cat", &pp_src_cat, "src_cat/I");

  pion_photons_tree->Branch("hit_pmt",      &pp_hit_pmt,      "hit_pmt/I");
  pion_photons_tree->Branch("end_boundary", &pp_end_boundary);
  pion_photons_tree->Branch("end_process",  &pp_end_process);

  // NEW PE-bridge branches
  pion_photons_tree->Branch("made_pe",      &pp_made_pe,      "made_pe/I");
  pion_photons_tree->Branch("pe_pmt",       &pp_pe_pmt,       "pe_pmt/I");
  pion_photons_tree->Branch("pe_time_ns",   &pp_pe_time_ns,   "pe_time_ns/F");
  pion_photons_tree->Branch("pe_lambda_nm", &pp_pe_lambda_nm, "pe_lambda_nm/F");
}

void PionPhotonsTree_Fill()
{
  if (pion_photons_tree) pion_photons_tree->Fill();
}

void PionPhotonsTree_Write()
{
  if (!pion_photons_tree) return;
  pion_photons_tree->Write("", TObject::kOverwrite);
}


// ================================================================
//  PE-bridge implementation (per-event map; auto-reset on new evt)
// ================================================================
namespace {
  struct PionPhotonPEInfo {
    int   pmt       = -1;
    float time_ns   = 0.f;
    float lambda_nm = 0.f;
  };

  static std::unordered_map<long long, PionPhotonPEInfo> g_pionPhotonPE;
  static int g_pionPhotonPE_evt = -999999;

  long long MakeKey(int evt, int trk)
  {
    return (static_cast<long long>(evt) << 32) ^
           static_cast<unsigned int>(trk);
  }
}

void PionPhotonsTree_ResetEvent(int evt)
{
  if (evt == g_pionPhotonPE_evt) return;
  g_pionPhotonPE_evt = evt;
  g_pionPhotonPE.clear();
}

void PionPhotonsTree_NotePE(int evt, int photon_trk,
                             int pe_pmt, float pe_time_ns, float pe_lambda_nm)
{
  PionPhotonsTree_ResetEvent(evt);
  PionPhotonPEInfo info;
  info.pmt       = pe_pmt;
  info.time_ns   = pe_time_ns;
  info.lambda_nm = pe_lambda_nm;
  g_pionPhotonPE[MakeKey(evt, photon_trk)] = info;
}

bool PionPhotonsTree_GetPE(int evt, int photon_trk,
                            int& pe_pmt, float& pe_time_ns, float& pe_lambda_nm)
{
  PionPhotonsTree_ResetEvent(evt);
  auto it = g_pionPhotonPE.find(MakeKey(evt, photon_trk));
  if (it == g_pionPhotonPE.end()) return false;

  pe_pmt       = it->second.pmt;
  pe_time_ns   = it->second.time_ns;
  pe_lambda_nm = it->second.lambda_nm;
  return true;
}