#include "WCSimPionStepsTree.hh"
#include "TFile.h"
#include "TTree.h"

// ---- global storage ----

TTree* pion_steps_tree = nullptr;

int         ps_evt           = -1;
int         ps_trk           = -1;
int         ps_pdg           = 0;
int         ps_src_cat       = 9;

float       ps_x_cm          = 0.f;
float       ps_y_cm          = 0.f;
float       ps_z_cm          = 0.f;
float       ps_t_ns          = 0.f;

float       ps_ke_pre_MeV    = 0.f;
float       ps_ke_post_MeV   = 0.f;

std::string ps_process       = "";

int         ps_n_secondaries = 0;

float ps_dir_pre_x=0,  ps_dir_pre_y=0,  ps_dir_pre_z=0;
float ps_dir_post_x=0, ps_dir_post_y=0, ps_dir_post_z=0;
// ----------------------------------------------------------------
void PionStepsTree_Book(TFile* fout)
{
  if (pion_steps_tree) return;
  if (!fout)           return;
  fout->cd();

  pion_steps_tree = new TTree("pion_steps",
      "Hadronic interaction / decay steps of pion tracks");

  pion_steps_tree->Branch("evt",           &ps_evt,           "evt/I");
  pion_steps_tree->Branch("trk",           &ps_trk,           "trk/I");
  pion_steps_tree->Branch("pdg",           &ps_pdg,           "pdg/I");
  pion_steps_tree->Branch("src_cat",       &ps_src_cat,       "src_cat/I");

  pion_steps_tree->Branch("x_cm",          &ps_x_cm,          "x_cm/F");
  pion_steps_tree->Branch("y_cm",          &ps_y_cm,          "y_cm/F");
  pion_steps_tree->Branch("z_cm",          &ps_z_cm,          "z_cm/F");
  pion_steps_tree->Branch("t_ns",          &ps_t_ns,          "t_ns/F");

  pion_steps_tree->Branch("ke_pre_MeV",    &ps_ke_pre_MeV,    "ke_pre_MeV/F");
  pion_steps_tree->Branch("ke_post_MeV",   &ps_ke_post_MeV,   "ke_post_MeV/F");

  pion_steps_tree->Branch("process",       &ps_process);

  pion_steps_tree->Branch("n_secondaries", &ps_n_secondaries, "n_secondaries/I");
  pion_steps_tree->Branch("dir_pre_x", &ps_dir_pre_x, "dir_pre_x/F");
  pion_steps_tree->Branch("dir_pre_y", &ps_dir_pre_y, "dir_pre_y/F");
  pion_steps_tree->Branch("dir_pre_z", &ps_dir_pre_z, "dir_pre_z/F");
  pion_steps_tree->Branch("dir_post_x", &ps_dir_post_x, "dir_post_x/F");
  pion_steps_tree->Branch("dir_post_y", &ps_dir_post_y, "dir_post_y/F");
  pion_steps_tree->Branch("dir_post_z", &ps_dir_post_z, "dir_post_z/F");
}


// ----------------------------------------------------------------
void PionStepsTree_Fill()
{
  if (pion_steps_tree) pion_steps_tree->Fill();
}


// ----------------------------------------------------------------
void PionStepsTree_Write()
{
  if (!pion_steps_tree) return;
  pion_steps_tree->Write("", TObject::kOverwrite);
}
