#include "WCSimIonizingStepsTree.hh"
#include "TFile.h"
#include "TTree.h"

TTree* ionizing_steps_tree = nullptr;

int         is_evt          = -1;
int         is_trk          = -1;
int         is_pdg          = 0;
int         is_parent_trk   = -1;
int         is_src_cat      = 9;
int         is_step         = 0;

float       is_x_cm         = 0.f;
float       is_y_cm         = 0.f;
float       is_z_cm         = 0.f;
float       is_t_ns         = 0.f;

float       is_ke_pre_MeV   = 0.f;
float       is_ke_post_MeV  = 0.f;
float       is_edep_MeV     = 0.f;

float       is_step_len_cm  = 0.f;
float       is_track_len_cm = 0.f;

float       is_beta_pre     = 0.f;
float       is_beta_post    = 0.f;
int         is_above_ckv    = 0;

int         is_n_cherenkov  = 0;
std::string is_process      = "";


void IonizingStepsTree_Book(TFile* fout)
{
  if (ionizing_steps_tree) return;
  if (!fout)               return;
  fout->cd();

  ionizing_steps_tree = new TTree("ionizing_steps",
      "Per-step record of the primary and its charged daughters");

  ionizing_steps_tree->Branch("evt",          &is_evt,          "evt/I");
  ionizing_steps_tree->Branch("trk",          &is_trk,          "trk/I");
  ionizing_steps_tree->Branch("pdg",          &is_pdg,          "pdg/I");
  ionizing_steps_tree->Branch("parent_trk",   &is_parent_trk,   "parent_trk/I");
  ionizing_steps_tree->Branch("src_cat",      &is_src_cat,      "src_cat/I");
  ionizing_steps_tree->Branch("step",         &is_step,         "step/I");

  ionizing_steps_tree->Branch("x_cm",         &is_x_cm,         "x_cm/F");
  ionizing_steps_tree->Branch("y_cm",         &is_y_cm,         "y_cm/F");
  ionizing_steps_tree->Branch("z_cm",         &is_z_cm,         "z_cm/F");
  ionizing_steps_tree->Branch("t_ns",         &is_t_ns,         "t_ns/F");

  ionizing_steps_tree->Branch("ke_pre_MeV",   &is_ke_pre_MeV,   "ke_pre_MeV/F");
  ionizing_steps_tree->Branch("ke_post_MeV",  &is_ke_post_MeV,  "ke_post_MeV/F");
  ionizing_steps_tree->Branch("edep_MeV",     &is_edep_MeV,     "edep_MeV/F");

  ionizing_steps_tree->Branch("step_len_cm",  &is_step_len_cm,  "step_len_cm/F");
  ionizing_steps_tree->Branch("track_len_cm", &is_track_len_cm, "track_len_cm/F");

  ionizing_steps_tree->Branch("beta_pre",     &is_beta_pre,     "beta_pre/F");
  ionizing_steps_tree->Branch("beta_post",    &is_beta_post,    "beta_post/F");
  ionizing_steps_tree->Branch("above_ckv",    &is_above_ckv,    "above_ckv/I");

  ionizing_steps_tree->Branch("n_cherenkov",  &is_n_cherenkov,  "n_cherenkov/I");
  ionizing_steps_tree->Branch("process",      &is_process);
}


void IonizingStepsTree_Fill()
{
  if (ionizing_steps_tree) ionizing_steps_tree->Fill();
}


void IonizingStepsTree_Write()
{
  if (!ionizing_steps_tree) return;
  ionizing_steps_tree->Write("", TObject::kOverwrite);
}
