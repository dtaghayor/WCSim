#include "WCSimAncestryMap.hh"
#include <cmath>

// Global definition (declared extern in the header)
std::unordered_map<int, AncestryInfo> g_trackAncestry;

// trackID -> true terminating process name (see header).
std::unordered_map<int, std::string> g_trackEndProc;

// Per-event index shared by all custom pion-analysis trees.
int g_wcsim_evt = -1;

// ----------------------------------------------------------------
const char* SourceCategoryName(SourceCategory cat)
{
  switch (cat) {
    case kPrimaryPion:  return "PrimaryPion";
    case kSecPionPlus:  return "SecPionPlus";
    case kSecPionMinus: return "SecPionMinus";
    case kPionZero:     return "PionZero";
    case kMuon:         return "Muon";
    case kMichelElec:   return "MichelElectron";
    case kDeltaRay:     return "DeltaRay";
    case kNuclear:      return "Nuclear";
    case kGamma:        return "Gamma";
    default:            return "Other";
  }
}

// ----------------------------------------------------------------
void AncestryMap_Register(int         trackID,
                           int         pdg,
                           int         parentID,
                           float       ke0_MeV,
                           float       x0_cm,
                           float       y0_cm,
                           float       z0_cm,
                           float       t0_ns,
                           const std::string& creatorProcess)
{
  AncestryInfo info;
  info.pdg            = pdg;
  info.parentID       = parentID;
  info.ke0_MeV        = ke0_MeV;
  info.x0_cm          = x0_cm;
  info.y0_cm          = y0_cm;
  info.z0_cm          = z0_cm;
  info.t0_ns          = t0_ns;
  info.creatorProcess = creatorProcess;
  g_trackAncestry[trackID] = info;
}

// ----------------------------------------------------------------
void AncestryMap_Clear()
{
  g_trackAncestry.clear();
  g_trackEndProc.clear();
}

// ----------------------------------------------------------------
// Internal helpers
// ----------------------------------------------------------------

static bool IsNuclearFragment(int pdg)
{
  const int a = std::abs(pdg);
  // Geant4 encodes heavy ions as 10LZZZAAAI
  // Proton = 2212, neutron = 2112, alpha = 1000020040, etc.
  return (a == 2212 || a == 2112 || a > 1000000000);
}

static bool IsIonisationCreated(const std::string& proc)
{
  // Delta-rays (knock-on electrons) are created by the ionisation
  // processes of heavier particles, or by electron ionisation itself.
  return (proc == "hIoni"   ||   // hadron ionisation
          proc == "muIoni"  ||   // muon  ionisation
          proc == "eIoni"   ||   // electron ionisation
          proc == "ionIoni" ||   // ion   ionisation
          proc == "hBrem"   ||   // hadron bremsstrahlung (e from brem shower)
          proc == "muBrem");     // muon  bremsstrahlung
}

// ----------------------------------------------------------------
SourceCategory AncestryMap_ClassifyPhoton(
    int          parentID,
    int&         outSrcTrk,
    int&         outSrcPDG,
    float&       outSrcKE0,
    std::string& outSrcCreator,
    int&         outAncTrk,
    int&         outAncPDG,
    float&       outAncKE0,
    std::string& outAncCreator)
{
  // --- Initialise outputs to safe defaults ---
  outSrcTrk = parentID;  outSrcPDG = 0;  outSrcKE0 = 0.f;  outSrcCreator = "UNKNOWN";
  outAncTrk = -1;        outAncPDG = 0;  outAncKE0 = 0.f;  outAncCreator = "UNKNOWN";

  // --- Fill direct-parent (src) info ---
  {
    auto it = g_trackAncestry.find(parentID);
    if (it == g_trackAncestry.end()) return kOther;
    outSrcPDG     = it->second.pdg;
    outSrcKE0     = it->second.ke0_MeV;
    outSrcCreator = it->second.creatorProcess;
  }

  // --- Walk the ancestry chain to find the interesting ancestor ---
  //
  // Strategy:
  //   • Stop immediately on π±, π⁰, μ±, Michel e±, nuclear fragment.
  //   • For a delta-ray e± (ionisation-created): stop here too.
  //   • For γ and conv-created e± (from γ→e+e-): keep walking upward
  //     to find what produced the γ (e.g. π⁰ decay).  Cache the γ as
  //     a fallback in case we reach the top without a better match.
  //   • At the top (parentID==0) with no match, return kOther.

  int   gammaCandTrk = -1;
  float gammaCandKE0 = 0.f;
  std::string gammaCandCreator;

  int cur = parentID;

  while (cur > 0) {
    auto it = g_trackAncestry.find(cur);
    if (it == g_trackAncestry.end()) break;

    const AncestryInfo& info = it->second;
    const int pdg    = info.pdg;
    const int abspdg = std::abs(pdg);

    // ---- π± ---------------------------------------------------------
    if (abspdg == 211) {
      outAncTrk = cur;  outAncPDG = pdg;
      outAncKE0 = info.ke0_MeV;  outAncCreator = info.creatorProcess;
      if (info.parentID == 0) return kPrimaryPion;
      return (pdg > 0) ? kSecPionPlus : kSecPionMinus;
    }

    // ---- π⁰ ---------------------------------------------------------
    if (pdg == 111) {
      outAncTrk = cur;  outAncPDG = pdg;
      outAncKE0 = info.ke0_MeV;  outAncCreator = info.creatorProcess;
      return kPionZero;
    }

    // ---- μ± ---------------------------------------------------------
    if (abspdg == 13) {
      outAncTrk = cur;  outAncPDG = pdg;
      outAncKE0 = info.ke0_MeV;  outAncCreator = info.creatorProcess;
      return kMuon;
    }

    // ---- Michel e± (from muon decay) --------------------------------
    if (abspdg == 11 && info.creatorProcess == "Decay") {
      // Only a true Michel electron if its immediate parent was a muon
      auto parentIt = g_trackAncestry.find(info.parentID);
      if (parentIt != g_trackAncestry.end() &&
          std::abs(parentIt->second.pdg) == 13) {
        outAncTrk = cur;  outAncPDG = pdg;
        outAncKE0 = info.ke0_MeV;  outAncCreator = info.creatorProcess;
        return kMichelElec;
      }
      // else: decay electron from pion or other source — keep walking up the chain
    }

    // ---- Delta-ray (ionisation-created e±) --------------------------
    if (abspdg == 11 && IsIonisationCreated(info.creatorProcess)) {
      outAncTrk = cur;  outAncPDG = pdg;
      outAncKE0 = info.ke0_MeV;  outAncCreator = info.creatorProcess;
      return kDeltaRay;
    }

    // ---- Nuclear fragment (p, n, heavy ion) -------------------------
    if (IsNuclearFragment(pdg)) {
      outAncTrk = cur;  outAncPDG = pdg;
      outAncKE0 = info.ke0_MeV;  outAncCreator = info.creatorProcess;
      return kNuclear;
    }

    // ---- γ: keep walking to find the γ's origin --------------------
    //  e.g.  γ (from π⁰ decay) → e+ → Cherenkov photon
    //  We want to report kPionZero, not kGamma, so we keep walking.
    //  Cache the γ as a fallback in case the chain ends here.
    if (pdg == 22) {
      if (gammaCandTrk < 0) {   // only save the first (lowest) γ in chain
        gammaCandTrk     = cur;
        gammaCandKE0     = info.ke0_MeV;
        gammaCandCreator = info.creatorProcess;
      }
      // continue walking upward
    }

    // ---- e± from γ conversion: keep walking (don't stop) -----------
    //  conv-created e± always have a γ parent, which may be from π⁰.
    //  We will pick up the γ above and continue to π⁰ or wherever.

    cur = info.parentID;
  }

  // Exhausted the chain.  Return best available result.
  if (gammaCandTrk >= 0) {
    outAncTrk     = gammaCandTrk;
    outAncPDG     = 22;
    outAncKE0     = gammaCandKE0;
    outAncCreator = gammaCandCreator;
    return kGamma;
  }

  return kOther;
}



// ================================================================
//  AncestryMap_ClassifyTrack implementation
//
//  Add this function to the bottom of WCSimAncestryMap.cc
//  (before the final closing brace if there is one, otherwise
//  just append it to the end of the file).
// ================================================================

SourceCategory AncestryMap_ClassifyTrack(
    int          trackID,
    int          pdg,
    int          parentID,
    const std::string& creatorProcess,
    int&         outSrcTrk,        // immediate parent track ID
    int&         outSrcPDG,        // immediate parent PDG
    float&       outSrcKE0,        // immediate parent birth KE [MeV]
    std::string& outSrcCreator,    // process that created the immediate parent
    int&         outAncTrk,        // interesting ancestor track ID
    int&         outAncPDG,        // interesting ancestor PDG
    float&       outAncKE0,        // interesting ancestor birth KE [MeV]
    std::string& outAncCreator)   // process that created interesting ancestor
{
  // Safe defaults
  outAncTrk = trackID;  outAncPDG = pdg;
  outAncKE0 = 0.f;      outAncCreator = creatorProcess;
  // Fill src* from immediate parent
  outSrcTrk     = parentID;
  outSrcPDG     = 0;
  outSrcKE0     = 0.f;
  outSrcCreator = "";
  {
    auto it = g_trackAncestry.find(parentID);
    if (it != g_trackAncestry.end()) {
      outSrcPDG     = it->second.pdg;
      outSrcKE0     = it->second.ke0_MeV;
      outSrcCreator = it->second.creatorProcess;
    }
  }
  // Fill outAncKE0 from ancestry map if available
  {
    auto it = g_trackAncestry.find(trackID);
    if (it != g_trackAncestry.end()) outAncKE0 = it->second.ke0_MeV;
  }

  const int abspdg = std::abs(pdg);

  // ----------------------------------------------------------------
  // PRIMARY PION: parentID == 0 means this IS the primary vertex
  // particle.  For pion beams this is always a π±.
  // ----------------------------------------------------------------
  if (parentID == 0) {
    if (abspdg == 211) return kPrimaryPion;
    if (abspdg == 13)  return kMuon;        // e.g. if simulating muon beam
    if (abspdg == 11)  return kOther;  // edge case: primary e±
    if (pdg == 22)     return kGamma;
    return kOther;
  }

  // ----------------------------------------------------------------
  // SECONDARY π±
  // ----------------------------------------------------------------
  if (abspdg == 211) {
    return (pdg > 0) ? kSecPionPlus : kSecPionMinus;
  }

  // ----------------------------------------------------------------
  // π⁰  (always secondary in pion beam; π⁰ has no charge so it
  //       cannot be the beam particle)
  // ----------------------------------------------------------------
  if (pdg == 111) {
    return kPionZero;
  }

  // ----------------------------------------------------------------
  // μ±  — from π decay (dominant) or any other source
  // ----------------------------------------------------------------
  if (abspdg == 13) {
    return kMuon;
  }

  // ----------------------------------------------------------------
  // e± from muon decay → Michel electron
  // Require the immediate parent to be a μ± to avoid mis-tagging
  // Dalitz electrons (π⁰ → e⁺e⁻γ) or pion leptonic decay electrons.
  // ----------------------------------------------------------------
  if (abspdg == 11 && creatorProcess == "Decay") {
    auto it = g_trackAncestry.find(parentID);
    if (it != g_trackAncestry.end() && std::abs(it->second.pdg) == 13) {
      // Report ancestor as the muon parent, not this electron
      outAncTrk     = parentID;
      outAncPDG     = it->second.pdg;
      outAncKE0     = it->second.ke0_MeV;
      outAncCreator = it->second.creatorProcess;
      return kMichelElec;
    }
    // Decay electron from non-muon (pion leptonic, Dalitz): fall through
    // to walk the ancestry chain and find the real ancestor below.
  }

  // ----------------------------------------------------------------
  // e± from ionisation → delta-ray
  // Report the delta-ray itself as the ancestor (it IS the interesting
  // object; its parent is just the track that knocked it out).
  // ----------------------------------------------------------------
  if (abspdg == 11 && IsIonisationCreated(creatorProcess)) {
    return kDeltaRay;
  }

  // ----------------------------------------------------------------
  // γ  — prompt (from π⁰→γγ, nCapture, bremsstrahlung, …)
  //  Walk up to find whether this γ descends from a π⁰.
  // ----------------------------------------------------------------
  if (pdg == 22) {
    int cur = parentID;
    while (cur > 0) {
      auto it = g_trackAncestry.find(cur);
      if (it == g_trackAncestry.end()) break;
      const int apd = std::abs(it->second.pdg);
      if (it->second.pdg == 111) {          // π⁰ ancestor found
        outAncTrk     = cur;
        outAncPDG     = 111;
        outAncKE0     = it->second.ke0_MeV;
        outAncCreator = it->second.creatorProcess;
        return kPionZero;
      }
      if (apd == 211) {                     // π± ancestor
        outAncTrk     = cur;
        outAncPDG     = it->second.pdg;
        outAncKE0     = it->second.ke0_MeV;
        outAncCreator = it->second.creatorProcess;
        return (it->second.parentID == 0) ? kPrimaryPion
               : (it->second.pdg > 0 ? kSecPionPlus : kSecPionMinus);
      }
      if (apd == 13) {                      // muon ancestor
        outAncTrk     = cur;
        outAncPDG     = it->second.pdg;
        outAncKE0     = it->second.ke0_MeV;
        outAncCreator = it->second.creatorProcess;
        return kMuon;
      }
      cur = it->second.parentID;
    }
    return kGamma;   // no recognised ancestor
  }

  // ----------------------------------------------------------------
  // Nuclear fragments: proton, neutron, heavy ions
  // ----------------------------------------------------------------
  if (IsNuclearFragment(pdg)) {
    return kNuclear;
  }

  // ----------------------------------------------------------------
  // e± or other particles not caught above: walk ancestry chain
  // to find the interesting ancestor (same logic as ClassifyPhoton).
  // ----------------------------------------------------------------
  {
    int cur = parentID;
    while (cur > 0) {
      auto it = g_trackAncestry.find(cur);
      if (it == g_trackAncestry.end()) break;
      const AncestryInfo& info = it->second;
      const int apd = std::abs(info.pdg);

      if (apd == 211) {
        outAncTrk = cur; outAncPDG = info.pdg;
        outAncKE0 = info.ke0_MeV; outAncCreator = info.creatorProcess;
        if (info.parentID == 0) return kPrimaryPion;
        return (info.pdg > 0) ? kSecPionPlus : kSecPionMinus;
      }
      if (info.pdg == 111) {
        outAncTrk = cur; outAncPDG = 111;
        outAncKE0 = info.ke0_MeV; outAncCreator = info.creatorProcess;
        return kPionZero;
      }
      if (apd == 13) {
        outAncTrk = cur; outAncPDG = info.pdg;
        outAncKE0 = info.ke0_MeV; outAncCreator = info.creatorProcess;
        return kMuon;
      }
      if (IsNuclearFragment(info.pdg)) {
        outAncTrk = cur; outAncPDG = info.pdg;
        outAncKE0 = info.ke0_MeV; outAncCreator = info.creatorProcess;
        return kNuclear;
      }
      cur = info.parentID;
    }
  }

  return kOther;
}