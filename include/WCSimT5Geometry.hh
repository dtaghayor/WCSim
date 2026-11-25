#ifndef WCSimT5Geometry_h
#define WCSimT5Geometry_h

// WCTE T5 trigger hodoscope: bar layout shared by the geometry
// (WCSimConstructBeamPipe.cc) and the GPS beam-spot acceptance
// (WCSimPrimaryGeneratorAction::ApplyGPSBeamSpot), so the two cannot drift.
//
// Bars are stacked along local y (bar 0 at the bottom), centred on the beam
// pipe axis, each centred at x = 0 with its length along local x.

#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include <cmath>

namespace WCSimT5 {

  constexpr G4int    nBars     = 8;
  constexpr G4double barHeight = 16.25 * CLHEP::mm;  // Y dimension
  constexpr G4double barLengths[nBars] = {           // X dimension
    41 * CLHEP::mm, 94 * CLHEP::mm, 112 * CLHEP::mm, 123 * CLHEP::mm,
    123 * CLHEP::mm, 112 * CLHEP::mm, 94 * CLHEP::mm, 41 * CLHEP::mm};

  // Index of the bar covering transverse point (x, y) relative to the pipe
  // axis, or -1 if the point misses every bar.
  inline G4int BarIndex(G4double x, G4double y)
  {
    const G4double fromBottom = y + 0.5 * nBars * barHeight;
    if (fromBottom < 0.) return -1;
    const G4int i = static_cast<G4int>(std::floor(fromBottom / barHeight));
    if (i >= nBars) return -1;
    return (std::abs(x) <= 0.5 * barLengths[i]) ? i : -1;
  }

  inline G4bool Contains(G4double x, G4double y) { return BarIndex(x, y) >= 0; }

} // namespace WCSimT5

#endif
