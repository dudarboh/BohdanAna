#ifndef TrackLengthEstimators_h
#define TrackLengthEstimators_h 1

#include "EVENT/ReconstructedParticle.h"
#include "IMPL/TrackStateImpl.h"
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "DDRec/Vector3D.h"
#include "EVENT/Track.h"
#include <vector>
#include <memory>
/**
 * Different estimators to calculate track lengths that are used by the TrackLengthProcessor.
 *
 * \author B. Dudar, DESY, 2022
*/
namespace TrackLengthDebugUtils{




//Track length from the IDR production version (1 September 2020)
double getTrackLengthIDR(EVENT::Track* track);

//Track length from the IDR production version (5 October 2020)
// but make it always positive with std::abs
double getTrackLengthIDR2(EVENT::Track* track);

//Track length from the IDR production version (5 October 2020)
// but make it always positive with std::abs
// and use track state at calorimeter for omega, tanL
double getTrackLengthIDR3(EVENT::Track* track);

//Track length from the IDR production version (8 September 2021)
// but make it always positive with std::abs
// and use track state at calorimeter for omega, tanL
// and fix phi flip bug
double getTrackLengthIDR4(EVENT::Track* track);

//Track length using Winni's formula for helix (8 September 2021).
double getTrackLengthWinni(EVENT::Track* track, double bField, MarlinTrk::IMarlinTrkSystem* trkSystem);

//Track length using Winni's formula for helix (18 March 2022).
// but improve hit ordering (use z, not rho) and try to fit forward in case backward fit fails
double getTrackLengthWinni2(EVENT::Track* track, double bField, MarlinTrk::IMarlinTrkSystem* trkSystem);

//Track length using helix formula w/o omega, only dz and tanL (30 January 2023).
double getTrackLengthUsingZ(EVENT::Track* track, double bField, MarlinTrk::IMarlinTrkSystem* trkSystem);

//Track length using helix formula w/o omega, only dz and tanL (30 January 2023).
// but using only momentum information from SimHits ("true level")
double getTrackLengthSimUsingZ(EVENT::Track* track, double bField, MarlinTrk::IMarlinTrkSystem* trkSystem);


}



#endif
