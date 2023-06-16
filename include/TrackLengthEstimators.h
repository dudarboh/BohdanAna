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

std::vector<EVENT::Track*> getSubTracks(EVENT::Track* track);

std::vector<IMPL::TrackStateImpl> getTrackStates(EVENT::ReconstructedParticle* pfo, double bField, MarlinTrk::IMarlinTrkSystem* trkSystem);

double getHelixNRevolutions(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2);


//Track length from the IDR production version (1 September 2020)
double getTrackLengthIDR(EVENT::Track* track);

// simple helix approxumation (SHA) option 1 from the thesis
double getTrackLengthSHA1(EVENT::Track* track);

// simple helix approxumation (SHA) option 2 from the thesis
double getTrackLengthSHA2(EVENT::Track* track);

// simple helix approxumation (SHA) option 3 from the thesis
double getTrackLengthSHA3(EVENT::Track* track);

// simple helix approxumation (SHA) option 4 from the thesis
double getTrackLengthSHA4(EVENT::Track* track);

// simple helix approxumation (SHA) option 5 from the thesis
double getTrackLengthSHA5(EVENT::Track* track);

// simple helix approxumation (SHA) option 6 from the thesis
double getTrackLengthSHA6(EVENT::Track* track);


// iterative Kalman Filter (IKF) option 1 from the thesis
double getTrackLengthIKF1(std::vector<IMPL::TrackStateImpl> trackStates);

// iterative Kalman Filter (IKF) option 2 from the thesis
double getTrackLengthIKF2(std::vector<IMPL::TrackStateImpl> trackStates);

// iterative Kalman Filter (IKF) option 3 from the thesis
double getTrackLengthIKF3(std::vector<IMPL::TrackStateImpl> trackStates);




#endif
