#ifndef TrackLengthDebugUtils_h
#define TrackLengthDebugUtils_h 1

#include "EVENT/ReconstructedParticle.h"
#include "IMPL/TrackStateImpl.h"
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "DDRec/Vector3D.h"
#include "EVENT/SimTrackerHit.h"
#include "IMPL/SimTrackerHitImpl.h"
#include "UTIL/PIDHandler.h"

#include <vector>
#include <memory>

namespace TrackLengthDebugUtils{
    bool sortByRho(EVENT::TrackerHit* a, EVENT::TrackerHit* b);
    IMPL::TrackStateImpl getTrackStateAtHit(MarlinTrk::IMarlinTrack* marlinTrk, EVENT::TrackerHit* hit);

    double getHelixLengthAlongZ(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2);

    std::vector<EVENT::Track*> getSubTracks(EVENT::Track* track);
    float getParameterFromPID(EVENT::ReconstructedParticle* pfo, UTIL::PIDHandler& pidHandler, std::string algorithmName, std::string parameterName);

}



#endif
