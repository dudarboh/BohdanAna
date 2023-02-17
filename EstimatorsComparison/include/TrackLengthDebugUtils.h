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

    // MarlinTrk::IMarlinTrkSystem* _trkSystem = MarlinTrk::Factory::createMarlinTrkSystem("DDKalTest", nullptr, "");
    // _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useQMS, true);
    // _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::usedEdx, true);
    // _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing, true);
    // _trkSystem->init();


    bool sortByRho(EVENT::TrackerHit* a, EVENT::TrackerHit* b);
    IMPL::TrackStateImpl getTrackStateAtHit(MarlinTrk::IMarlinTrack* marlinTrk, EVENT::TrackerHit* hit);

    double getHelixArcLength(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2);

    double getHelixArcLengthOld(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2);

    double getHelixLengthAlongZ(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2);

    double getHelixNRevolutions(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2);


    std::vector<EVENT::Track*> getSubTracks(EVENT::Track* track);
    EVENT::SimTrackerHit* getSimTrackerHit(EVENT::LCEvent* evt, EVENT::TrackerHit* hit);
    std::vector<EVENT::SimTrackerHit*> convertHitsToSimHits(EVENT::LCEvent* evt, const std::vector<EVENT::Track*>& tracks, MarlinTrk::IMarlinTrkSystem* trkSystem, double bField);
    float getParameterFromPID(EVENT::ReconstructedParticle* pfo, UTIL::PIDHandler& pidHandler, std::string algorithmName, std::string parameterName);

}



#endif
