#ifndef BohdanUtils_h
#define BohdanUtils_h 1

#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Track.h"
#include "EVENT/TrackState.h"
#include "EVENT/MCParticle.h"
#include "IMPL/TrackStateImpl.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "UTIL/LCRelationNavigator.h"
#include "UTIL/PIDHandler.h"
#include "UTIL/ILDConf.h"
#include "UTIL/Operators.h" // for debuging <<
#include "DDRec/Vector3D.h"
#include "EVENT/SimTrackerHit.h"
#include <string>
#include <vector>

int parseLine(char* line);

int getVirtualMemoryUsage();

int getPhysicalMemoryUsage();

std::vector<EVENT::Track*> getSubTracks(EVENT::Track* track);

float getParameterFromPID(EVENT::ReconstructedParticle* pfo, UTIL::PIDHandler& pidHandler, std::string algorithmName, std::string parameterName);

EVENT::TrackerHit* getSETHit(EVENT::Track* track);

const EVENT::TrackState* getTrackStateAtCalorimeter(EVENT::Track* track);

IMPL::TrackStateImpl getTrackStateAtHit(MarlinTrk::IMarlinTrack* marlinTrack, EVENT::TrackerHit* hit);

EVENT::MCParticle* getMC(EVENT::ReconstructedParticle* pfo, const UTIL::LCRelationNavigator& pfo2mc);

float getECALBarelRMin();

float getECALEndcapZMin();

dd4hep::rec::Vector3D getPhotonAtCalorimeter(EVENT::MCParticle* mc);

EVENT::SimTrackerHit* getSimTrackerHit(EVENT::TrackerHit* hit, const UTIL::LCRelationNavigator& navToSimTrackerHits);

unsigned long interpolateHexColor(unsigned long startColor, unsigned long endColor, float ratio);

bool isSETHit(const EVENT::TrackerHit* hit);

#endif