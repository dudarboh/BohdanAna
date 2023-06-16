#ifndef BohdanUtils_h
#define BohdanUtils_h 1

#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Track.h"
#include "EVENT/TrackState.h"
#include "IMPL/TrackStateImpl.h"
#include "MarlinTrk/IMarlinTrack.h"

#include "UTIL/PIDHandler.h"

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

void printTrackStateLong(IMPL::TrackStateImpl ts);

void printTrackStateShort(IMPL::TrackStateImpl ts);

#endif