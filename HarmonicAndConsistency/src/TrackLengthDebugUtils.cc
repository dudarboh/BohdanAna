#include "TrackLengthDebugUtils.h"
#include "marlin/VerbosityLevels.h"
#include "marlinutil/CalorimeterHitType.h"
#include "HelixClass.h"
#include "CLHEP/Random/Randomize.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "MarlinTrk/MarlinTrkUtils.h"
#include "UTIL/ILDConf.h"
#include "UTIL/LCRelationNavigator.h"
#include "IMPL/TrackImpl.h"
#include <cmath>
#include <algorithm>
#include <limits>

using namespace EVENT;
using namespace UTIL;
using namespace IMPL;
using std::vector;
using std::numeric_limits;
using std::pair;
using dd4hep::rec::Vector3D;
using MarlinTrk::IMarlinTrack;
using MarlinTrk::IMarlinTrkSystem;
using CLHEP::RandGauss;


double TrackLengthDebugUtils::getHelixLengthAlongZ(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2){
    double tanL = ts1.getTanLambda();
    double z1 = ts1.getReferencePoint()[2] + ts1.getZ0();
    double z2 = ts2.getReferencePoint()[2] + ts2.getZ0();

    return std::abs( (z2-z1)/tanL ) * std::sqrt( 1.+tanL*tanL );
}


bool TrackLengthDebugUtils::sortByRho(EVENT::TrackerHit* a, EVENT::TrackerHit* b){
    Vector3D posA( a->getPosition() );
    Vector3D posB( b->getPosition() );
    return posA.rho() < posB.rho();
}


IMPL::TrackStateImpl TrackLengthDebugUtils::getTrackStateAtHit(MarlinTrk::IMarlinTrack* marlinTrk, EVENT::TrackerHit* hit){
    TrackStateImpl ts;
    double chi2Dummy;
    int ndfDummy;
    marlinTrk->getTrackState(hit, ts, chi2Dummy, ndfDummy);
    return ts;
}


std::vector<EVENT::Track*> TrackLengthDebugUtils::getSubTracks(EVENT::Track* track){
    vector<Track*> subTracks;
    // track itself contains VXD+FTD+SIT+TPC hits of the first curl.
    subTracks.push_back(track);

    int nSubTracks = track->getTracks().size();
    if (nSubTracks <= 1) return subTracks;

    int nTPCHits = track->getSubdetectorHitNumbers()[(ILDDetID::TPC)*2-1];
    int nSubTrack0Hits = track->getTracks()[0]->getTrackerHits().size();
    int nSubTrack1Hits = track->getTracks()[1]->getTrackerHits().size();

    //OPTIMIZE: this is not reliable, but I don't see any other way at the moment.
    //Read documentation in the header file for details.
    int startIdx;
    if( std::abs(nTPCHits - nSubTrack0Hits) <= 1  ) startIdx = 1;
    else if ( std::abs(nTPCHits - nSubTrack1Hits) <= 1 ) startIdx = 2;
    else{
        //FIXME: This happens very rarily (0.01%) for unknown reasons, so we just, skip adding subTracks...
        streamlog_out(WARNING)<<"Can't understand which subTrack is responsible for the first TPC hits! Skip adding subTracks."<<std::endl;
        return subTracks;
    }
    for(int j=startIdx; j < nSubTracks; ++j) subTracks.push_back( track->getTracks()[j] );
    return subTracks;
}

float TrackLengthDebugUtils::getParameterFromPID(EVENT::ReconstructedParticle* pfo, UTIL::PIDHandler& pidHandler, std::string algorithmName, std::string parameterName){
    int algorithmID = pidHandler.getAlgorithmID(algorithmName);
    const ParticleID& pfoPID = pidHandler.getParticleID(pfo, algorithmID);
    const std::vector<float>& parameters = pfoPID.getParameters();
    int parIdx = pidHandler.getParameterIndex(algorithmID, parameterName);
    return parameters[parIdx]; 
}
