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


double TrackLengthDebugUtils::getHelixArcLength(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2){
    double omega = ts1.getOmega();
    double z1 = ts1.getReferencePoint()[2] + ts1.getZ0();
    double z2 = ts2.getReferencePoint()[2] + ts2.getZ0();
    double dPhi = std::abs( ts2.getPhi() - ts1.getPhi() );
    if (dPhi > M_PI) dPhi = 2*M_PI - dPhi;

    return std::sqrt( std::pow(dPhi/omega, 2) + std::pow(z2-z1, 2) );
}


double TrackLengthDebugUtils::getHelixArcLengthOld(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2){
    double omega = ts1.getOmega();
    double tanL = ts1.getTanLambda();
    double dPhi = std::abs( ts2.getPhi() - ts1.getPhi() );
    if (dPhi > M_PI) dPhi = 2*M_PI - dPhi;

    return std::abs(dPhi/omega) * std::sqrt( 1. + std::pow(tanL, 2) );
}


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

EVENT::SimTrackerHit* TrackLengthDebugUtils::getSimTrackerHit(EVENT::LCEvent* evt, EVENT::TrackerHit* hit){
    // we will have to check all collections.. Oh my god why..
    const std::vector<std::string>* colNames = evt->getCollectionNames();
    bool hasVXDCollection = (std::find(colNames->begin(), colNames->end(), "VXDTrackerHitRelations") != colNames->end());
    bool hasTPCCollection = (std::find(colNames->begin(), colNames->end(), "TPCTrackerHitRelations") != colNames->end());
    bool hasSITCollection = (std::find(colNames->begin(), colNames->end(), "SITTrackerHitRelations") != colNames->end());
    bool hasSETCollection = (std::find(colNames->begin(), colNames->end(), "SETTrackerHitRelations") != colNames->end());
    bool hasFTDPixelCollection = (std::find(colNames->begin(), colNames->end(), "FTDPixelTrackerHitRelations") != colNames->end());
    bool hasFTDStripCollection = (std::find(colNames->begin(), colNames->end(), "FTDStripTrackerHitRelations") != colNames->end());
    if (hasTPCCollection){
        LCRelationNavigator nav( evt->getCollection("TPCTrackerHitRelations") );
        if (nav.getRelatedToObjects(hit).size() > 0){
            const vector<float>& weights = nav.getRelatedToWeights(hit);
            int max_i = std::max_element(weights.begin(), weights.end()) - weights.begin();
            SimTrackerHit* const simHit = static_cast<SimTrackerHit* const> (nav.getRelatedToObjects(hit)[max_i]);
            std::cout<<"Return TPC sim hit "<<simHit<<" to the reco hit "<<hit<<" with weight "<<weights[max_i]<<"    TYPE:    "<<hit->getType()<<std::endl;
            return simHit;
        }        
    }
    if (hasVXDCollection){
        LCRelationNavigator nav( evt->getCollection("VXDTrackerHitRelations") );
        if (nav.getRelatedToObjects(hit).size() > 0){
            const vector<float>& weights = nav.getRelatedToWeights(hit);
            int max_i = std::max_element(weights.begin(), weights.end()) - weights.begin();
            SimTrackerHit* simHit = static_cast<SimTrackerHit*> (nav.getRelatedToObjects(hit)[max_i]);
            std::cout<<"Return VXD sim hit "<<simHit<<" to the reco hit "<<hit<<" with weight "<<weights[max_i]<<"    TYPE:    "<<hit->getType()<<std::endl;
            return simHit;
        }
    }
    if (hasSITCollection){
        LCRelationNavigator nav( evt->getCollection("SITTrackerHitRelations") );
        if (nav.getRelatedToObjects(hit).size() > 0){
            const vector<float>& weights = nav.getRelatedToWeights(hit);
            int max_i = std::max_element(weights.begin(), weights.end()) - weights.begin();
            SimTrackerHit* simHit = static_cast<SimTrackerHit*> (nav.getRelatedToObjects(hit)[max_i]);
            std::cout<<"Return SIT sim hit "<<simHit<<" to the reco hit "<<hit<<" with weight "<<weights[max_i]<<"    TYPE:    "<<hit->getType()<<std::endl;
            return simHit;
        }
    }
    if (hasSETCollection){
        LCRelationNavigator nav( evt->getCollection("SETTrackerHitRelations") );
        if (nav.getRelatedToObjects(hit).size() > 0){
            const vector<float>& weights = nav.getRelatedToWeights(hit);
            int max_i = std::max_element(weights.begin(), weights.end()) - weights.begin();
            SimTrackerHit* simHit = static_cast<SimTrackerHit*> (nav.getRelatedToObjects(hit)[max_i]);
            std::cout<<"Return SET sim hit "<<simHit<<" to the reco hit "<<hit<<" with weight "<<weights[max_i]<<"    TYPE:    "<<hit->getType()<<std::endl;
            return simHit;
        }
    }
    if (hasFTDPixelCollection){
        LCRelationNavigator nav( evt->getCollection("FTDPixelTrackerHitRelations") );
        if (nav.getRelatedToObjects(hit).size() > 0){
            const vector<float>& weights = nav.getRelatedToWeights(hit);
            int max_i = std::max_element(weights.begin(), weights.end()) - weights.begin();
            SimTrackerHit* simHit = static_cast<SimTrackerHit*> (nav.getRelatedToObjects(hit)[max_i]);
            std::cout<<"Return FTD PIXEL sim hit "<<simHit<<" to the reco hit "<<hit<<" with weight "<<weights[max_i]<<"    TYPE:    "<<hit->getType()<<std::endl;
            return simHit;
        }
    }
    if (hasFTDStripCollection){
        LCRelationNavigator nav( evt->getCollection("FTDStripTrackerHitRelations") );
        if (nav.getRelatedToObjects(hit).size() > 0){
            const vector<float>& weights = nav.getRelatedToWeights(hit);
            int max_i = std::max_element(weights.begin(), weights.end()) - weights.begin();
            SimTrackerHit* simHit = static_cast<SimTrackerHit*> (nav.getRelatedToObjects(hit)[max_i]);
            std::cout<<"Return FTD STRIP sim hit "<<simHit<<" to the reco hit "<<hit<<" with weight "<<weights[max_i]<<"    TYPE:    "<<hit->getType()<<std::endl;
            return simHit;
        }
    }
    std::cout<<"Couldn't find MC sim hit in any of the collections for the hit "<<hit<<"    TYPE:    "<<hit->getType()<<std::endl;
    return nullptr;
}

double TrackLengthDebugUtils::getHelixNRevolutions(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2){
    double omega = ts1.getOmega();
    double tanL = ts1.getTanLambda();
    double z1 = ts1.getReferencePoint()[2] + ts1.getZ0();
    double z2 = ts2.getReferencePoint()[2] + ts2.getZ0();

    // helix length projected on xy
    double circHelix = std::abs( (z2-z1)/tanL );
    double circFull = 2*M_PI/std::abs(omega);

    return circHelix/circFull;
}

std::vector<EVENT::SimTrackerHit*> TrackLengthDebugUtils::convertHitsToSimHits(EVENT::LCEvent* evt, const std::vector<EVENT::Track*>& tracks, MarlinTrk::IMarlinTrkSystem* trkSystem, double bField){
    std::vector<SimTrackerHit*> simHits;
    int nTracks = tracks.size();
    for(int i=0; i<nTracks; ++i){
        Track* track = tracks[i];
        vector <TrackerHit*> hits = track->getTrackerHits();
        std::sort(hits.begin(), hits.end(), sortByRho);

        // setup initial dummy covariance matrix
        vector<float> covMatrix(15);
        // initialize variancesk
        covMatrix[0]  = 1e+06; //sigma_d0^2
        covMatrix[2]  = 100.; //sigma_phi0^2
        covMatrix[5]  = 0.00001; //sigma_omega^2
        covMatrix[9]  = 1e+06; //sigma_z0^2
        covMatrix[14] = 100.; //sigma_tanl^2
        double maxChi2PerHit = 100.;
        std::unique_ptr<IMarlinTrack> marlinTrk( trkSystem->createTrack() );
        TrackImpl refittedTrack;

        //Need to initialize trackState at last hit
        TrackStateImpl preFit = *track->getTrackState(TrackState::AtLastHit);
        preFit.setCovMatrix( covMatrix );
        int errorFit = MarlinTrk::createFinalisedLCIOTrack(marlinTrk.get(), hits, &refittedTrack, IMarlinTrack::backward, &preFit, bField, maxChi2PerHit);
        //if fit fails, try also fit forward
        if (errorFit != 0){
            streamlog_out(DEBUG8)<<"Fit backward fails! Trying to fit forward for "<<i+1<<" subTrack in this PFO!"<<std::endl;

            marlinTrk.reset( trkSystem->createTrack() );
            errorFit = MarlinTrk::createFinalisedLCIOTrack(marlinTrk.get(), hits, &refittedTrack, IMarlinTrack::forward, &preFit, bField, maxChi2PerHit);
        }
        if (errorFit != 0){
            streamlog_out(WARNING)<<"Fit fails in both directions. Skipping "<<i+1<<" subTrack in this PFO!"<<std::endl;
            continue;
        }
        //here hits are sorted by rho=(x^2+y^2) in the fit direction. forward - increasing rho, backward - decreasing rho
        vector< pair<TrackerHit*, double> > hitsInFit;
        marlinTrk->getHitsInFit(hitsInFit);

        //Find which way to loop over the array of hits. We need to loop in the direction of the track.
        bool loopForward = true;
        double zFirst = std::abs( hitsInFit.front().first->getPosition()[2] );
        double zLast = std::abs( hitsInFit.back().first->getPosition()[2] );

        // OPTIMIZE: 10 mm is just a round number. With very small z difference it is more robust to use rho, to be sure z difference is not caused by tpc Z resolution or multiple scattering
        if ( std::abs(zLast - zFirst) > 10. ){
            if ( zLast < zFirst ) loopForward = false;
            streamlog_out(DEBUG8)<<"Using Z to define loop direction over subTrack hits."<<std::endl;
            streamlog_out(DEBUG8)<<"subTrack "<<i+1<<" zFirst: "<<hitsInFit.front().first->getPosition()[2]<<" zLast: "<<hitsInFit.back().first->getPosition()[2]<<" loop forward: "<<loopForward<<std::endl;
        }
        else{
            double rhoFirst = std::hypot( hitsInFit.front().first->getPosition()[0], hitsInFit.front().first->getPosition()[1] );
            double rhoLast = std::hypot( hitsInFit.back().first->getPosition()[0], hitsInFit.back().first->getPosition()[1] );
            if ( rhoLast < rhoFirst ) loopForward = false;
            streamlog_out(DEBUG8)<<"Track is very perpendicular (dz < 10 mm). Using rho to define loop direction over subTrack hits."<<std::endl;
            streamlog_out(DEBUG8)<<"subTrack "<<i+1<<" zFirst: "<<hitsInFit.front().first->getPosition()[2]<<" zLast: "<<hitsInFit.back().first->getPosition()[2]<<std::endl;
            streamlog_out(DEBUG8)<<"subTrack "<<i+1<<" rhoFirst: "<<rhoFirst<<" rhoLast: "<<rhoLast<<" loop forward: "<<loopForward<<std::endl;
        }

        int nHitsInFit = hitsInFit.size();
        // NOTE: although we use z to understand track direction, hits are still sorted by rho
        if (loopForward){
            for( int j=0; j<nHitsInFit; ++j ){
                if( simHits.empty() ) simHits.push_back( getSimTrackerHit(evt, hitsInFit[j].first) );
                simHits.push_back( getSimTrackerHit(evt, hitsInFit[j].first) );
            }
        }
        else{
            for( int j=nHitsInFit-1; j>=0; --j ){
                if( simHits.empty() ) simHits.push_back( getSimTrackerHit(evt, hitsInFit[j].first) );
                simHits.push_back( getSimTrackerHit(evt, hitsInFit[j].first) );            }
        }

        // OPTIMIZE: if last subtrack fit fails in both directions we don't add track state at the ECal.
        if (i == nTracks - 1){
            // SET hit is not present in hitsInFit as it is composite hit from strips
            // Add ts at the SET hit manualy which fitter returns with reffited track
            // If LastHit != SET hit, then we duplicate previous track state at last TPC hit
            // isn't pretty, but shouldn't affect the track length
            simHits.push_back( simHits.back() );
            simHits.push_back( simHits.back() );
        }
    }
    // one can maybe use hits of refittedTrack, but they include also hits that had failed in the fit
    // code would look cleaner, but using hits that are failed in fit probably would have worse performance..
    // needs to be checked.
    return simHits;
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
