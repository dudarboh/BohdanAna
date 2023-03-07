#include "TrackLengthEstimators.h"
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


double TrackLengthDebugUtils::getTrackLengthIDR(EVENT::Track* track){
    const TrackState* tsIP = track->getTrackState( TrackState::AtIP );
    const TrackState* tsCalo = track->getTrackState( TrackState::AtCalorimeter );
    
    float phiIP = tsIP->getPhi();
    float phiCalo = tsCalo->getPhi();
    float omega = tsIP->getOmega();
    float tanL = tsIP->getTanLambda();    
    return (phiIP-phiCalo)*(1/omega)*sqrt(1+tanL*tanL);
}

double TrackLengthDebugUtils::getTrackLengthIDR2(EVENT::Track* track){
    return std::abs( getTrackLengthIDR(track) );
}

double TrackLengthDebugUtils::getTrackLengthIDR3(EVENT::Track* track){
    const TrackState* tsIP = track->getTrackState( TrackState::AtIP );
    const TrackState* tsCalo = track->getTrackState( TrackState::AtCalorimeter );

    float omega = tsCalo->getOmega();
    float tanL = tsCalo->getTanLambda();        
    float phiIP = tsIP->getPhi();
    float phiCalo = tsCalo->getPhi();
    return std::abs( (phiIP-phiCalo)/omega )*std::sqrt(1+tanL*tanL);
}

double TrackLengthDebugUtils::getTrackLengthIDR4(EVENT::Track* track){
    const TrackState* tsIP = track->getTrackState( TrackState::AtIP );
    const TrackState* tsCalo = track->getTrackState( TrackState::AtCalorimeter );
    
    float omega = std::abs( tsCalo->getOmega() );
    float tanL = std::abs( tsCalo->getTanLambda() );
    double dPhi = std::abs( tsCalo->getPhi() - tsIP->getPhi() );
    if (dPhi > M_PI) dPhi = 2*M_PI - dPhi;
    return dPhi/omega*std::sqrt(1+tanL*tanL);
}


double TrackLengthDebugUtils::getTrackLengthWinni(EVENT::Track* track, double bField, MarlinTrk::IMarlinTrkSystem* trkSystem){
    std::vector<Track*> subTracks = getSubTracks(track);

    // extract track state per every hit
    std::vector<TrackStateImpl> trackStates;
    int nTracks = subTracks.size();
    int nGoodFits = 0;  
    for(int i=0; i<nTracks; ++i){
        std::vector <TrackerHit*> hits = subTracks[i]->getTrackerHits();
        std::sort(hits.begin(), hits.end(), sortByRho);

        // setup initial dummy covariance matrix
        vector<float> covMatrix(15);
        // initialize variances
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
        if (errorFit == 0) nGoodFits++;
        if (errorFit != 0) continue;

        vector< pair<TrackerHit*, double> > hitsInFit;
        marlinTrk->getHitsInFit(hitsInFit);

        int nHitsInFit = hitsInFit.size();
        // if first subTrack
        if (nGoodFits == 1){
            trackStates.push_back(*(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtIP)) ));

            //add hits in increasing rho for the FIRST subTrack!!!!!
            for( int j=nHitsInFit-1; j>=0; --j ){
                TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                trackStates.push_back(ts);
            }
        }
        else{
            // check which hit is closer to the last hit of previous fit.
            // and iterate starting from the closest
            Vector3D innerHit ( hitsInFit.back().first->getPosition() );
            Vector3D outerHit ( hitsInFit.front().first->getPosition() );
            Vector3D prevHit ( trackStates.back().getReferencePoint() );

            if ( (innerHit - prevHit).r() < (outerHit - prevHit).r() ){
                for( int j=nHitsInFit-1; j>=0; --j ){
                    //iterate in increasing rho
                    TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                    trackStates.push_back(ts);
                }
            }
            else{
                for( int j=0; j<nHitsInFit; ++j ){
                    //iterate in decreasing rho
                    TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                    trackStates.push_back(ts);
                }
            }
        }
        //if last subTrack
        if (i == nTracks - 1){
            // SET hit is not present in hitsInFit as it is composite hit from strips
            // Add ts at the SET hit manualy which fitter returns with reffited track
            // If LastHit != SET hit, then we duplicate previous track state at last TPC hit
            // isn't pretty, but shouldn't affect the track length
            trackStates.push_back( *(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtLastHit)) ) );
            trackStates.push_back( *(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtCalorimeter) ) ) );
        }
    }

    double trackLength = 0.;
    int nTrackStates = trackStates.size();
    if (nTrackStates <= 1) return trackLength;
    for( int j=1; j < nTrackStates; ++j ){
        //we check which track length formula to use
        double nTurns = getHelixNRevolutions( trackStates[j-1], trackStates[j] );
        double arcLength;
        // we cannot calculate arc length for more than pi revolution using delta phi. Use formula with only z
        if ( nTurns <= 0.5 ) arcLength = getHelixArcLength( trackStates[j-1], trackStates[j] );
        else arcLength = getHelixLengthAlongZ( trackStates[j-1], trackStates[j] );
        trackLength += arcLength;
    }
    return trackLength;
}


double TrackLengthDebugUtils::getTrackLengthWinni2(EVENT::Track* track, double bField, MarlinTrk::IMarlinTrkSystem* trkSystem){
    std::vector<Track*> subTracks = getSubTracks(track);

    // extract track state per every hit
    std::vector<TrackStateImpl> trackStates;
    int nTracks = subTracks.size();
    for(int i=0; i<nTracks; ++i){
        std::vector <TrackerHit*> hits = subTracks[i]->getTrackerHits();
        std::sort(hits.begin(), hits.end(), sortByRho);

        // setup initial dummy covariance matrix
        vector<float> covMatrix(15);
        // initialize variances
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
            marlinTrk.reset( trkSystem->createTrack() );
            errorFit = MarlinTrk::createFinalisedLCIOTrack(marlinTrk.get(), hits, &refittedTrack, IMarlinTrack::forward, &preFit, bField, maxChi2PerHit);
        }
        if (errorFit != 0) continue;

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
        }
        else{
            double rhoFirst = std::hypot( hitsInFit.front().first->getPosition()[0], hitsInFit.front().first->getPosition()[1] );
            double rhoLast = std::hypot( hitsInFit.back().first->getPosition()[0], hitsInFit.back().first->getPosition()[1] );
            if ( rhoLast < rhoFirst ) loopForward = false;
        }

        int nHitsInFit = hitsInFit.size();
        // if first subTrack
        if ( trackStates.empty() ) trackStates.push_back(*(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtIP)) ));
        if (loopForward){
            //add hits in increasing rho for the FIRST subTrack!!!!!
            for( int j=0; j<nHitsInFit; ++j ){
                TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                trackStates.push_back(ts);
            }
        }
        else{
            for( int j=nHitsInFit-1; j>=0; --j ){
                TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                trackStates.push_back(ts);
            }
        }
        if (i == nTracks - 1){
            // SET hit is not present in hitsInFit as it is composite hit from strips
            // Add ts at the SET hit manualy which fitter returns with reffited track
            // If LastHit != SET hit, then we duplicate previous track state at last TPC hit
            // isn't pretty, but shouldn't affect the track length
            trackStates.push_back( *(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtLastHit)) ) );
            trackStates.push_back( *(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtCalorimeter) ) ) );
        }
    }

    double trackLength = 0.;
    int nTrackStates = trackStates.size();
    if (nTrackStates <= 1) return trackLength;
    for( int j=1; j < nTrackStates; ++j ){
        //we check which track length formula to use
        double nTurns = getHelixNRevolutions( trackStates[j-1], trackStates[j] );
        double arcLength;
        // we cannot calculate arc length for more than pi revolution using delta phi. Use formula with only z
        if ( nTurns <= 0.5 ) arcLength = getHelixArcLength( trackStates[j-1], trackStates[j] );
        else arcLength = getHelixLengthAlongZ( trackStates[j-1], trackStates[j] );
        trackLength += arcLength;
    }
    return trackLength;
}


double TrackLengthDebugUtils::getTrackLengthUsingZ(EVENT::Track* track, double bField, MarlinTrk::IMarlinTrkSystem* trkSystem){
    std::vector<Track*> subTracks = getSubTracks(track);

    // extract track state per every hit
    std::vector<TrackStateImpl> trackStates;
    int nTracks = subTracks.size();
    for(int i=0; i<nTracks; ++i){
        std::vector <TrackerHit*> hits = subTracks[i]->getTrackerHits();
        std::sort(hits.begin(), hits.end(), sortByRho);

        // setup initial dummy covariance matrix
        vector<float> covMatrix(15);
        // initialize variances
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
            marlinTrk.reset( trkSystem->createTrack() );
            errorFit = MarlinTrk::createFinalisedLCIOTrack(marlinTrk.get(), hits, &refittedTrack, IMarlinTrack::forward, &preFit, bField, maxChi2PerHit);
        }
        if (errorFit != 0) continue;

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
        }
        else{
            double rhoFirst = std::hypot( hitsInFit.front().first->getPosition()[0], hitsInFit.front().first->getPosition()[1] );
            double rhoLast = std::hypot( hitsInFit.back().first->getPosition()[0], hitsInFit.back().first->getPosition()[1] );
            if ( rhoLast < rhoFirst ) loopForward = false;
        }

        int nHitsInFit = hitsInFit.size();
        // if first subTrack
        if ( trackStates.empty() ) trackStates.push_back(*(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtIP)) ));
        if (loopForward){
            //add hits in increasing rho for the FIRST subTrack!!!!!
            for( int j=0; j<nHitsInFit; ++j ){
                TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                trackStates.push_back(ts);
            }
        }
        else{
            for( int j=nHitsInFit-1; j>=0; --j ){
                TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                trackStates.push_back(ts);
            }
        }
        if (i == nTracks - 1){
            // SET hit is not present in hitsInFit as it is composite hit from strips
            // Add ts at the SET hit manualy which fitter returns with reffited track
            // If LastHit != SET hit, then we duplicate previous track state at last TPC hit
            // isn't pretty, but shouldn't affect the track length
            trackStates.push_back( *(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtLastHit)) ) );
            trackStates.push_back( *(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtCalorimeter) ) ) );
        }
    }

    double trackLength = 0.;
    int nTrackStates = trackStates.size();
    if (nTrackStates <= 1) return trackLength;
    for( int j=1; j < nTrackStates; ++j ) trackLength += getHelixLengthAlongZ( trackStates[j-1], trackStates[j] );
    return trackLength;
}



double TrackLengthDebugUtils::getTrackLengthUsingZ2(EVENT::ReconstructedParticle* pfo, double bField, MarlinTrk::IMarlinTrkSystem* trkSystem){
    if ( pfo->getTracks().empty() ) return 0.;
    Track* track = pfo->getTracks()[0];

    std::vector<Track*> subTracks = getSubTracks(track);

    // extract track state per every hit
    std::vector<TrackStateImpl> trackStates;
    int nTracks = subTracks.size();
    TrackImpl lastGoodRefittedTrack;
    for(int i=0; i<nTracks; ++i){
        std::vector <TrackerHit*> hits = subTracks[i]->getTrackerHits();
        std::sort(hits.begin(), hits.end(), sortByRho);

        // setup initial dummy covariance matrix
        vector<float> covMatrix(15);
        // initialize variances
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
            marlinTrk.reset( trkSystem->createTrack() );
            errorFit = MarlinTrk::createFinalisedLCIOTrack(marlinTrk.get(), hits, &refittedTrack, IMarlinTrack::forward, &preFit, bField, maxChi2PerHit);
        }
        if (errorFit != 0) continue;
        lastGoodRefittedTrack = refittedTrack;

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
        }
        else{
            double rhoFirst = std::hypot( hitsInFit.front().first->getPosition()[0], hitsInFit.front().first->getPosition()[1] );
            double rhoLast = std::hypot( hitsInFit.back().first->getPosition()[0], hitsInFit.back().first->getPosition()[1] );
            if ( rhoLast < rhoFirst ) loopForward = false;
        }

        int nHitsInFit = hitsInFit.size();
        // if first subTrack
        if ( trackStates.empty() ) trackStates.push_back(*(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtIP)) ));
        if (loopForward){
            //add hits in increasing rho for the FIRST subTrack!!!!!
            for( int j=0; j<nHitsInFit; ++j ){
                TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                trackStates.push_back(ts);
            }
        }
        else{
            for( int j=nHitsInFit-1; j>=0; --j ){
                TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                trackStates.push_back(ts);
            }
        }
        if (i == nTracks - 1){
            // SET hit is not present in hitsInFit as it is composite hit from strips
            // Add ts at the SET hit manualy which fitter returns with reffited track
            // If LastHit != SET hit, then we duplicate previous track state at last TPC hit
            // isn't pretty, but shouldn't affect the track length
            trackStates.push_back( *(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtLastHit)) ) );
        }
    }

    const TrackStateImpl* tsCalo = static_cast<const TrackStateImpl*> (lastGoodRefittedTrack.getTrackState(TrackState::AtCalorimeter) );
    if ( pfo->getClusters().size() > 0 && tsCalo != nullptr ) trackStates.push_back( *(tsCalo) );

    double trackLength = 0.;
    int nTrackStates = trackStates.size();
    if (nTrackStates <= 1) return trackLength;
    for( int j=1; j < nTrackStates; ++j ) trackLength += getHelixLengthAlongZ( trackStates[j-1], trackStates[j] );
    return trackLength;
}



double TrackLengthDebugUtils::getTrackLengthUsingZ3(EVENT::ReconstructedParticle* pfo, double bField, MarlinTrk::IMarlinTrkSystem* trkSystem){
    if ( pfo->getTracks().empty() ) return 0.;
    Track* track = pfo->getTracks()[0];

    std::vector<Track*> subTracks = getSubTracks(track);

    // extract track state per every hit
    std::vector<TrackStateImpl> trackStates;
    int nTracks = subTracks.size();
    TrackImpl lastGoodRefittedTrack;
    for(int i=0; i<nTracks; ++i){
        std::vector <TrackerHit*> hits = subTracks[i]->getTrackerHits();
        std::sort(hits.begin(), hits.end(), sortByRho);

        // setup initial dummy covariance matrix
        vector<float> covMatrix(15);
        // initialize variances
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
            marlinTrk.reset( trkSystem->createTrack() );
            errorFit = MarlinTrk::createFinalisedLCIOTrack(marlinTrk.get(), hits, &refittedTrack, IMarlinTrack::forward, &preFit, bField, maxChi2PerHit);
        }
        if (errorFit != 0) continue;
        lastGoodRefittedTrack = refittedTrack;

        //here hits are sorted by rho=(x^2+y^2) in the fit direction. forward - increasing rho, backward - decreasing rho
        vector< pair<TrackerHit*, double> > hitsInFit;
        marlinTrk->getHitsInFit(hitsInFit);

        //Find which way to loop over the array of hits. We need to loop in the direction of the track.
        double zFirst = std::abs( hitsInFit.front().first->getPosition()[2] );
        double zLast = std::abs( hitsInFit.back().first->getPosition()[2] );
        bool loopForward = ( zLast >= zFirst );
        int nHitsInFit = hitsInFit.size();
        // if first subTrack
        if ( trackStates.empty() ) trackStates.push_back(*(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtIP)) ));
        if (loopForward){
            //add hits in increasing rho for the FIRST subTrack!!!!!
            for( int j=0; j<nHitsInFit; ++j ){
                TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                trackStates.push_back(ts);
            }
        }
        else{
            for( int j=nHitsInFit-1; j>=0; --j ){
                TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                trackStates.push_back(ts);
            }
        }
        if (i == nTracks - 1){
            // SET hit is not present in hitsInFit as it is composite hit from strips
            // Add ts at the SET hit manualy which fitter returns with reffited track
            // If LastHit != SET hit, then we duplicate previous track state at last TPC hit
            // isn't pretty, but shouldn't affect the track length
            trackStates.push_back( *(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtLastHit)) ) );
        }
    }

    const TrackStateImpl* tsCalo = static_cast<const TrackStateImpl*> (lastGoodRefittedTrack.getTrackState(TrackState::AtCalorimeter) );
    if ( pfo->getClusters().size() > 0 && tsCalo != nullptr ) trackStates.push_back( *(tsCalo) );

    double trackLength = 0.;
    int nTrackStates = trackStates.size();
    if (nTrackStates <= 1) return trackLength;
    for( int j=1; j < nTrackStates; ++j ) trackLength += getHelixLengthAlongZ( trackStates[j-1], trackStates[j] );
    return trackLength;
}





//FIXME
double TrackLengthDebugUtils::getTrackLengthSimUsingZ(EVENT::Track* track, double bField, MarlinTrk::IMarlinTrkSystem* trkSystem){
    std::vector<Track*> subTracks = getSubTracks(track);

    // extract track state per every hit
    std::vector<TrackStateImpl> trackStates;
    int nTracks = subTracks.size();
    for(int i=0; i<nTracks; ++i){
        std::vector <TrackerHit*> hits = subTracks[i]->getTrackerHits();
        std::sort(hits.begin(), hits.end(), sortByRho);

        // setup initial dummy covariance matrix
        vector<float> covMatrix(15);
        // initialize variances
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
            marlinTrk.reset( trkSystem->createTrack() );
            errorFit = MarlinTrk::createFinalisedLCIOTrack(marlinTrk.get(), hits, &refittedTrack, IMarlinTrack::forward, &preFit, bField, maxChi2PerHit);
        }
        if (errorFit != 0) continue;

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
        }
        else{
            double rhoFirst = std::hypot( hitsInFit.front().first->getPosition()[0], hitsInFit.front().first->getPosition()[1] );
            double rhoLast = std::hypot( hitsInFit.back().first->getPosition()[0], hitsInFit.back().first->getPosition()[1] );
            if ( rhoLast < rhoFirst ) loopForward = false;
        }

        int nHitsInFit = hitsInFit.size();
        // if first subTrack
        if ( trackStates.empty() ) trackStates.push_back(*(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtIP)) ));
        if (loopForward){
            //add hits in increasing rho for the FIRST subTrack!!!!!
            for( int j=0; j<nHitsInFit; ++j ){
                TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                trackStates.push_back(ts);
            }
        }
        else{
            for( int j=nHitsInFit-1; j>=0; --j ){
                TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                trackStates.push_back(ts);
            }
        }
        if (i == nTracks - 1){
            // SET hit is not present in hitsInFit as it is composite hit from strips
            // Add ts at the SET hit manualy which fitter returns with reffited track
            // If LastHit != SET hit, then we duplicate previous track state at last TPC hit
            // isn't pretty, but shouldn't affect the track length
            trackStates.push_back( *(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtLastHit)) ) );
            trackStates.push_back( *(static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtCalorimeter) ) ) );
        }
    }

    double trackLength = 0.;
    int nTrackStates = trackStates.size();
    if (nTrackStates <= 1) return trackLength;
    for( int j=1; j < nTrackStates; ++j ) trackLength += getHelixLengthAlongZ( trackStates[j-1], trackStates[j] );
    return trackLength;
}
