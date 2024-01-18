#include "TrackLength.h"
#include "BohdanUtils.h"

#include "marlin/VerbosityLevels.h"
#include "MarlinTrk/IMarlinTrack.h"

#include "MarlinTrk/MarlinTrkUtils.h"
#include "UTIL/TrackTools.h"
#include "IMPL/TrackImpl.h"

using namespace EVENT;
using dd4hep::rec::Vector3D;

std::vector<HitState> getTrackStates(ReconstructedParticle* pfo, double bField, MarlinTrk::IMarlinTrkSystem* trkSystem, const UTIL::LCRelationNavigator& navToSimTrackerHits){
    // Refit the track and extract track state at every tracker hit along the track
    std::vector<HitState> trackStates;
    if ( pfo->getTracks().empty() ) return trackStates;
    std::vector<Track*> subTracks = getSubTracks( pfo->getTracks()[0] );

    IMPL::TrackImpl lastGoodRefittedTrack;

    auto sortByRho = [](TrackerHit* a, TrackerHit* b) -> bool {
        Vector3D posA( a->getPosition() ), posB( b->getPosition() );
        return posA.rho() < posB.rho();
    };

    streamlog_out(DEBUG7)<<"PFOs track has "<<subTracks.size()<<" subTracks."<<std::endl;
    for(size_t i=0; i<subTracks.size(); ++i){
        std::vector <TrackerHit*> hits = subTracks[i]->getTrackerHits();

        streamlog_out(DEBUG7)<<"Subtrack "<<i+1<<" has "<<hits.size()<<" hits."<<std::endl;
        std::sort(hits.begin(), hits.end(), sortByRho);

        // setup initial dummy covariance matrix
        std::vector<float> covMatrix(15);
        // initialize variances
        covMatrix[0]  = 1e+06; //sigma_d0^2
        covMatrix[2]  = 100.; //sigma_phi0^2
        covMatrix[5]  = 0.00001; //sigma_omega^2
        covMatrix[9]  = 1e+06; //sigma_z0^2
        covMatrix[14] = 100.; //sigma_tanl^2
        double maxChi2PerHit = 100.;
        std::unique_ptr<MarlinTrk::IMarlinTrack> marlinTrk( trkSystem->createTrack() );
        IMPL::TrackImpl refittedTrack;

        //Need to initialize trackState at last hit
        IMPL::TrackStateImpl preFit = *(subTracks[i]->getTrackState(TrackState::AtLastHit));
        preFit.setCovMatrix( covMatrix );
        int errorFit = MarlinTrk::createFinalisedLCIOTrack(marlinTrk.get(), hits, &refittedTrack, MarlinTrk::IMarlinTrack::backward, &preFit, bField, maxChi2PerHit);
        //if fit fails, try also fit forward
        if (errorFit != 0){
            streamlog_out(DEBUG7)<<"Fit backward fails! Trying to fit forward for "<<i+1<<" subTrack in this PFO!"<<std::endl;
            marlinTrk.reset( trkSystem->createTrack() );
            errorFit = MarlinTrk::createFinalisedLCIOTrack(marlinTrk.get(), hits, &refittedTrack, MarlinTrk::IMarlinTrack::forward, &preFit, bField, maxChi2PerHit);
        }
        if (errorFit != 0){
            streamlog_out(WARNING)<<"Fit fails in both directions. Skipping "<<i+1<<" subTrack in this PFO!"<<std::endl;
            continue;
        }
        lastGoodRefittedTrack = refittedTrack;

        //here hits are sorted by rho=(x^2+y^2) in the fit direction. forward - increasing rho, backward - decreasing rho
        std::vector< std::pair<TrackerHit*, double> > hitsInFit;
        marlinTrk->getHitsInFit(hitsInFit);

        //Find which way to loop over the array of hits. We need to loop in the direction of the track.
        bool loopForward = true;
        double zFirst = std::abs( hitsInFit.front().first->getPosition()[2] );
        double zLast = std::abs( hitsInFit.back().first->getPosition()[2] );

        // OPTIMIZE: 10 mm is just a round number. With very small z difference it is more robust to use rho, to be sure z difference is not caused by tpc Z resolution or multiple scattering
        if ( std::abs(zLast - zFirst) > 10. ){
            if ( zLast < zFirst ) loopForward = false;
            streamlog_out(DEBUG7)<<"Using Z to define loop direction over subTrack hits."<<std::endl;
            streamlog_out(DEBUG7)<<"subTrack "<<i+1<<" zFirst: "<<hitsInFit.front().first->getPosition()[2]<<" zLast: "<<hitsInFit.back().first->getPosition()[2]<<" loop forward: "<<loopForward<<std::endl;
        }
        else{
            double rhoFirst = std::hypot( hitsInFit.front().first->getPosition()[0], hitsInFit.front().first->getPosition()[1] );
            double rhoLast = std::hypot( hitsInFit.back().first->getPosition()[0], hitsInFit.back().first->getPosition()[1] );
            if ( rhoLast < rhoFirst ) loopForward = false;
            streamlog_out(DEBUG7)<<"Track is very perpendicular (dz < 10 mm). Using rho to define loop direction over subTrack hits."<<std::endl;
            streamlog_out(DEBUG7)<<"subTrack "<<i+1<<" zFirst: "<<hitsInFit.front().first->getPosition()[2]<<" zLast: "<<hitsInFit.back().first->getPosition()[2]<<std::endl;
            streamlog_out(DEBUG7)<<"subTrack "<<i+1<<" rhoFirst: "<<rhoFirst<<" rhoLast: "<<rhoLast<<" loop forward: "<<loopForward<<std::endl;
        }

        int nHitsInFit = hitsInFit.size();
        // if first successfully fitted subTrack add IP track state
        if ( trackStates.empty() ){
            HitState state;
            state.ts = *(static_cast<const IMPL::TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtIP)) );
            // no hit/simhit for the extrapolated track state at the IP
            trackStates.push_back(state);
        }
        // NOTE: although we use z to understand subTrack's direction, subTrack's hits are still sorted by rho
        if (loopForward){
            for( int j=0; j<nHitsInFit; ++j ){
                HitState state;
                state.ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                state.hit = hitsInFit[j].first;
                state.simHit = getSimTrackerHit(state.hit, navToSimTrackerHits);
                streamlog_out(DEBUG6)<<"Added state from hit at ("<<state.hit->getPosition()[0]<<", "<<state.hit->getPosition()[1]<<", "<<state.hit->getPosition()[2]<<")"<<std::endl;
                streamlog_out(DEBUG6)<<state.ts<<std::endl;
                trackStates.push_back(state);
            }
        }
        else{
            for( int j=nHitsInFit-1; j>=0; --j ){
                HitState state;
                state.ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                state.hit = hitsInFit[j].first;
                state.simHit = getSimTrackerHit(state.hit, navToSimTrackerHits);
                streamlog_out(DEBUG6)<<"Added state from hit at ("<<state.hit->getPosition()[0]<<", "<<state.hit->getPosition()[1]<<", "<<state.hit->getPosition()[2]<<")"<<std::endl;
                streamlog_out(DEBUG6)<<state.ts<<std::endl;
                trackStates.push_back(state);
            }
        }
    }

    const EVENT::TrackState* tsCaloBugged = lastGoodRefittedTrack.getTrackState(TrackState::AtCalorimeter);
    if ( pfo->getClusters().size() > 0 && tsCaloBugged != nullptr ){
        // d0 and z0 of the track state at calo MUST be 0. This is a bug! Fix manually here.
        IMPL::TrackStateImpl tsCalo = *(dynamic_cast<const IMPL::TrackStateImpl*> (tsCaloBugged));
        tsCalo.setD0(0.);
        tsCalo.setZ0(0.);
        HitState state;
        state.ts = tsCalo;
        // no hit/simhit for the extrapolated track state at the ECAL
        streamlog_out(DEBUG7)<<"Added state tsCalo ("<<(state.ts).getReferencePoint()[0]<<", "<<(state.ts).getReferencePoint()[1]<<", "<<(state.ts).getReferencePoint()[2]<<")"<<std::endl;
        streamlog_out(DEBUG7)<<state.ts<<std::endl;
        trackStates.push_back( state );
    }
    return trackStates;
}

double getHelixLength(Vector3D p_start, double z_start, double z_end, double bField){
    double c_factor = 0.299792458;
    double pt = p_start.rho();
    double pz = p_start.z();
    return (z_end - z_start)/pz * std::sqrt(  std::pow( pt/(c_factor*bField), 2) + pz*pz  );
}

double getTrackLengthSHA(Track* track, int location=TrackState::AtCalorimeter, TrackLengthOption option=TrackLengthOption::zedLambda){
    const TrackState* tsIP = track->getTrackState( TrackState::AtIP );
    const TrackState* tsCalo = getTrackStateAtCalorimeter(track);
    return getHelixLength( tsIP, tsCalo, location, option );
}


TrackLengthResult getTrackLengthIKF(const std::vector<IMPL::TrackStateImpl>& trackStates, double bField, TrackLengthOption option){
    TrackLengthResult result;
    int nTrackStates = trackStates.size();
    streamlog_out(DEBUG7)<<"Calculating track length for"<<nTrackStates<<" track states"<<std::endl;
    if (nTrackStates <= 1) return result;

    for( int i=1; i < nTrackStates-1; ++i ){
        double arcLength = getHelixLength( &trackStates[i-1], &trackStates[i], TrackState::AtIP, option );
        std::array<double, 3> momArr = UTIL::getTrackMomentum( &(trackStates[i-1]), bField);
        Vector3D mom(momArr[0], momArr[1], momArr[2]);
        result.trackLengthToSET += arcLength;
        result.harmonicMomToSET += arcLength/mom.r2();
        streamlog_out(DEBUG7)<<"Arc length "<<i<<": "<<arcLength<<"  sum: "<<result.trackLengthToSET<<std::endl;
    }

    //now calculate to the Ecal one more step
    double lastArcLength = getHelixLength( &trackStates[nTrackStates - 2], &trackStates[nTrackStates - 1], TrackState::AtIP, option );
    std::array<double, 3> lastMomArr = UTIL::getTrackMomentum( &(trackStates[nTrackStates - 2]), bField );
    Vector3D lastMom(lastMomArr[0], lastMomArr[1], lastMomArr[2]);
    result.trackLengthToEcal = result.trackLengthToSET + lastArcLength;
    result.harmonicMomToEcal = result.harmonicMomToSET + lastArcLength/lastMom.r2();
    streamlog_out(DEBUG7)<<"Last Arc length "<<": "<<lastArcLength<<"  sum: "<<result.trackLengthToEcal<<std::endl;

    // don't forget to do the last step to properly calculate harmonic momentum!
    result.harmonicMomToSET = std::sqrt(result.trackLengthToSET/result.harmonicMomToSET);
    result.harmonicMomToEcal = std::sqrt(result.trackLengthToEcal/result.harmonicMomToEcal);

    return result;
}


////////////////////////////////////////////////////////////////
double getHelixNRevolutions(const TrackState& ts1, const TrackState& ts2){
    double omega = ts1.getOmega();
    double tanL = ts1.getTanLambda();
    double z1 = ts1.getReferencePoint()[2] + ts1.getZ0();
    double z2 = ts2.getReferencePoint()[2] + ts2.getZ0();

    // helix length projected on xy
    double circHelix = std::abs( (z2-z1)/tanL );
    double circFull = 2*M_PI/std::abs(omega);

    return circHelix/circFull;
}


double getTrackLengthIDR(Track* track){
    const TrackState* tsIP = track->getTrackState( TrackState::AtIP );
    const TrackState* tsCalo = track->getTrackState( TrackState::AtCalorimeter );
    
    float phiIP = tsIP->getPhi();
    float phiCalo = tsCalo->getPhi();
    float omega = tsIP->getOmega();
    float tanL = tsIP->getTanLambda();    
    return (phiIP-phiCalo)*(1/omega)*sqrt(1+tanL*tanL);
}

