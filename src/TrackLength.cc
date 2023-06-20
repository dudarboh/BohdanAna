#include "TrackLength.h"
#include "BohdanUtils.h"

#include "marlin/VerbosityLevels.h"
#include "MarlinTrk/IMarlinTrack.h"

#include "MarlinTrk/MarlinTrkUtils.h"
#include "UTIL/TrackTools.h"
#include "IMPL/TrackImpl.h"

using namespace EVENT;
using dd4hep::rec::Vector3D;

std::vector<IMPL::TrackStateImpl> getTrackStates(ReconstructedParticle* pfo, double bField, MarlinTrk::IMarlinTrkSystem* trkSystem){
    // Refit the track and extract track state at every tracker hit along the track
    std::vector<IMPL::TrackStateImpl> trackStates;
    if ( pfo->getTracks().empty() ) return trackStates;
    std::vector<Track*> subTracks = getSubTracks( pfo->getTracks()[0] );

    IMPL::TrackImpl lastGoodRefittedTrack;

    auto sortByRho = [](TrackerHit* a, TrackerHit* b) -> bool {
        Vector3D posA( a->getPosition() ), posB( b->getPosition() );
        return posA.rho() < posB.rho();
    };

    streamlog_out(DEBUG8)<<"PFOs track has "<<subTracks.size()<<" subTracks."<<std::endl;
    for(size_t i=0; i<subTracks.size(); ++i){
        std::vector <TrackerHit*> hits = subTracks[i]->getTrackerHits();

        streamlog_out(DEBUG8)<<"Subtrack "<<i+1<<" has "<<hits.size()<<" hits."<<std::endl;
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
            streamlog_out(DEBUG8)<<"Fit backward fails! Trying to fit forward for "<<i+1<<" subTrack in this PFO!"<<std::endl;
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
        // if first successfully fitted subTrack add IP track state
        if ( trackStates.empty() ) trackStates.push_back(*(static_cast<const IMPL::TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtIP)) ));

        // NOTE: although we use z to understand subTrack's direction, subTrack's hits are still sorted by rho
        if (loopForward){
            for( int j=0; j<nHitsInFit; ++j ){
                IMPL::TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                trackStates.push_back(ts);
            }
        }
        else{
            for( int j=nHitsInFit-1; j>=0; --j ){
                IMPL::TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                trackStates.push_back(ts);
            }
        }
    }

    const IMPL::TrackStateImpl* tsCalo = static_cast<const IMPL::TrackStateImpl*> (lastGoodRefittedTrack.getTrackState(TrackState::AtCalorimeter) );
    if ( pfo->getClusters().size() > 0 && tsCalo != nullptr ) trackStates.push_back( *(tsCalo) );
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


TrackLengthResult getTrackLengthIKF(std::vector<IMPL::TrackStateImpl> trackStates, double bField, TrackLengthOption option){
    TrackLengthResult result;
    int nTrackStates = trackStates.size();
    if (nTrackStates <= 1) return result;

    for( int i=1; i < nTrackStates-1; ++i ){
        double arcLength = getHelixLength( &trackStates[i-1], &trackStates[i], TrackState::AtIP, option );
        std::array<double, 3> momArr = UTIL::getTrackMomentum( &(trackStates[i-1]), bField);
        Vector3D mom(momArr[0], momArr[1], momArr[2]);
        result.trackLengthToSET += arcLength;
        result.harmonicMomToSET += arcLength/mom.r2();
    }

    //now calculate to the Ecal one more step
    double lastArcLength = getHelixLength( &trackStates[nTrackStates - 2], &trackStates[nTrackStates - 1], TrackState::AtIP, option );
    std::array<double, 3> lastMomArr = UTIL::getTrackMomentum( &(trackStates[nTrackStates - 2]), bField );
    Vector3D lastMom(lastMomArr[0], lastMomArr[1], lastMomArr[2]);
    result.trackLengthToEcal = result.trackLengthToSET + lastArcLength;
    result.harmonicMomToEcal = result.harmonicMomToSET + lastArcLength/lastMom.r2();

    // don't forget to properly calculate harmonic momentum!
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

