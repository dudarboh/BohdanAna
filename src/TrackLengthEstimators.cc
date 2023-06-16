#include "TrackLengthEstimators.h"
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
using std::unique_ptr;
using dd4hep::rec::Vector3D;
using MarlinTrk::IMarlinTrack;
using MarlinTrk::IMarlinTrkSystem;
using CLHEP::RandGauss;


std::vector<EVENT::Track*> getSubTracks(EVENT::Track* track){
    vector<Track*> subTracks;
    // add track itself, which contains VXD+FTD+SIT+TPC hits of the first curl.
    subTracks.push_back(track);

    int nSubTracks = track->getTracks().size();
    if (nSubTracks <= 1) return subTracks;

    UTIL::BitField64 encoder( UTIL::LCTrackerCellID::encoding_string() ) ;
    auto isTPCHit = [&encoder](TrackerHit* hit) -> bool {
        encoder.setValue( hit->getCellID0() ) ;
        int subdet = encoder[ UTIL::LCTrackerCellID::subdet() ];
        return subdet == UTIL::ILDDetID::TPC;
    };

    int indexOfFirstTPCCurl = 0;
    for(int i = 0; i < nSubTracks; ++i){
        Track* subTrack = track->getTracks()[i];
        auto hits = subTrack->getTrackerHits();
        if ( std::find_if(hits.begin(), hits.end(), isTPCHit) != hits.end() ){
            indexOfFirstTPCCurl = i;
            break;
        }
    }

    for(int j=indexOfFirstTPCCurl+1; j < nSubTracks; ++j) subTracks.push_back( track->getTracks()[j] );
    return subTracks;
}

vector<TrackStateImpl> getTrackStates(ReconstructedParticle* pfo, double bField, MarlinTrk::IMarlinTrkSystem* trkSystem){
    vector<TrackStateImpl> trackStates;
    if ( pfo->getTracks().empty() ) return trackStates;
    vector<Track*> subTracks = getSubTracks( pfo->getTracks()[0] );

    // extract track state per every hit
    TrackImpl lastGoodRefittedTrack;

    auto getTrackStateAtHit = [](MarlinTrk::IMarlinTrack* marlinTrack, TrackerHit* hit) -> TrackStateImpl {
        TrackStateImpl ts;
        double chi2Dummy;
        int ndfDummy;
        marlinTrack->getTrackState(hit, ts, chi2Dummy, ndfDummy);
        return ts;
    };

    auto printStateLong = [](TrackStateImpl ts){
        auto part = ts;
        std::stringstream tmp;
        //out << std::scientific << std::setprecision (2) << std::showpos;
        std::cout << std::noshowpos;
        std::cout << std::setw(41) << std::setfill('-') << std::right << "-- TrackState ---" << std::setfill('-') << std::setw(29) << "-" << std::endl;

        tmp.str("") ;
        switch( part.getLocation() ){
        case EVENT::TrackState::AtOther         :   tmp <<  "AtOther"        ;   break ;
        case EVENT::TrackState::AtIP            :   tmp <<  "AtIP"           ;   break ;
        case EVENT::TrackState::AtFirstHit      :   tmp <<  "AtFirstHit"     ;   break ;
        case EVENT::TrackState::AtLastHit       :   tmp <<  "AtLastHit"      ;   break ;
        case EVENT::TrackState::AtCalorimeter   :   tmp <<  "AtCalorimeter " ;   break ;
        case EVENT::TrackState::AtVertex        :   tmp <<  "AtVertex"       ;   break ;
        }
        std::cout << std::setw(30) << std::setfill(' ') << std::left << "Location" << std::right << std::setw(40) << tmp.str() << std::endl;
        tmp.str("") ;
        tmp << std::dec << std::setfill('0') << std::setw(8) << part.id();
        std::cout << std::scientific << std::setprecision(6) ;
        std::cout << std::setw(30) << std::setfill(' ') << std::left << "Id"          << std::right << std::setw(40) << tmp.str() << std::endl;
        std::cout << std::setw(30) << std::setfill(' ') << std::left << "D0"          << std::right << std::setw(40) << part.getD0() << std::endl;
        std::cout << std::setw(30) << std::setfill(' ') << std::left << "Phi"         << std::right << std::setw(40) << part.getPhi() << std::endl;
        std::cout << std::setw(30) << std::setfill(' ') << std::left << "Omega"       << std::right << std::setw(40) << part.getOmega() << std::endl;
        std::cout << std::setw(30) << std::setfill(' ') << std::left << "Z0"          << std::right << std::setw(40) << part.getZ0() << std::endl;
        std::cout << std::setw(30) << std::setfill(' ') << std::left << "Tan Lambda"  << std::right << std::setw(40) << part.getTanLambda() << std::endl;
        tmp.str("");
        tmp  << std::dec << part.getReferencePoint()[0] << ", " << part.getReferencePoint()[1]  << ", " << part.getReferencePoint()[2]; 
        std::cout << std::setw(30) << std::setfill(' ') << std::left << "ReferencePoint" << std::right << std::setw(40) << tmp.str() << std::endl;
        std::cout << "Cov matrix:" << std::showpos << std::scientific << std::setprecision(6) << std::setw(15) << std::setfill(' ')  ;
        // print cov matrix as lower triangle matrix 
        for( unsigned l=0 , N=part.getCovMatrix().size(), ncolumns = 1 , nele =1 ; l <N ; ++l , ++nele) {
        std::cout << part.getCovMatrix()[l];
        if(! ( (nele) % ncolumns ) ){ 
            nele = 0 ;
            ++ncolumns ;
            std::cout << std::endl << "             " ;
        } else {
            std::cout << ", ";
        } 
        }
        std::cout << std::noshowpos;
        std::cout<<std::endl;
    };

    auto printState = [](TrackStateImpl ts){
        std::cout<<std::endl;
        std::cout<<"Location:    "<<ts.getLocation()<<std::endl;
        std::cout<<"D0 Z0:    "<<ts.getD0()<<"  "<<ts.getZ0()<<std::endl;
        std::cout<<"Phi:    "<<ts.getPhi()<<std::endl;
        std::cout<<"Omega:    "<<ts.getOmega()<<std::endl;
        std::cout<<"TanLambda:    "<<ts.getTanLambda()<<std::endl;
        std::cout<<"Ref:    "<<ts.getReferencePoint()[0]<<"  "<<ts.getReferencePoint()[1]<<"  "<<ts.getReferencePoint()[2]<<std::endl;
    };


    for(int i=0; i<subTracks.size(); ++i){
        vector <TrackerHit*> hits = subTracks[i]->getTrackerHits();

        auto sortByRho = [](TrackerHit* a, TrackerHit* b) -> bool {
            Vector3D posA( a->getPosition() ), posB( b->getPosition() );
            return posA.rho() < posB.rho();
        };

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
        unique_ptr<IMarlinTrack> marlinTrk( trkSystem->createTrack() );
        TrackImpl refittedTrack;

        //Need to initialize trackState at last hit
        TrackStateImpl preFit = *(subTracks[i]->getTrackState(TrackState::AtLastHit));
        preFit.setCovMatrix( covMatrix );
        bool fitDirection = IMarlinTrack::backward;
        std::cout<<"Calling MarlinTrk::createFinalisedLCIOTrack(); "<<std::endl;
        int errorFit = MarlinTrk::createFinalisedLCIOTrack(marlinTrk.get(), hits, &refittedTrack, fitDirection, &preFit, bField, maxChi2PerHit);
        //if fit fails, try also fit forward
        if (errorFit != 0){
            marlinTrk.reset( trkSystem->createTrack() );
            fitDirection = IMarlinTrack::forward;
            errorFit = MarlinTrk::createFinalisedLCIOTrack(marlinTrk.get(), hits, &refittedTrack, fitDirection, &preFit, bField, maxChi2PerHit);
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
        // if first subTrack add state at the IP
        if ( trackStates.empty() ){
            const TrackStateImpl* tsIP = static_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtIP));
            trackStates.push_back(*(tsIP));
        }
        if (loopForward){
            //add hits in increasing rho order
            for( int j=0; j<nHitsInFit; ++j ){
                TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                trackStates.push_back(ts);
            }
        }
        else{
            //in decreasing rho order
            for( int j=nHitsInFit-1; j>=0; --j ){
                TrackStateImpl ts = getTrackStateAtHit(marlinTrk.get(), hitsInFit[j].first);
                trackStates.push_back(ts);
            }
        }
    }

    const TrackStateImpl* tsCalo = static_cast<const TrackStateImpl*> (lastGoodRefittedTrack.getTrackState(TrackState::AtCalorimeter) );
    if ( pfo->getClusters().size() > 0 && tsCalo != nullptr ) trackStates.push_back( *(tsCalo) );
    return trackStates;

}


double getHelixNRevolutions(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2){
    double omega = ts1.getOmega();
    double tanL = ts1.getTanLambda();
    double z1 = ts1.getReferencePoint()[2] + ts1.getZ0();
    double z2 = ts2.getReferencePoint()[2] + ts2.getZ0();

    // helix length projected on xy
    double circHelix = std::abs( (z2-z1)/tanL );
    double circFull = 2*M_PI/std::abs(omega);

    return circHelix/circFull;
}


double getTrackLengthIDR(EVENT::Track* track){
    const TrackState* tsIP = track->getTrackState( TrackState::AtIP );
    const TrackState* tsCalo = track->getTrackState( TrackState::AtCalorimeter );
    
    float phiIP = tsIP->getPhi();
    float phiCalo = tsCalo->getPhi();
    float omega = tsIP->getOmega();
    float tanL = tsIP->getTanLambda();    
    return (phiIP-phiCalo)*(1/omega)*sqrt(1+tanL*tanL);
}

double getTrackLengthSHA1(EVENT::Track* track){
    const TrackState* tsIP = track->getTrackState( TrackState::AtIP );
    const TrackState* tsCalo = nullptr;
    if (track->getTracks().size() <= 2) {
        tsCalo = track->getTrackState( TrackState::AtCalorimeter );
    }
    else{
        Track* lastSubTrack = track->getTracks().back();
        tsCalo = lastSubTrack->getTrackState( TrackState::AtCalorimeter );
    }

    float omega = std::abs( tsIP->getOmega() );
    float tanL = std::abs( tsIP->getTanLambda() );
    double dPhi = std::abs( tsCalo->getPhi() - tsIP->getPhi() );
    if (dPhi > M_PI) dPhi = 2*M_PI - dPhi;
    return dPhi/omega*std::sqrt(1+tanL*tanL);
}

double getTrackLengthSHA2(EVENT::Track* track){
    const TrackState* tsIP = track->getTrackState( TrackState::AtIP );
    const TrackState* tsCalo = nullptr;
    if (track->getTracks().size() <= 2) {
        tsCalo = track->getTrackState( TrackState::AtCalorimeter );
    }
    else{
        Track* lastSubTrack = track->getTracks().back();
        tsCalo = lastSubTrack->getTrackState( TrackState::AtCalorimeter );
    }

    float tanL = std::abs( tsIP->getTanLambda() );

    double zIP = tsIP->getReferencePoint()[2] + tsIP->getZ0();
    double zCalo = tsCalo->getReferencePoint()[2] + tsCalo->getZ0();
    double dz = std::abs( zCalo - zIP);

    return dz/tanL*std::sqrt(1+tanL*tanL);
}

double getTrackLengthSHA3(EVENT::Track* track){
    const TrackState* tsIP = track->getTrackState( TrackState::AtIP );
    const TrackState* tsCalo = nullptr;
    if (track->getTracks().size() <= 2) {
        tsCalo = track->getTrackState( TrackState::AtCalorimeter );
    }
    else{
        Track* lastSubTrack = track->getTracks().back();
        tsCalo = lastSubTrack->getTrackState( TrackState::AtCalorimeter );
    }

    float omega = std::abs( tsIP->getOmega() );
    double zIP = tsIP->getReferencePoint()[2] + tsIP->getZ0();
    double zCalo = tsCalo->getReferencePoint()[2] + tsCalo->getZ0();
    double dz = std::abs( zCalo - zIP);

    double dPhi = std::abs( tsCalo->getPhi() - tsIP->getPhi() );
    if (dPhi > M_PI) dPhi = 2*M_PI - dPhi;
    return std::sqrt(dPhi*dPhi/(omega*omega)+dz*dz);
}

double getTrackLengthSHA4(EVENT::Track* track){
    const TrackState* tsIP = track->getTrackState( TrackState::AtIP );
    const TrackState* tsCalo = nullptr;
    if (track->getTracks().size() <= 2) {
        tsCalo = track->getTrackState( TrackState::AtCalorimeter );
    }
    else{
        Track* lastSubTrack = track->getTracks().back();
        tsCalo = lastSubTrack->getTrackState( TrackState::AtCalorimeter );
    }

    float omega = std::abs( tsCalo->getOmega() );
    float tanL = std::abs( tsCalo->getTanLambda() );
    double dPhi = std::abs( tsCalo->getPhi() - tsIP->getPhi() );
    if (dPhi > M_PI) dPhi = 2*M_PI - dPhi;
    return dPhi/omega*std::sqrt(1+tanL*tanL);
}

double getTrackLengthSHA5(EVENT::Track* track){
    const TrackState* tsIP = track->getTrackState( TrackState::AtIP );
    const TrackState* tsCalo = nullptr;

    if (track->getTracks().size() <= 2) {
        tsCalo = track->getTrackState( TrackState::AtCalorimeter );
    }
    else{
        Track* lastSubTrack = track->getTracks().back();
        tsCalo = lastSubTrack->getTrackState( TrackState::AtCalorimeter );
    }

    float tanL = std::abs( tsCalo->getTanLambda() );

    double zIP = tsIP->getReferencePoint()[2] + tsIP->getZ0();
    double zCalo = tsCalo->getReferencePoint()[2] + tsCalo->getZ0();
    double dz = std::abs( zCalo - zIP);

    return dz/tanL*std::sqrt(1+tanL*tanL);
}


double getTrackLengthSHA6(EVENT::Track* track){
    const TrackState* tsIP = track->getTrackState( TrackState::AtIP );
    const TrackState* tsCalo = nullptr;

    if (track->getTracks().size() <= 2) {
        tsCalo = track->getTrackState( TrackState::AtCalorimeter );
    }
    else{
        Track* lastSubTrack = track->getTracks().back();
        tsCalo = lastSubTrack->getTrackState( TrackState::AtCalorimeter );
    }

    float omega = std::abs( tsCalo->getOmega() );

    double zIP = tsIP->getReferencePoint()[2] + tsIP->getZ0();
    double zCalo = tsCalo->getReferencePoint()[2] + tsCalo->getZ0();
    double dz = std::abs( zCalo - zIP);

    double dPhi = std::abs( tsCalo->getPhi() - tsIP->getPhi() );
    if (dPhi > M_PI) dPhi = 2*M_PI - dPhi;
    return std::sqrt(dPhi*dPhi/(omega*omega)+dz*dz);
}


double getTrackLengthIKF1(std::vector<TrackStateImpl> trackStates){
    double trackLength = 0.;
    int nTrackStates = trackStates.size();
    if (nTrackStates <= 1) return trackLength;

    auto getHelixLength = [](const EVENT::TrackState& ts1, const EVENT::TrackState& ts2){
        float omega = std::abs( ts1.getOmega() );
        float tanL = std::abs( ts1.getTanLambda() );
        double dPhi = std::abs( ts2.getPhi() - ts1.getPhi() );
        if (dPhi > M_PI) dPhi = 2*M_PI - dPhi;
        return dPhi/omega*std::sqrt(1+tanL*tanL);
    };

    for( int j=1; j < nTrackStates; ++j ) trackLength += getHelixLength( trackStates[j-1], trackStates[j] );
    return trackLength;
}


double getTrackLengthIKF2(std::vector<TrackStateImpl> trackStates){
    double trackLength = 0.;
    int nTrackStates = trackStates.size();
    if (nTrackStates <= 1) return trackLength;

    auto getHelixLength = [](const EVENT::TrackState& ts1, const EVENT::TrackState& ts2){
        double tanL = ts1.getTanLambda();
        double z1 = ts1.getReferencePoint()[2] + ts1.getZ0();
        double z2 = ts2.getReferencePoint()[2] + ts2.getZ0();
        return std::abs( (z2-z1)/tanL ) * std::sqrt( 1.+tanL*tanL );
    };

    for( int j=1; j < nTrackStates; ++j ) trackLength += getHelixLength( trackStates[j-1], trackStates[j] );
    return trackLength;
}


double getTrackLengthIKF3(std::vector<TrackStateImpl> trackStates){
    double trackLength = 0.;
    int nTrackStates = trackStates.size();
    if (nTrackStates <= 1) return trackLength;

    auto getHelixLength = [](const EVENT::TrackState& ts1, const EVENT::TrackState& ts2){
        double omega = std::abs( ts1.getOmega() );
        double z1 = ts1.getReferencePoint()[2] + ts1.getZ0();
        double z2 = ts2.getReferencePoint()[2] + ts2.getZ0();
        double dz = std::abs( z2 - z1);

        double dPhi = std::abs( ts2.getPhi() - ts1.getPhi() );
        if (dPhi > M_PI) dPhi = 2*M_PI - dPhi;
        return std::sqrt(dPhi*dPhi/(omega*omega)+dz*dz);
    };

    for( int j=1; j < nTrackStates; ++j ) trackLength += getHelixLength( trackStates[j-1], trackStates[j] );
    return trackLength;
}
