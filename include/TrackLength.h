#ifndef TrackLength_h
#define TrackLength_h 1

#include "EVENT/ReconstructedParticle.h"
#include "EVENT/SimTrackerHit.h"
#include "IMPL/TrackStateImpl.h"
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "DDRec/Vector3D.h"
#include "EVENT/Track.h"
#include <vector>
#include "UTIL/LCRelationNavigator.h"

enum class TrackLengthOption{
    phiLambda,
    phiZed,
    zedLambda, // this one latest and greatest by a margin!
};

struct TrackLengthResult{
    double trackLengthToSET{};
    double harmonicMomToSET{};
    double trackLengthToEcal{};
    double harmonicMomToEcal{};
};

struct HitState{
    // Struct which contains a combination of the simulated, reconstructed tracker hit and the track state of the track at this reconstructed hit position
    IMPL::TrackStateImpl ts;
    EVENT::TrackerHit* hit = nullptr;
    EVENT::SimTrackerHit* simHit = nullptr;
};

template <typename TrackStateLike>
double getHelixLength(const TrackStateLike* ts1, const TrackStateLike* ts2, int location, TrackLengthOption option){
    float omega{}, tanL{};
    switch (location){
        case EVENT::TrackState::AtIP :
            omega = std::abs( ts1->getOmega() );
            tanL = std::abs( ts1->getTanLambda() );
        break;

        case EVENT::TrackState::AtCalorimeter :
            omega = std::abs( ts2->getOmega() );
            tanL = std::abs( ts2->getTanLambda() );
        break;
    }

    double dPhi = std::abs( ts2->getPhi() - ts1->getPhi() );
    if (dPhi > M_PI) dPhi = 2*M_PI - dPhi;

    double zIP = ts1->getReferencePoint()[2] + ts1->getZ0();
    double zCalo = ts2->getReferencePoint()[2] + ts2->getZ0();
    double dz = std::abs( zCalo - zIP );

    switch (option){
        case TrackLengthOption::phiLambda :
            return dPhi/omega*std::sqrt(1+tanL*tanL);

        case TrackLengthOption::phiZed :
            return std::sqrt(dPhi*dPhi/(omega*omega)+dz*dz);

        case TrackLengthOption::zedLambda :
            return dz/tanL*std::sqrt(1+tanL*tanL);
    }
    //something is wrong, return unphysical value
    return -1;
}


//theoretical formula integrating helix curve from z_start to z_end using momentum p_start at z_start.
double getHelixLength(dd4hep::rec::Vector3D p_start, double z_start, double z_end, double bField);

// simple helix approximation (SHA)
double getTrackLengthSHA(EVENT::Track* track, int location, TrackLengthOption option);


std::vector<HitState> getTrackStates(EVENT::ReconstructedParticle* pfo, double bField, MarlinTrk::IMarlinTrkSystem* trkSystem, const UTIL::LCRelationNavigator& navToSimTrackerHits);

// iterative Kalman Filter (IKF)
TrackLengthResult getTrackLengthIKF(const std::vector<IMPL::TrackStateImpl>& trackStates, double bField, TrackLengthOption option);


////////////////////////////////////////////////////////////////

double getHelixNRevolutions(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2);

//Track length from the IDR production version (1 September 2020)
// It doesn't account for the phi singularity producing sometimes wrong results
// It also doesn't take the absolute value and returns signed track length, which might be bug or a feature...
double getTrackLengthIDR(EVENT::Track* track);


#endif
