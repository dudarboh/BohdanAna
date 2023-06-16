#include "BohdanDrawing.h"

#include "EVENT/TrackState.h"
#include "UTIL/TrackTools.h"
#include "marlinutil/DDMarlinCED.h"


void drawPFO(EVENT::ReconstructedParticle* pfo, IMPL::TrackStateImpl tsStdReco, IMPL::TrackStateImpl tsEasy){
    std::vector<EVENT::Track*> tracks;
    if ( not pfo->getTracks().empty() )
        tracks = getSubTracks(pfo->getTracks()[0]);
    for(auto* track: tracks){
        auto hits = track->getTrackerHits();
        for (auto* hit : hits){
            auto pos = hit->getPosition();
            int type = 0; // point
            int size = 4; // larger point
            unsigned long color = 0x000000;
            ced_hit(pos[0], pos[1], pos[2], type, size, color);
        }
    }
    //plot ecal state for fun
    const TrackState* tsCalo = tracks.back()->getTrackState( TrackState::AtCalorimeter );
    auto posCalo = tsCalo->getReferencePoint();
    ced_hit(posCalo[0], posCalo[1], posCalo[2] + tsCalo->getZ0(), 0, 4, 0x000000); // track state at calorimeterer


    std::vector<EVENT::Cluster*> clusters = pfo->getClusters();
    for(auto* cluster: clusters){
        auto hits = cluster->getCalorimeterHits();
        for (auto* hit: hits){
            auto pos = hit->getPosition();
            int type = 0; // point
            int size = 4; // larger point
            unsigned long color = 0x000000;
            ced_hit(pos[0], pos[1], pos[2], type, size, color);
        }
    }

    auto plotHelixFromTrackState = [](TrackStateImpl ts, unsigned long color){
        double omega = ts.getOmega();
        // double tanL = ts.getTanLambda();
        double phi = ts.getPhi();
        double d0 = ts.getD0();
        double z0 = ts.getZ0();
        auto ref = ts.getReferencePoint();
        std::array<double, 3> mom = UTIL::getTrackMomentum(&ts, 3.5);

        //
        int charge = omega > 0 ? 1 : -1;
        double x = -d0*std::sin(phi) + ref[0];
        double y = d0*std::cos(phi) + ref[1];
        double z = z0 + ref[2];


        int type = 0; // point
        int size = 8; // larger point
        ced_hit(x, y, z, type, size, color);
        size = 2;
        DDMarlinCED::drawHelix( 3.5, charge, x, y, z, mom[0], mom[1], mom[2], type, size, color, 0, 2000, 2600, 0);
        DDMarlinCED::drawHelix( -3.5, -charge, x, y, z, mom[0], mom[1], mom[2], type, size, color, 0, 2000, 2600, 0);
    };

    plotHelixFromTrackState(tsStdReco, 0xfc0505);
    printTrackStateLong(tsStdReco);
    plotHelixFromTrackState(tsEasy, 0x0905fc);
    printTrackStateLong(tsEasy);

}
