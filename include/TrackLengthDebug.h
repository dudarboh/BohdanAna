#ifndef TrackLengthDebug_h
#define TrackLengthDebug_h 1

#include "marlin/Processor.h"
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "EVENT/Track.h"
#include "IMPL/TrackStateImpl.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "DD4hep/Detector.h"
#include "EVENT/SimTrackerHit.h"
#include <string>
#include <vector>
#include "EventDisplayer.h"

struct TOFResults{
    double tof0 = 0.;
    double tof10 = 0.;
    double tof20 = 0.;
    double tof30 = 0.;
    double tof40 = 0.;
    double tof50 = 0.;
    double tof60 = 0.;
    double tof70 = 0.;
    double tof80 = 0.;
    double tof90 = 0.;
    double tof100 = 0.;
    int layer = -1;
    double dToTrack = -1;
};


class TrackLengthDebug : public marlin::Processor, EventDisplayer {
    friend class EventDisplayer;
    public:
        TrackLengthDebug(const TrackLengthDebug&) = delete;
        TrackLengthDebug& operator=(const TrackLengthDebug&) = delete;

        marlin::Processor* newProcessor() { return new TrackLengthDebug; }

        TrackLengthDebug();
        void init();
        void processEvent(EVENT::LCEvent* evt);
        void end();
        void printout();
    private:
        int _nEvent{};
        double _bField{};
        
        std::unique_ptr<TFile> _file;
        std::unique_ptr<TTree> _tree;
        int _pdg;
        double _momentumIP;
        double _momentumCalo;
        double _momentumHM;

        TOFResults _tofClosest;

        double _trackLengthSHA1;
        double _trackLengthSHA2;
        double _trackLengthSHA3;
        double _trackLengthSHA4;
        double _trackLengthSHA5;
        double _trackLengthSHA6;

        double _trackLengthIKF1;
        double _trackLengthIKF2Bug;
        double _trackLengthIKF2;
        double _trackLengthIKF3;

        MarlinTrk::IMarlinTrkSystem* _trkSystem = nullptr;
};

#endif
