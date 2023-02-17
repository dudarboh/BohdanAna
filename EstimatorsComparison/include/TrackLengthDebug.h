#ifndef TrackLengthDebug_h
#define TrackLengthDebug_h 1

#include "EventDisplayer.h"

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

class TrackLengthDebug : public marlin::Processor{
    public:
        TrackLengthDebug(const TrackLengthDebug&) = delete;
        TrackLengthDebug& operator=(const TrackLengthDebug&) = delete;

        marlin::Processor* newProcessor() { return new TrackLengthDebug; }

        TrackLengthDebug();
        void init();
        void processEvent(EVENT::LCEvent* evt);
        void end();
    private:
        int _nEvent{};
        double _bField{};
        
        std::unique_ptr<TFile> _file;
        std::unique_ptr<TTree> _tree;
        int _pdg;
        double _momentum;
        double _tof;
        double _trackLengthIDR;
        double _trackLengthIDR2;
        double _trackLengthIDR3;
        double _trackLengthIDR4;
        double _trackLengthWinni;
        double _trackLengthWinni2;
        double _trackLengthUsingZ;
        double _trackLengthSimUsingZ;

        MarlinTrk::IMarlinTrkSystem* _trkSystem = nullptr;
};

#endif
