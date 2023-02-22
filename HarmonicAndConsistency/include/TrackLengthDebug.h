#ifndef TrackLengthDebug_h
#define TrackLengthDebug_h 1

#include "marlin/Processor.h"
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "EVENT/Track.h"
#include "IMPL/TrackStateImpl.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "DD4hep/Detector.h"
#include "EVENT/SimTrackerHit.h"
#include <string>
#include <vector>
#include "EVENT/ReconstructedParticle.h"
#include "DDRec/Vector3D.h"
#include "EVENT/Track.h"
#include <memory>

struct Vars{
    double trackLength = 0.;
    double momentumHM = 0.;
    int nBadPhi = 0;
    int nBadZ = 0;
    int nSegments = 0;
};

class TrackLengthDebug : public marlin::Processor{
    public:
        TrackLengthDebug(const TrackLengthDebug&) = delete;
        TrackLengthDebug& operator=(const TrackLengthDebug&) = delete;

        marlin::Processor* newProcessor() { return new TrackLengthDebug; }

        TrackLengthDebug();
        void init();
        void processEvent(EVENT::LCEvent* evt);
        void end();
        Vars getTrackLengthUsingZ(EVENT::Track* track, double bField, MarlinTrk::IMarlinTrkSystem* trkSystem);

    private:
        int _nEvent{};
        double _bField{};
        
        std::unique_ptr<TFile> _file;
        std::unique_ptr<TTree> _tree;
        int _pdg;
        double _momentumIP;
        double _momentumCalo;
        double _tof;
        double _momentumHM;
        double _trackLength;
        int _nBadPhi;
        int _nBadZ;
        int _nSegments;
        MarlinTrk::IMarlinTrkSystem* _trkSystem = nullptr;
};

#endif
