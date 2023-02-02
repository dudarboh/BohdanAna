#ifndef TrackLengthDebug_h
#define TrackLengthDebug_h 1

#include "EventDisplayer.h"

#include "marlin/Processor.h"
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "EVENT/Track.h"
#include "IMPL/TrackStateImpl.h"
#include "UTIL/PIDHandler.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "DD4hep/Detector.h"
#include "EVENT/SimTrackerHit.h"
#include <string>
#include <vector>

/**
Marlin processor that calculates harmonic mean momentum and track length of the track.
\author B. Dudar, DESY, 2022
*/
class TrackLengthDebug : public marlin::Processor, EventDisplayer {
    friend class EventDisplayer;
    public:
        /**
        Copy constructor.
        */
        TrackLengthDebug(const TrackLengthDebug&) = delete;

        /**
        Copy assignment operator.
        */
        TrackLengthDebug& operator=(const TrackLengthDebug&) = delete;

        /**
        Method required by the Marlin to register processor in the global scope.
        */
        marlin::Processor* newProcessor() { return new TrackLengthDebug; }

        /**
        Default constructor.
        Registers steering parameters from the xml steering file.
        */
        TrackLengthDebug();

        /** Called at the begin of the job before anything is read.
        Extracts geometry details and initializes Kalman Filter System.
        */
        void init();

        /** Called for every event - the working horse.
        Calculates momentum and track length and writes them into PIDHandler of the input collection object.
        */
        void processEvent(EVENT::LCEvent* evt);

        void end();


        double getTrackLengthDefault(const std::vector<IMPL::TrackStateImpl>& trackStates);
        double getTrackLengthTanL(const std::vector<IMPL::TrackStateImpl>& trackStates);
        double getTrackLengthZ(const std::vector<IMPL::TrackStateImpl>& trackStates);
        float getParameterFromPID(EVENT::ReconstructedParticle* pfo, UTIL::PIDHandler& pidHandler, std::string algorithmName, std::string parameterName);

        void prepareRootTree();
        // void drawPFO(EVENT::ReconstructedParticle* pfo);



    private:
        /** Stores ReconstructedParticleCollection steering parameter.
        */
        std::string _pfoCollectionName{};
        /** Stores current event number.
        */
        int _nEvent{};

        /** Stores names of the output parameters.
        These are "trackLengthToSET", "trackLengthToEcal", "momentumHMToSET", "momentumHMToEcal".
        */
        std::vector<std::string> _outputParNames{};


        /** Kalman Filter System.
        \note Release notes of MarlinTrk v02-00:
        \note users should no longer delete the IMarlinTrkSystem pointer in their code
        */
        MarlinTrk::IMarlinTrkSystem* _trkSystem = nullptr;

        /** Stores z component of the magnetic field at the origin in Tesla.
        */
        double _bField{};
        
        std::unique_ptr<TFile> _file;
        std::unique_ptr<TTree> _tree;
        double _momentum;
        double _tof;
        double _trackLengthDefault;
        double _trackLengthTanL;
        double _trackLengthZ;
        double _massDefault;
        double _massTanL;
        double _massZ;

        int _pdg;
        dd4hep::rec::Vector3D _momIp;
        double _d0Ip;
        double _z0Ip;
        double _omegaIp;
        double _tanLIp;
        double _phiIp;
        dd4hep::rec::Vector3D _momCalo;
        double _d0Calo;
        double _z0Calo;
        double _omegaCalo;
        double _tanLCalo;
        double _phiCalo;
        double _zCalo;

        TGraph* _hD0;
        TGraph* _hZ0;
        TGraph* _hOmega;
        TGraph* _hTanL;
        TGraph* _hPhi;
        TGraph* _hZ;
        TGraph* _hPt;
        TGraph* _hPz;
        TGraph* _hP;
        TGraph* _hTrkLenDefault;
        TGraph* _hTrkLenTanL;
        TGraph* _hTrkLenZ;

        TGraph* _hMcOmega;
        TGraph* _hMcTanL;
        TGraph* _hMcPhi;
        TGraph* _hMcZ;
        TGraph* _hMcPt;
        TGraph* _hMcPz;
        TGraph* _hMcP;
        TGraph* _hMcTrkLenDefault;
        TGraph* _hMcTrkLenTanL;
        TGraph* _hMcTrkLenZ;

};

#endif
