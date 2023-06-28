#ifndef BohdanDrawing_h
#define BohdanDrawing_h 1

#include "BohdanUtils.h"

#include "EVENT/ReconstructedParticle.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/MCParticle.h"
#include "IMPL/TrackStateImpl.h"
#include "DDRec/Vector3D.h"
#include "TStyle.h"

TStyle* getMyStyle();
void displayPFO(EVENT::ReconstructedParticle* pfo, IMPL::TrackStateImpl tsStdReco, IMPL::TrackStateImpl tsEasy);
void displayFTDSimHits(EVENT::LCEvent* evt);
void plotCanvas(EVENT::Cluster* cluster, dd4hep::rec::Vector3D posAtEcal, dd4hep::rec::Vector3D momAtEcal, EVENT::MCParticle* mc);


#endif