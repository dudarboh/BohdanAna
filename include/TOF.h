#ifndef TOF_h
#define TOF_h 1

#include "EVENT/CalorimeterHit.h"
#include "EVENT/Cluster.h"
#include "EVENT/Track.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/MCParticle.h"
#include "DDRec/Vector3D.h"
#include <vector>
#include <utility>

std::vector<EVENT::CalorimeterHit*> selectFrankEcalHits( EVENT::Cluster* cluster, dd4hep::rec::Vector3D posAtEcal, dd4hep::rec::Vector3D momAtEcal, int maxEcalLayer);

std::pair<int, double> getTofClosest( EVENT::Cluster* cluster, dd4hep::rec::Vector3D posAtEcal, double timeResolution);

double getTofFrankAvg( const std::vector<EVENT::CalorimeterHit*>& selectedHits, dd4hep::rec::Vector3D posAtEcal, double timeResolution);

double getTofFrankFit( const std::vector<EVENT::CalorimeterHit*>& selectedHits, dd4hep::rec::Vector3D posAtEcal, double timeResolution);

double getTofSET(EVENT::Track* track, double timeResolution);

double getTofPhotonTrue(EVENT::MCParticle* mc);


#endif