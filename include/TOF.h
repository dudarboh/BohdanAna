#ifndef TOF_h
#define TOF_h 1

#include "EVENT/CalorimeterHit.h"
#include "EVENT/Cluster.h"
#include "EVENT/Track.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/MCParticle.h"
#include "UTIL/LCRelationNavigator.h"
#include "DDRec/Vector3D.h"
#include <vector>
#include <utility>

std::vector<EVENT::CalorimeterHit*> selectFrankEcalHits( EVENT::Cluster* cluster, const dd4hep::rec::Vector3D& posAtEcal, const dd4hep::rec::Vector3D& momAtEcal, int maxEcalLayer );

EVENT::CalorimeterHit* getClosestHit( EVENT::Cluster* cluster, const dd4hep::rec::Vector3D& posAtEcal );

float getHitTof( EVENT::CalorimeterHit* hit, const dd4hep::rec::Vector3D& posAtEcal, float timeResolution) ;
int getHitCaloType( EVENT::CalorimeterHit* hit );
int getHitCaloID( EVENT::CalorimeterHit* hit );
int getHitCaloLayout( EVENT::CalorimeterHit* hit );
int getHitCaloLayer( EVENT::CalorimeterHit* hit );
EVENT::MCParticle* getHitEarliestMC( EVENT::CalorimeterHit* hit, const UTIL::LCRelationNavigator& navToSimCalorimeterHits );

float getTofFrankAvg( const std::vector<EVENT::CalorimeterHit*>& selectedHits, const dd4hep::rec::Vector3D& posAtEcal, float timeResolution );

float getTofFrankFit( const std::vector<EVENT::CalorimeterHit*>& selectedHits, const dd4hep::rec::Vector3D& posAtEcal, float timeResolution );

std::tuple<float, float> getTofSET(EVENT::Track* track, float timeResolution);

float getTofPhotonTrue(EVENT::MCParticle* mc );


#endif