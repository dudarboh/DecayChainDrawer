#ifndef TOFDecayChainDrawer_h
#define TOFDecayChainDrawer_h 1

#include "marlin/Processor.h"
#include "DD4hep/Detector.h"
#include "UTIL/LCRelationNavigator.h"
#include "EVENT/MCParticle.h"
#include "EVENT/Vertex.h"



class DecayChainDrawer : public marlin::Processor {
    public:
        DecayChainDrawer();
        marlin::Processor* newProcessor() {return new DecayChainDrawer;}
        void init();
        void processEvent(LCEvent* event);
        void end(){};

        void fillDecayChainUp(EVENT::MCParticle* mc, std::vector<EVENT::MCParticle*>& mcs);
        void fillDecayChainDown(EVENT::MCParticle* mc, std::vector<EVENT::MCParticle*>& mcs);
        EVENT::MCParticle* getMcMaxTrackWeight(EVENT::ReconstructedParticle* pfo, UTIL::LCRelationNavigator nav);
        std::vector<EVENT::MCParticle*> getVertexDecayChain(EVENT::Vertex* vertex, UTIL::LCRelationNavigator navRecoToMc);
        std::map<int, std::string> getPdgNamesMap();
        
        bool isInHadronization(EVENT::MCParticle* mc);

        //map of pointer to index in the collection
        std::map<MCParticle*, int> _p2idx;
        std::map<MCParticle*, int> _p2vtx;
        std::map<MCParticle*, bool> _p2hadronization;
        std::map<MCParticle*, double> _p2distance;
        std::map<MCParticle*, double> _p2pt;
        std::map<MCParticle*, double> _p2pz;
        std::map<int, std::string> _pdg2str;
        std::vector <std::string> _vtxColors = {"yellow", "yellow4", "yellowgreen", "orange", "orange4", "lightpink", "lightcoral", "lightcyan", "lightslateblue", "lightseagreen"};

        int _nEvent{};
};


#endif
