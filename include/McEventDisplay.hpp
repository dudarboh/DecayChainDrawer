#ifndef TOFMcEventDisplay_h
#define TOFMcEventDisplay_h 1

#include "marlin/Processor.h"
#include "DD4hep/Detector.h"
#include "UTIL/LCRelationNavigator.h"
#include "EVENT/MCParticle.h"
#include "EVENT/Vertex.h"



class McEventDisplay : public marlin::Processor {
    public:
        McEventDisplay();
        marlin::Processor* newProcessor() {return new McEventDisplay;}
        void init();
        void processEvent(LCEvent* event);
        void end();

        void fillDecayChainUp(EVENT::MCParticle* mc, std::vector<EVENT::MCParticle*>& mcs);
        void fillDecayChainDown(EVENT::MCParticle* mc, std::vector<EVENT::MCParticle*>& mcs);
        EVENT::MCParticle* getMcMaxTrackWeight(EVENT::ReconstructedParticle* pfo, UTIL::LCRelationNavigator nav);
        std::vector<EVENT::MCParticle*> getVertexDecayChain(EVENT::Vertex* vertex, UTIL::LCRelationNavigator navRecoToMc);

        void createDotDiagram(LCEvent* event);
        std::map<int, std::string> getPdgNamesMap();

        void drawMcParticle(EVENT::MCParticle* mc);

        void debugPrint();


        dd4hep::Detector& _detector = dd4hep::Detector::getInstance();
        double _bField;
        //map of pointer to index in the collection
        std::map<MCParticle*, int> _p2idx;
        std::map<MCParticle*, int> _p2vtx;
        std::map<MCParticle*, double> _p2distance;
        std::map<MCParticle*, double> _p2pt;
        std::map<MCParticle*, double> _p2pz;
        std::map<int, std::string> _pdg2str;
        std::vector <std::string> _vtxColors = {"yellow", "yellow4", "yellowgreen", "orange", "orange4", "lightpink", "lightcoral", "lightcyan", "lightslateblue", "lightseagreen"};

        int _nEvent;
};


#endif
