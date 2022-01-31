#include "DecayChainDrawer.hpp"
#include "ColorMap.hpp"

#include "marlinutil/DDMarlinCED.h"
#include "marlinutil/GeometryUtil.h"
#include "marlinutil/MarlinUtil.h"
#include "DDRec/Vector3D.h"

using namespace std;
using dd4hep::rec::Vector3D;

DecayChainDrawer aDecayChainDrawer;

DecayChainDrawer::DecayChainDrawer() : Processor("DecayChainDrawer"){}

void DecayChainDrawer::init(){
    _pdg2str = getPdgNamesMap();
}


void DecayChainDrawer::processEvent(LCEvent* event){
    std::cout<<++_nEvent<<std::endl;
    _p2idx.clear();
    _p2vtx.clear();
    _p2hadronization.clear();
    _p2distance.clear();
    _p2pt.clear();
    _p2pz.clear();

    LCCollection* mcCol = event->getCollection("MCParticle");
    LCCollection* vertices = event->getCollection("BuildUpVertex");
    LCRelationNavigator navRecoToMc( event->getCollection("RecoMCTruthLink") );

    MCParticle* mc0 = static_cast<MCParticle*> (mcCol->getElementAt(0));
    Vector3D ip( mc0->getVertex() );
    int nVertices = vertices->getNumberOfElements();
    if (nVertices == 0) return;


    // Fill maps
    for(int i=0; i< mcCol->getNumberOfElements(); ++i){
        MCParticle* mc = static_cast<MCParticle*> (mcCol->getElementAt(i));
        _p2idx[mc] = i;
        _p2distance[mc] = (Vector3D( mc->getVertex() ) - ip).r();
        _p2vtx[mc] = 0;
        _p2pt[mc] = Vector3D(mc->getMomentum()).trans();
        _p2pz[mc] = Vector3D(mc->getMomentum()).z();

        // is MCParticle in the hadronization decay?
        _p2hadronization[mc] = isInHadronization(mc);

        for(int j=0; j<nVertices; ++j){
            Vertex* vertex = static_cast<Vertex*> (vertices->getElementAt(j));
            vector<MCParticle*> vertexDecayChain = getVertexDecayChain(vertex, navRecoToMc);
            bool foundInTheVertexDecay = std::find(vertexDecayChain.begin(), vertexDecayChain.end(), mc) != vertexDecayChain.end();
            if (foundInTheVertexDecay) _p2vtx[mc] = j+1;
        }
    }

    // draw all MC Particles with their relations:
    std::stringstream nodes;
    std::stringstream labels;
    for(int i=0; i < mcCol->getNumberOfElements(); ++i){
        MCParticle* mc = static_cast<MCParticle*> (mcCol->getElementAt(i));
        if (! ( _p2vtx[mc] !=0 || (mc->getGeneratorStatus() != 0 && _p2hadronization[mc] )) ) continue;

        for( auto daughter : mc->getDaughters() ){
            if (! ( _p2vtx[daughter] !=0 || (daughter->getGeneratorStatus() != 0 && _p2hadronization[daughter] )) ) continue;
            nodes<<"    "<<_p2idx[mc]<<"->"<<_p2idx[daughter]<<";"<<endl;
        }

        int pdg = mc->getPDG();
        bool hasName = _pdg2str.find(pdg) != _pdg2str.end();
        string label;
        if (hasName) label = _pdg2str[pdg];
        else label = std::to_string(pdg);

        labels<<_p2idx[mc]<<"[label=<"<<label<<"<BR/>"<<std::fixed<<std::setprecision(2)<<_p2distance[mc]<<" mm<BR/>"<<_p2pt[mc]<<" | "<<_p2pz[mc]<<" GeV"<<">";
        if (_p2vtx[mc] != 0 ) labels<<" style=\"filled\" fillcolor=\""<<_vtxColors[_p2vtx[mc]-1]<<"\"";
        labels<<"];"<<endl;
    }

    //create dot file with a header
    std::ofstream outfile;
    outfile.open("test.dot");
    outfile<<"digraph {"<<endl;
    outfile<<"    rankdir=TB;"<<endl;
    outfile<<nodes.str()<<endl;
    outfile<<labels.str()<<endl;

    outfile<<"}"<<endl;
    system("rm -f test.svg && dot -Tsvg test.dot > test.svg");
    system("xdg-open test.svg");  
}


void DecayChainDrawer::fillDecayChainUp(EVENT::MCParticle* mc, std::vector<EVENT::MCParticle*>& decayChain){
    decayChain.push_back(mc);
    // stop iterating up at hadronization
    if(mc->getPDG() == 92) return;
    const vector<MCParticle*> parents = mc->getParents();
    for(auto parent : parents){
        bool foundInTheList = std::find(decayChain.begin(), decayChain.end(), parent) != decayChain.end();
        if ( !foundInTheList ) fillDecayChainUp(parent, decayChain);
    }
}

void DecayChainDrawer::fillDecayChainDown(EVENT::MCParticle* mc, std::vector<EVENT::MCParticle*>& decayChain){
    // stop iterating down at particles not created by generator
    // if(mc->getGeneratorStatus() == 0) return;
    decayChain.push_back(mc);
    const vector<MCParticle*> daughters = mc->getDaughters();
    for(auto daughter : daughters){
        bool foundInTheList = std::find(decayChain.begin(), decayChain.end(), daughter) != decayChain.end();
        if ( !foundInTheList ) fillDecayChainDown(daughter, decayChain);
    }
}

EVENT::MCParticle* DecayChainDrawer::getMcMaxTrackWeight(EVENT::ReconstructedParticle* pfo, UTIL::LCRelationNavigator nav){
    const vector<LCObject*>& mcs = nav.getRelatedToObjects(pfo);
    const vector<float>& weights = nav.getRelatedToWeights(pfo);
    //get index of highest TRACK weight MC particle
    int i = std::max_element(weights.begin(), weights.end(), [](float a, float b){return (int(a)%10000)/1000. < (int(b)%10000)/1000.;}) - weights.begin();
    return static_cast<MCParticle*> ( mcs[i] );
}


std::vector<EVENT::MCParticle*> DecayChainDrawer::getVertexDecayChain(EVENT::Vertex* vertex, UTIL::LCRelationNavigator navRecoToMc){
    vector<MCParticle*> decayChain;
    std::vector <ReconstructedParticle*> pfos = vertex->getAssociatedParticle()->getParticles();
    for(auto pfo : pfos){
        MCParticle* mc = getMcMaxTrackWeight(pfo, navRecoToMc);
        fillDecayChainUp(mc, decayChain);
    }
    return decayChain;    
}


bool DecayChainDrawer::isInHadronization(EVENT::MCParticle* mc){
    if ( mc->getPDG() == 92 ) return true;
    const vector<MCParticle*> parents = mc->getParents();
    for(auto parent : parents){
        if ( isInHadronization(parent) ) return true;
    }
    return false;
}



std::map<int, std::string> DecayChainDrawer::getPdgNamesMap(){
    std::map<int, std::string> pdg2str;

    pdg2str[92] = "Hadronization";
    // leptons
    pdg2str[11] = "e<SUP>-</SUP>";
    pdg2str[-11] = "e<SUP>+</SUP>";
    pdg2str[12] = "&nu;<SUB>e</SUB>";
    pdg2str[-12] = "&nu;<SUB>e</SUB>";
    pdg2str[13] = "&mu;<SUP>-</SUP>";
    pdg2str[-13] = "&mu;<SUP>+</SUP>";
    pdg2str[14] = "&nu;<SUB>&mu;</SUB>";
    pdg2str[-14] = "&nu;<SUB>&mu;</SUB>";
    pdg2str[15] = "&tau;<SUP>-</SUP>";
    pdg2str[-15] = "&tau;<SUP>+</SUP>";
    pdg2str[16] = "&nu;<SUB>&tau;</SUB>";
    pdg2str[-16] = "&nu;<SUB>&tau;</SUB>";
    //Gauge and Higgs Boson
    pdg2str[22] = "&gamma;";
    pdg2str[23] = "Z<SUP>0</SUP>";
    pdg2str[24] = "W<SUP>+</SUP>";
    pdg2str[-24] = "W<SUP>-</SUP>";
    // Higgs?
    // light mesons I=1
    pdg2str[111] = "&pi;<SUP>0</SUP>";
    pdg2str[-111] = "&pi;<SUP>0</SUP>";
    pdg2str[211] = "&pi;<SUP>+</SUP>";
    pdg2str[-211] = "&pi;<SUP>-</SUP>";
    pdg2str[9000111] = "a<SUB>0</SUB>(980)<SUP>0</SUP>";
    pdg2str[-9000111] = "a<SUB>0</SUB>(980)<SUP>0</SUP>";
    pdg2str[9000211] = "a<SUB>0</SUB>(980)<SUP>+</SUP>";
    pdg2str[-9000211] = "a<SUB>0</SUB>(980)<SUP>-</SUP>";
    pdg2str[100111] = "&pi;(1300)<SUP>0</SUP>";
    pdg2str[-100111] = "&pi;(1300)<SUP>0</SUP>";
    pdg2str[100211] = "&pi;(1300)<SUP>+</SUP>";
    pdg2str[-100211] = "&pi;(1300)<SUP>-</SUP>";
    pdg2str[10111] = "a<SUB>0</SUB>(1450)<SUP>0</SUP>";
    pdg2str[-10111] = "a<SUB>0</SUB>(1450)<SUP>0</SUP>";
    pdg2str[10211] = "a<SUB>0</SUB>(1450)<SUP>+</SUP>";
    pdg2str[-10211] = "a<SUB>0</SUB>(1450)<SUP>-</SUP>";
    pdg2str[9010111] = "&pi;(1800)<SUP>0</SUP>";
    pdg2str[-9010111] = "&pi;(1800)<SUP>0</SUP>";
    pdg2str[9010211] = "&pi;(1800)<SUP>+</SUP>";
    pdg2str[-9010211] = "&pi;(1800)<SUP>-</SUP>";
    pdg2str[113] = "&rho;(770)<SUP>0</SUP>";
    pdg2str[-113] = "&rho;(770)<SUP>0</SUP>";
    pdg2str[213] = "&rho;(770)<SUP>+</SUP>";
    pdg2str[-213] = "&rho;(770)<SUP>-</SUP>";
    pdg2str[10113] = "b<SUB>1</SUB>(1235)<SUP>0</SUP>";
    pdg2str[-10113] = "b<SUB>1</SUB>(1235)<SUP>0</SUP>";
    pdg2str[10213] = "b<SUB>1</SUB>(1235)<SUP>+</SUP>";
    pdg2str[-10213] = "b<SUB>1</SUB>(1235)<SUP>-</SUP>";
    pdg2str[20113] = "a<SUB>1</SUB>(1260)<SUP>0</SUP>";
    pdg2str[-20113] = "a<SUB>1</SUB>(1260)<SUP>0</SUP>";
    pdg2str[20213] = "a<SUB>1</SUB>(1260)<SUP>+</SUP>";
    pdg2str[-20213] = "a<SUB>1</SUB>(1260)<SUP>-</SUP>";
    pdg2str[9000113] = "&pi;<SUB>1</SUB>(1400)<SUP>0</SUP>";
    pdg2str[-9000113] = "&pi;<SUB>1</SUB>(1400)<SUP>0</SUP>";
    pdg2str[9000213] = "&pi;<SUB>1</SUB>(1400)<SUP>+</SUP>";
    pdg2str[-9000213] = "&pi;<SUB>1</SUB>(1400)<SUP>-</SUP>";
    pdg2str[100113] = "&rho;(1450)<SUP>0</SUP>";
    pdg2str[-100113] = "&rho;(1450)<SUP>0</SUP>";
    pdg2str[100213] = "&rho;(1450)<SUP>+</SUP>";
    pdg2str[-100213] = "&rho;(1450)<SUP>-</SUP>";
    pdg2str[9010113] = "&pi;<SUB>1</SUB>(1600)<SUP>0</SUP>";
    pdg2str[-9010113] = "&pi;<SUB>1</SUB>(1600)<SUP>0</SUP>";
    pdg2str[9010213] = "&pi;<SUB>1</SUB>(1600)<SUP>+</SUP>";
    pdg2str[-9010213] = "&pi;<SUB>1</SUB>(1600)<SUP>-</SUP>";
    pdg2str[9020113] = "a<SUB>1</SUB>(1640)<SUP>0</SUP>";
    pdg2str[-9020113] = "a<SUB>1</SUB>(1640)<SUP>0</SUP>";
    pdg2str[9020213] = "a<SUB>1</SUB>(1640)<SUP>+</SUP>";
    pdg2str[-9020213] = "a<SUB>1</SUB>(1640)<SUP>-</SUP>";
    pdg2str[30113] = "&rho;(1700)<SUP>0</SUP>";
    pdg2str[-30113] = "&rho;(1700)<SUP>0</SUP>";
    pdg2str[30213] = "&rho;(1700)<SUP>+</SUP>";
    pdg2str[-30213] = "&rho;(1700)<SUP>-</SUP>";
    pdg2str[9030113] = "&rho;(1900)<SUP>0</SUP>";
    pdg2str[-9030113] = "&rho;(1900)<SUP>0</SUP>";
    pdg2str[9030213] = "&rho;(1900)<SUP>+</SUP>";
    pdg2str[-9030213] = "&rho;(1900)<SUP>-</SUP>";
    pdg2str[9040113] = "&rho;(2150)<SUP>0</SUP>";
    pdg2str[-9040113] = "&rho;(2150)<SUP>0</SUP>";
    pdg2str[9040213] = "&rho;(2150)<SUP>+</SUP>";
    pdg2str[-9040213] = "&rho;(2150)<SUP>-</SUP>";
    pdg2str[115] = "a<SUB>2</SUB>(1320)<SUP>0</SUP>";
    pdg2str[-115] = "a<SUB>2</SUB>(1320)<SUP>0</SUP>";
    pdg2str[215] = "a<SUB>2</SUB>(1320)<SUP>+</SUP>";
    pdg2str[-215] = "a<SUB>2</SUB>(1320)<SUP>-</SUP>";
    pdg2str[10115] = "&pi;<SUB>2</SUB>(1670)<SUP>0</SUP>";
    pdg2str[-10115] = "&pi;<SUB>2</SUB>(1670)<SUP>0</SUP>";
    pdg2str[10215] = "&pi;<SUB>2</SUB>(1670)<SUP>+</SUP>";
    pdg2str[-10215] = "&pi;<SUB>2</SUB>(1670)<SUP>-</SUP>";
    pdg2str[9000115] = "a<SUB>2</SUB>(1700)<SUP>0</SUP>";
    pdg2str[-9000115] = "a<SUB>2</SUB>(1700)<SUP>0</SUP>";
    pdg2str[9000215] = "a<SUB>2</SUB>(1700)<SUP>+</SUP>";
    pdg2str[-9000215] = "a<SUB>2</SUB>(1700)<SUP>-</SUP>";
    pdg2str[9010115] = "&pi;<SUB>2</SUB>(2100)<SUP>0</SUP>";
    pdg2str[-9010115] = "&pi;<SUB>2</SUB>(2100)<SUP>0</SUP>";
    pdg2str[9010215] = "&pi;<SUB>2</SUB>(2100)<SUP>+</SUP>";
    pdg2str[-9010215] = "&pi;<SUB>2</SUB>(2100)<SUP>-</SUP>";
    pdg2str[117] = "&rho;<SUB>3</SUB>(1690)<SUP>0</SUP>";
    pdg2str[-117] = "&rho;<SUB>3</SUB>(1690)<SUP>0</SUP>";
    pdg2str[217] = "&rho;<SUB>3</SUB>(1690)<SUP>+</SUP>";
    pdg2str[-217] = "&rho;<SUB>3</SUB>(1690)<SUP>-</SUP>";
    pdg2str[9000117] = "&rho;<SUB>3</SUB>(1990)<SUP>0</SUP>";
    pdg2str[-9000117] = "&rho;<SUB>3</SUB>(1990)<SUP>0</SUP>";
    pdg2str[9000217] = "&rho;<SUB>3</SUB>(1990)<SUP>+</SUP>";
    pdg2str[-9000217] = "&rho;<SUB>3</SUB>(1990)<SUP>-</SUP>";
    pdg2str[9010117] = "&rho;<SUB>3</SUB>(2250)<SUP>0</SUP>";
    pdg2str[-9010117] = "&rho;<SUB>3</SUB>(2250)<SUP>0</SUP>";
    pdg2str[9010217] = "&rho;<SUB>3</SUB>(2250)<SUP>+</SUP>";
    pdg2str[-9010217] = "&rho;<SUB>3</SUB>(2250)<SUP>-</SUP>";
    pdg2str[119] = "a<SUB>4</SUB>(2040)<SUP>0</SUP>";
    pdg2str[-119] = "a<SUB>4</SUB>(2040)<SUP>0</SUP>";
    pdg2str[219] = "a<SUB>4</SUB>(2040)<SUP>+</SUP>";
    pdg2str[-219] = "a<SUB>4</SUB>(2040)<SUP>-</SUP>";


    // light mesons I=0
    pdg2str[221] = "&eta;";
    pdg2str[-221] = "&eta;";
    pdg2str[331] = "&eta;\'";
    pdg2str[-331] = "&eta;\'";
    pdg2str[9000221] = "f<SUB>0</SUB>(600)";
    pdg2str[-9000221] = "f<SUB>0</SUB>(600)";
    pdg2str[9010221] = "f<SUB>0</SUB>(980)";
    pdg2str[-9010221] = "f<SUB>0</SUB>(980)";
    pdg2str[100221] = "&eta;(1295)";
    pdg2str[-100221] = "&eta;(1295)";
    pdg2str[10221] = "f<SUB>0</SUB>(1370)";
    pdg2str[-10221] = "f<SUB>0</SUB>(1370)";
    pdg2str[9020221] = "&eta;(1405)";
    pdg2str[-9020221] = "&eta;(1405)";
    pdg2str[100331] = "&eta;(1475)";
    pdg2str[-100331] = "&eta;(1475)";
    pdg2str[9030221] = "f<SUB>0</SUB>(1500)";
    pdg2str[-9030221] = "f<SUB>0</SUB>(1500)";
    pdg2str[10331] = "f<SUB>0</SUB>(1710)";
    pdg2str[-10331] = "f<SUB>0</SUB>(1710)";
    pdg2str[9040221] = "&eta;(1760)";
    pdg2str[-9040221] = "&eta;(1760)";
    pdg2str[9050221] = "f<SUB>0</SUB>(2020)";
    pdg2str[-9050221] = "f<SUB>0</SUB>(2020)";
    pdg2str[9060221] = "f<SUB>0</SUB>(2100)";
    pdg2str[-9060221] = "f<SUB>0</SUB>(2100)";
    pdg2str[9070221] = "f<SUB>0</SUB>(2200)";
    pdg2str[-9070221] = "f<SUB>0</SUB>(2200)";
    pdg2str[9080221] = "&eta;(2225)";
    pdg2str[-9080221] = "&eta;(2225)";
    pdg2str[223] = "&omega;";
    pdg2str[-223] = "&omega;";
    pdg2str[333] = "&phi;";
    pdg2str[-333] = "&phi;";
    pdg2str[10223] = "h<SUB>1</SUB>(1170)";
    pdg2str[-10223] = "h<SUB>1</SUB>(1170)";
    pdg2str[20223] = "f<SUB>1</SUB>(1285)";
    pdg2str[-20223] = "f<SUB>1</SUB>(1285)";
    pdg2str[10333] = "h<SUB>1</SUB>(1380)";
    pdg2str[-10333] = "h<SUB>1</SUB>(1380)";
    pdg2str[20333] = "f<SUB>1</SUB>(1420)";
    pdg2str[-20333] = "f<SUB>1</SUB>(1420)";
    pdg2str[100223] = "&omega;(1420)";
    pdg2str[-100223] = "&omega;(1420)";
    pdg2str[9000223] = "f<SUB>1</SUB>(1510)";
    pdg2str[-9000223] = "f<SUB>1</SUB>(1510)";
    pdg2str[9010223] = "h<SUB>1</SUB>(1595)";
    pdg2str[-9010223] = "h<SUB>1</SUB>(1595)";
    pdg2str[30223] = "&omega;(1650)";
    pdg2str[-30223] = "&omega;(1650)";
    pdg2str[100333] = "&phi;(1680)";
    pdg2str[-100333] = "&phi;(1680)";
    pdg2str[225] = "f<SUB>2</SUB>(1270)";
    pdg2str[-225] = "f<SUB>2</SUB>(1270)";
    pdg2str[9000225] = "f<SUB>2</SUB>(1430)";
    pdg2str[-9000225] = "f<SUB>2</SUB>(1430)";
    pdg2str[335] = "f<SUB>2</SUB>\'(1525)";
    pdg2str[-335] = "f<SUB>2</SUB>\'(1525)";
    pdg2str[9010225] = "f<SUB>2</SUB>(1565)";
    pdg2str[-9010225] = "f<SUB>2</SUB>(1565)";
    pdg2str[9020225] = "f<SUB>2</SUB>(1640)";
    pdg2str[-9020225] = "f<SUB>2</SUB>(1640)";
    pdg2str[10225] = "&eta;<SUB>2</SUB>(1645)";
    pdg2str[-10225] = "&eta;<SUB>2</SUB>(1645)";
    pdg2str[9030225] = "f<SUB>2</SUB>(1810)";
    pdg2str[-9030225] = "f<SUB>2</SUB>(1810)";
    pdg2str[10335] = "&eta;<SUB>2</SUB>(1870)";
    pdg2str[-10335] = "&eta;<SUB>2</SUB>(1870)";
    pdg2str[9040225] = "f<SUB>2</SUB>(1910)";
    pdg2str[-9040225] = "f<SUB>2</SUB>(1910)";
    pdg2str[9050225] = "f<SUB>2</SUB>(1950)";
    pdg2str[-9050225] = "f<SUB>2</SUB>(1950)";
    pdg2str[9060225] = "f<SUB>2</SUB>(2010)";
    pdg2str[-9060225] = "f<SUB>2</SUB>(2010)";
    pdg2str[9070225] = "f<SUB>2</SUB>(2150)";
    pdg2str[-9070225] = "f<SUB>2</SUB>(2150)";
    pdg2str[9080225] = "f<SUB>2</SUB>(2300)";
    pdg2str[-9080225] = "f<SUB>2</SUB>(2300)";
    pdg2str[9090225] = "f<SUB>2</SUB>(2340)";
    pdg2str[-9090225] = "f<SUB>2</SUB>(2340)";
    pdg2str[227] = "&omega;<SUB>3</SUB>(1670)";
    pdg2str[-227] = "&omega;<SUB>3</SUB>(1670)";
    pdg2str[337] = "&phi;<SUB>3</SUB>(1850)";
    pdg2str[-337] = "&phi;<SUB>3</SUB>(1850)";
    pdg2str[229] = "f<SUB>4</SUB>(2050)";
    pdg2str[-229] = "f<SUB>4</SUB>(2050)";
    pdg2str[9000229] = "f<SUB>J</SUB>(2220)";
    pdg2str[-9000229] = "f<SUB>J</SUB>(2220)";
    pdg2str[9010229] = "f<SUB>4</SUB>(2300)";
    pdg2str[-9010229] = "f<SUB>4</SUB>(2300)";
    
    // strange mesons
    pdg2str[130] = "K<SUB>L</SUB><SUP>0</SUP>";
    pdg2str[-130] = "K<SUB>L</SUB><SUP>0</SUP>";
    pdg2str[310] = "K<SUB>S</SUB><SUP>0</SUP>";
    pdg2str[-310] = "K<SUB>S</SUB><SUP>0</SUP>";
    pdg2str[311] = "K<SUP>0</SUP>";
    pdg2str[-311] = "K<SUP>0</SUP>";
    pdg2str[321] = "K<SUP>+</SUP>";
    pdg2str[-321] = "K<SUP>-</SUP>";
    pdg2str[9000311] = "K<SUB>0</SUB><SUP>*</SUP>(800)<SUP>0</SUP>";
    pdg2str[-9000311] = "K<SUB>0</SUB><SUP>*</SUP>(800)<SUP>0</SUP>";
    pdg2str[9000321] = "K<SUB>0</SUB><SUP>*</SUP>(800)<SUP>+</SUP>";
    pdg2str[-9000321] = "K<SUB>0</SUB><SUP>*</SUP>(800)<SUP>-</SUP>";
    pdg2str[10311] = "K<SUB>0</SUB><SUP>*</SUP>(1430)<SUP>0</SUP>";
    pdg2str[-10311] = "K<SUB>0</SUB><SUP>*</SUP>(1430)<SUP>0</SUP>";
    pdg2str[10321] = "K<SUB>0</SUB><SUP>*</SUP>(1430)<SUP>+</SUP>";
    pdg2str[-10321] = "K<SUB>0</SUB><SUP>*</SUP>(1430)<SUP>-</SUP>";
    pdg2str[100311] = "K(1460)<SUP>0</SUP>";
    pdg2str[-100311] = "K(1460)<SUP>0</SUP>";
    pdg2str[100321] = "K(1460)<SUP>+</SUP>";
    pdg2str[-100321] = "K(1460)<SUP>-</SUP>";
    pdg2str[9010311] = "K(1830)<SUP>0</SUP>";
    pdg2str[-9010311] = "K(1830)<SUP>0</SUP>";
    pdg2str[9010321] = "K(1830)<SUP>+</SUP>";
    pdg2str[-9010321] = "K(1830)<SUP>-</SUP>";
    pdg2str[9020311] = "K<SUB>0</SUB><SUP>*</SUP>(1950)<SUP>0</SUP>";
    pdg2str[-9020311] = "K<SUB>0</SUB><SUP>*</SUP>(1950)<SUP>0</SUP>";
    pdg2str[9020321] = "K<SUB>0</SUB><SUP>*</SUP>(1950)<SUP>+</SUP>";
    pdg2str[-9020321] = "K<SUB>0</SUB><SUP>*</SUP>(1950)<SUP>-</SUP>";
    pdg2str[313] = "K<SUP>*</SUP>(892)<SUP>0</SUP>";
    pdg2str[-313] = "K<SUP>*</SUP>(892)<SUP>0</SUP>";
    pdg2str[323] = "K<SUP>*</SUP>(892)<SUP>+</SUP>";
    pdg2str[-323] = "K<SUP>*</SUP>(892)<SUP>-</SUP>";
    pdg2str[10313] = "K<SUB>1</SUB>(1270)<SUP>0</SUP>";
    pdg2str[-10313] = "K<SUB>1</SUB>(1270)<SUP>0</SUP>";
    pdg2str[10323] = "K<SUB>1</SUB>(1270)<SUP>+</SUP>";
    pdg2str[-10323] = "K<SUB>1</SUB>(1270)<SUP>-</SUP>";
    pdg2str[20313] = "K<SUB>1</SUB>(1400)<SUP>0</SUP>";
    pdg2str[-20313] = "K<SUB>1</SUB>(1400)<SUP>0</SUP>";
    pdg2str[20323] = "K<SUB>1</SUB>(1400)<SUP>+</SUP>";
    pdg2str[-20323] = "K<SUB>1</SUB>(1400)<SUP>-</SUP>";
    pdg2str[100313] = "K<SUP>*</SUP>(1410)<SUP>0</SUP>";
    pdg2str[-100313] = "K<SUP>*</SUP>(1410)<SUP>0</SUP>";
    pdg2str[100323] = "K<SUP>*</SUP>(1410)<SUP>+</SUP>";
    pdg2str[-100323] = "K<SUP>*</SUP>(1410)<SUP>-</SUP>";
    pdg2str[9000313] = "K<SUB>1</SUB>(1650)<SUP>0</SUP>";
    pdg2str[-9000313] = "K<SUB>1</SUB>(1650)<SUP>0</SUP>";
    pdg2str[9000323] = "K<SUB>1</SUB>(1650)<SUP>+</SUP>";
    pdg2str[-9000323] = "K<SUB>1</SUB>(1650)<SUP>-</SUP>";
    pdg2str[30313] = "K<SUP>*</SUP>(1680)<SUP>0</SUP>";
    pdg2str[-30313] = "K<SUP>*</SUP>(1680)<SUP>0</SUP>";
    pdg2str[30323] = "K<SUP>*</SUP>(1680)<SUP>+</SUP>";
    pdg2str[-30323] = "K<SUP>*</SUP>(1680)<SUP>-</SUP>";
    pdg2str[315] = "K<SUB>2</SUB><SUP>*</SUP>(1430)<SUP>0</SUP>";
    pdg2str[-315] = "K<SUB>2</SUB><SUP>*</SUP>(1430)<SUP>0</SUP>";
    pdg2str[325] = "K<SUB>2</SUB><SUP>*</SUP>(1430)<SUP>+</SUP>";
    pdg2str[-325] = "K<SUB>2</SUB><SUP>*</SUP>(1430)<SUP>-</SUP>";
    pdg2str[9000315] = "K<SUB>2</SUB>(1580)<SUP>0</SUP>";
    pdg2str[-9000315] = "K<SUB>2</SUB>(1580)<SUP>0</SUP>";
    pdg2str[9000325] = "K<SUB>2</SUB>(1580)<SUP>+</SUP>";
    pdg2str[-9000325] = "K<SUB>2</SUB>(1580)<SUP>-</SUP>";
    pdg2str[10315] = "K<SUB>2</SUB>(1770)<SUP>0</SUP>";
    pdg2str[-10315] = "K<SUB>2</SUB>(1770)<SUP>0</SUP>";
    pdg2str[10325] = "K<SUB>2</SUB>(1770)<SUP>+</SUP>";
    pdg2str[-10325] = "K<SUB>2</SUB>(1770)<SUP>-</SUP>";
    pdg2str[20315] = "K<SUB>2</SUB>(1820)<SUP>0</SUP>";
    pdg2str[-20315] = "K<SUB>2</SUB>(1820)<SUP>0</SUP>";
    pdg2str[20325] = "K<SUB>2</SUB>(1820)<SUP>+</SUP>";
    pdg2str[-20325] = "K<SUB>2</SUB>(1820)<SUP>-</SUP>";
    pdg2str[9010315] = "K<SUB>2</SUB><SUP>*</SUP>(1980)<SUP>0</SUP>";
    pdg2str[-9010315] = "K<SUB>2</SUB><SUP>*</SUP>(1980)<SUP>0</SUP>";
    pdg2str[9010325] = "K<SUB>2</SUB><SUP>*</SUP>(1980)<SUP>+</SUP>";
    pdg2str[-9010325] = "K<SUB>2</SUB><SUP>*</SUP>(1980)<SUP>-</SUP>";
    pdg2str[9020315] = "K<SUB>2</SUB>(2250)<SUP>0</SUP>";
    pdg2str[-9020315] = "K<SUB>2</SUB>(2250)<SUP>0</SUP>";
    pdg2str[9020325] = "K<SUB>2</SUB>(2250)<SUP>+</SUP>";
    pdg2str[-9020325] = "K<SUB>2</SUB>(2250)<SUP>-</SUP>";
    pdg2str[317] = "K<SUB>3</SUB><SUP>*</SUP>(1780)<SUP>0</SUP>";
    pdg2str[-317] = "K<SUB>3</SUB><SUP>*</SUP>(1780)<SUP>0</SUP>";
    pdg2str[327] = "K<SUB>3</SUB><SUP>*</SUP>(1780)<SUP>+</SUP>";
    pdg2str[-327] = "K<SUB>3</SUB><SUP>*</SUP>(1780)<SUP>-</SUP>";
    pdg2str[9010317] = "K<SUB>3</SUB>(2320)<SUP>0</SUP>";
    pdg2str[-9010317] = "K<SUB>3</SUB>(2320)<SUP>0</SUP>";
    pdg2str[9010327] = "K<SUB>3</SUB>(2320)<SUP>+</SUP>";
    pdg2str[-9010327] = "K<SUB>3</SUB>(2320)<SUP>-</SUP>";
    pdg2str[319] = "K<SUB>4</SUB><SUP>*</SUP>(2045)<SUP>0</SUP>";
    pdg2str[-319] = "K<SUB>4</SUB><SUP>*</SUP>(2045)<SUP>0</SUP>";
    pdg2str[329] = "K<SUB>4</SUB><SUP>*</SUP>(2045)<SUP>+</SUP>";
    pdg2str[-329] = "K<SUB>4</SUB><SUP>*</SUP>(2045)<SUP>-</SUP>";
    pdg2str[9000319] = "K<SUB>4</SUB>(2500)<SUP>0</SUP>";
    pdg2str[-9000319] = "K<SUB>4</SUB>(2500)<SUP>0</SUP>";
    pdg2str[9000329] = "K<SUB>4</SUB>(2500)<SUP>+</SUP>";
    pdg2str[-9000329] = "K<SUB>4</SUB>(2500)<SUP>-</SUP>";


    //charmed mesons
    pdg2str[411] = "D<SUP>+</SUP>";
    pdg2str[-411] = "D<SUP>-</SUP>";
    pdg2str[421] = "D<SUP>0</SUP>";
    pdg2str[-421] = "D<SUP>0</SUP>";
    pdg2str[10411] = "D<SUB>0</SUB><SUP>*</SUP>(2400)<SUP>+</SUP>";
    pdg2str[-10411] = "D<SUB>0</SUB><SUP>*</SUP>(2400)<SUP>-</SUP>";
    pdg2str[10421] = "D<SUB>0</SUB><SUP>*</SUP>(2400)<SUP>0</SUP>";
    pdg2str[-10421] = "D<SUB>0</SUB><SUP>*</SUP>(2400)<SUP>0</SUP>";
    pdg2str[413] = "D<SUP>*</SUP>(2010)<SUP>+</SUP>";
    pdg2str[-413] = "D<SUP>*</SUP>(2010)<SUP>-</SUP>";
    pdg2str[423] = "D<SUP>*</SUP>(2007)<SUP>0</SUP>";
    pdg2str[-423] = "D<SUP>*</SUP>(2007)<SUP>0</SUP>";
    pdg2str[10413] = "D<SUB>1</SUB>(2420)<SUP>+</SUP>";
    pdg2str[-10413] = "D<SUB>1</SUB>(2420)<SUP>-</SUP>";
    pdg2str[10423] = "D<SUB>1</SUB>(2420)<SUP>0</SUP>";
    pdg2str[-10423] = "D<SUB>1</SUB>(2420)<SUP>0</SUP>";
    pdg2str[20413] = "D<SUB>1</SUB>(H)<SUP>+</SUP>";
    pdg2str[-20413] = "D<SUB>1</SUB>(H)<SUP>-</SUP>";
    pdg2str[20423] = "D<SUB>1</SUB>(2430)<SUP>0</SUP>";
    pdg2str[-20423] = "D<SUB>1</SUB>(2430)<SUP>0</SUP>";
    pdg2str[415] = "D<SUB>2</SUB><SUP>*</SUP>(2460)<SUP>+</SUP>";
    pdg2str[-415] = "D<SUB>2</SUB><SUP>*</SUP>(2460)<SUP>-</SUP>";
    pdg2str[425] = "D<SUB>2</SUB><SUP>*</SUP>(2460)<SUP>0</SUP>";
    pdg2str[-425] = "D<SUB>2</SUB><SUP>*</SUP>(2460)<SUP>0</SUP>";
    pdg2str[431] = "D<SUB>s</SUB><SUP>+</SUP>";
    pdg2str[-431] = "D<SUB>s</SUB><SUP>-</SUP>";
    pdg2str[10431] = "D<SUB>s0</SUB><SUP>*</SUP>(2317)<SUP>+</SUP>";
    pdg2str[-10431] = "D<SUB>s0</SUB><SUP>*</SUP>(2317)<SUP>-</SUP>";
    pdg2str[433] = "D<SUB>s</SUB><SUP>*+</SUP>";
    pdg2str[-433] = "D<SUB>s</SUB><SUP>*-</SUP>";
    pdg2str[10433] = "D<SUB>s1</SUB>(2536)<SUP>+</SUP>";
    pdg2str[-10433] = "D<SUB>s1</SUB>(2536)<SUP>-</SUP>";
    pdg2str[20433] = "D<SUB>s1</SUB>(2460)<SUP>+</SUP>";
    pdg2str[-20433] = "D<SUB>s1</SUB>(2460)<SUP>-</SUP>";
    pdg2str[435] = "D<SUB>s2</SUB><SUP>*</SUP>(2573)<SUP>+</SUP>";
    pdg2str[-435] = "D<SUB>s2</SUB><SUP>*</SUP>(2573)<SUP>-</SUP>";

    //bottom mesons
    pdg2str[511] = "B<SUP>0</SUP>";
    pdg2str[-511] = "B<SUP>0</SUP>";
    pdg2str[521] = "B<SUP>+</SUP>";
    pdg2str[-521] = "B<SUP>-</SUP>";
    pdg2str[10511] = "B<SUB>0</SUB><SUP>*0</SUP>";
    pdg2str[-10511] = "B<SUB>0</SUB><SUP>*0</SUP>";
    pdg2str[10521] = "B<SUB>0</SUB><SUP>*+</SUP>";
    pdg2str[-10521] = "B<SUB>0</SUB><SUP>*-</SUP>";
    pdg2str[513] = "B<SUP>*0</SUP>";
    pdg2str[-513] = "B<SUP>*0</SUP>";
    pdg2str[523] = "B<SUP>*+</SUP>";
    pdg2str[-523] = "B<SUP>*-</SUP>";
    pdg2str[10513] = "B<SUB>1</SUB>(L)<SUP>0</SUP>";
    pdg2str[-10513] = "B<SUB>1</SUB>(L)<SUP>0</SUP>";
    pdg2str[10523] = "B<SUB>1</SUB>(L)<SUP>+</SUP>";
    pdg2str[-10523] = "B<SUB>1</SUB>(L)<SUP>-</SUP>";
    pdg2str[20513] = "B<SUB>1</SUB>(H)<SUP>0</SUP>";
    pdg2str[-20513] = "B<SUB>1</SUB>(H)<SUP>0</SUP>";
    pdg2str[20523] = "B<SUB>1</SUB>(H)<SUP>+</SUP>";
    pdg2str[-20523] = "B<SUB>1</SUB>(H)<SUP>-</SUP>";
    pdg2str[515] = "B<SUB>2</SUB><SUP>*0</SUP>";
    pdg2str[-515] = "B<SUB>2</SUB><SUP>*0</SUP>";
    pdg2str[525] = "B<SUB>2</SUB><SUP>*+</SUP>";
    pdg2str[-525] = "B<SUB>2</SUB><SUP>*-</SUP>";
    pdg2str[531] = "B<SUB>s</SUB><SUP>0</SUP>";
    pdg2str[-531] = "B<SUB>s</SUB><SUP>0</SUP>";
    pdg2str[10531] = "B<SUB>s0</SUB><SUP>*0</SUP>";
    pdg2str[-10531] = "B<SUB>s0</SUB><SUP>*0</SUP>";
    pdg2str[533] = "B<SUB>s</SUB><SUP>*0</SUP>";
    pdg2str[-533] = "B<SUB>s</SUB><SUP>*0</SUP>";
    pdg2str[10533] = "B<SUB>s1</SUB>(L)<SUP>0</SUP>";
    pdg2str[-10533] = "B<SUB>s1</SUB>(L)<SUP>0</SUP>";
    pdg2str[20533] = "B<SUB>s1</SUB>(H)<SUP>0</SUP>";
    pdg2str[-20533] = "B<SUB>s1</SUB>(H)<SUP>0</SUP>";
    pdg2str[535] = "B<SUB>s2</SUB><SUP>*0</SUP>";
    pdg2str[-535] = "B<SUB>s2</SUB><SUP>*0</SUP>";
    pdg2str[541] = "B<SUB>c</SUB><SUP>+</SUP>";
    pdg2str[-541] = "B<SUB>c</SUB><SUP>-</SUP>";
    pdg2str[10541] = "B<SUB>c0</SUB><SUP>*+</SUP>";
    pdg2str[-10541] = "B<SUB>c0</SUB><SUP>*-</SUP>";
    pdg2str[543] = "B<SUB>c</SUB><SUP>*+</SUP>";
    pdg2str[-543] = "B<SUB>c</SUB><SUP>*-</SUP>";
    pdg2str[10543] = "B<SUB>c1</SUB>(L)<SUP>+</SUP>";
    pdg2str[-10543] = "B<SUB>c1</SUB>(L)<SUP>-</SUP>";
    pdg2str[20543] = "B<SUB>c1</SUB>(H)<SUP>+</SUP>";
    pdg2str[-20543] = "B<SUB>c1</SUB>(H)<SUP>-</SUP>";
    pdg2str[545] = "B<SUB>c2</SUB><SUP>+</SUP>";
    pdg2str[-545] = "B<SUB>c2</SUB><SUP>-</SUP>";

    //cc mesons
    pdg2str[441] = "&eta;<SUB>c</SUB>(1S)";
    pdg2str[-441] = "&eta;<SUB>c</SUB>(1S)";
    pdg2str[10441] = "&chi;<SUB>c0</SUB>(1P)";
    pdg2str[-10441] = "&chi;<SUB>c0</SUB>(1P)";
    pdg2str[100441] = "&eta;<SUB>c</SUB>(2S)";
    pdg2str[-100441] = "&eta;<SUB>c</SUB>(2S)";
    pdg2str[443] = "J/&psi;(1S)";
    pdg2str[-443] = "J/&psi;(1S)";
    pdg2str[10443] = "h<SUB>c</SUB>(1P)";
    pdg2str[-10443] = "h<SUB>c</SUB>(1P)";
    pdg2str[20443] = "&chi;<SUB>c1</SUB>(1P)";
    pdg2str[-20443] = "&chi;<SUB>c1</SUB>(1P)";
    pdg2str[100443] = "&psi;(2S)";
    pdg2str[-100443] = "&psi;(2S)";
    pdg2str[30443] = "&psi;(3770)";
    pdg2str[-30443] = "&psi;(3770)";
    pdg2str[9000443] = "&psi;(4040)";
    pdg2str[-9000443] = "&psi;(4040)";
    pdg2str[9010443] = "&psi;(4160)";
    pdg2str[-9010443] = "&psi;(4160)";
    pdg2str[9020443] = "&psi;(4415)";
    pdg2str[-9020443] = "&psi;(4415)";
    pdg2str[445] = "&chi;<SUB>c2</SUB>(1P)";
    pdg2str[-445] = "&chi;<SUB>c2</SUB>(1P)";
    pdg2str[100445] = "&chi;<SUB>c2</SUB>(2P)";
    pdg2str[-100445] = "&chi;<SUB>c2</SUB>(2P)";

    //light baryons
    pdg2str[2212] = "p<SUP>+</SUP>";
    pdg2str[-2212] = "p<SUP>-</SUP>";
    pdg2str[2112] = "n";
    pdg2str[-2112] = "n";
    pdg2str[2224] = "&Delta;<SUP>++</SUP>";
    pdg2str[-2224] = "&Delta;<SUP>--</SUP>";
    pdg2str[2214] = "&Delta;<SUP>+</SUP>";
    pdg2str[-2214] = "&Delta;<SUP>+</SUP>";
    pdg2str[2114] = "&Delta;<SUP>0</SUP>";
    pdg2str[-2114] = "&Delta;<SUP>0</SUP>";
    pdg2str[1114] = "&Delta;<SUP>-</SUP>";
    pdg2str[-1114] = "&Delta;<SUP>-</SUP>";
    //strange baryons
    pdg2str[3122] = "&Lambda;";
    pdg2str[-3122] = "&Lambda;";
    pdg2str[3222] = "&Sigma;<SUP>+</SUP>";
    pdg2str[-3222] = "&Sigma;<SUP>+</SUP>";
    pdg2str[3212] = "&Sigma;<SUP>0</SUP>";
    pdg2str[-3212] = "&Sigma;<SUP>0</SUP>";
    pdg2str[3112] = "&Sigma;<SUP>-</SUP>";
    pdg2str[-3112] = "&Sigma;<SUP>-</SUP>";
    pdg2str[3224] = "&Sigma;<SUP>*+</SUP>";
    pdg2str[-3224] = "&Sigma;<SUP>*+</SUP>";
    pdg2str[3214] = "&Sigma;<SUP>*0</SUP>";
    pdg2str[-3214] = "&Sigma;<SUP>*0</SUP>";
    pdg2str[3114] = "&Sigma;<SUP>*-</SUP>";
    pdg2str[-3114] = "&Sigma;<SUP>*-</SUP>";
    pdg2str[3322] = "&Xi;<SUP>0</SUP>";
    pdg2str[-3322] = "&Xi;<SUP>0</SUP>";
    pdg2str[3312] = "&Xi;<SUP>-</SUP>";
    pdg2str[-3312] = "&Xi;<SUP>-</SUP>";
    pdg2str[3324] = "&Xi;<SUP>*0</SUP>";
    pdg2str[-3324] = "&Xi;<SUP>*0</SUP>";
    pdg2str[3314] = "&Xi;<SUP>*-</SUP>";
    pdg2str[-3314] = "&Xi;<SUP>*-</SUP>";
    pdg2str[3334] = "&Omega;<SUP>-</SUP>";
    pdg2str[-3334] = "&Omega;<SUP>-</SUP>";
    //charmed baryons
    pdg2str[4122] = "&Lambda;<SUB>c</SUB><SUP>+</SUP>";
    pdg2str[-4122] = "&Lambda;<SUB>c</SUB><SUP>-</SUP>";
    pdg2str[4222] = "&Sigma;<SUB>c</SUB><SUP>++</SUP>";
    pdg2str[-4222] = "&Sigma;<SUB>c</SUB><SUP>--</SUP>";
    pdg2str[4212] = "&Sigma;<SUB>c</SUB><SUP>+</SUP>";
    pdg2str[-4212] = "&Sigma;<SUB>c</SUB><SUP>-</SUP>";
    pdg2str[4112] = "&Sigma;<SUB>c</SUB><SUP>0</SUP>";
    pdg2str[-4112] = "&Sigma;<SUB>c</SUB><SUP>0</SUP>";
    pdg2str[4224] = "&Sigma;<SUB>c</SUB><SUP>*++</SUP>";
    pdg2str[-4224] = "&Sigma;<SUB>c</SUB><SUP>*--</SUP>";
    pdg2str[4214] = "&Sigma;<SUB>c</SUB><SUP>*+</SUP>";
    pdg2str[-4214] = "&Sigma;<SUB>c</SUB><SUP>*-</SUP>";
    pdg2str[4114] = "&Sigma;<SUB>c</SUB><SUP>*0</SUP>";
    pdg2str[-4114] = "&Sigma;<SUB>c</SUB><SUP>*0</SUP>";
    pdg2str[4232] = "&Xi;<SUB>c</SUB><SUP>+</SUP>";
    pdg2str[-4232] = "&Xi;<SUB>c</SUB><SUP>-</SUP>";
    pdg2str[4132] = "&Xi;<SUB>c</SUB><SUP>0</SUP>";
    pdg2str[-4132] = "&Xi;<SUB>c</SUB><SUP>0</SUP>";
    pdg2str[4322] = "&Xi;<SUB>c</SUB><SUP>\'+</SUP>";
    pdg2str[-4322] = "&Xi;<SUB>c</SUB><SUP>\'-</SUP>";
    pdg2str[4312] = "&Xi;<SUB>c</SUB><SUP>\'0</SUP>";
    pdg2str[-4312] = "&Xi;<SUB>c</SUB><SUP>\'0</SUP>";
    pdg2str[4324] = "&Xi;<SUB>c</SUB><SUP>*+</SUP>";
    pdg2str[-4324] = "&Xi;<SUB>c</SUB><SUP>*-</SUP>";
    pdg2str[4314] = "&Xi;<SUB>c</SUB><SUP>*0</SUP>";
    pdg2str[-4314] = "&Xi;<SUB>c</SUB><SUP>*0</SUP>";
    pdg2str[4332] = "&Omega;<SUB>c</SUB><SUP>0</SUP>";
    pdg2str[-4332] = "&Omega;<SUB>c</SUB><SUP>0</SUP>";
    pdg2str[4334] = "&Omega;<SUB>c</SUB><SUP>*0</SUP>";
    pdg2str[-4334] = "&Omega;<SUB>c</SUB><SUP>*0</SUP>";
    pdg2str[4412] = "&Xi;<SUB>cc</SUB><SUP>+</SUP>";
    pdg2str[-4412] = "&Xi;<SUB>cc</SUB><SUP>-</SUP>";
    pdg2str[4422] = "&Xi;<SUB>cc</SUB><SUP>++</SUP>";
    pdg2str[-4422] = "&Xi;<SUB>cc</SUB><SUP>--</SUP>";
    pdg2str[4414] = "&Xi;<SUB>cc</SUB><SUP>*+</SUP>";
    pdg2str[-4414] = "&Xi;<SUB>cc</SUB><SUP>*-</SUP>";
    pdg2str[4424] = "&Xi;<SUB>cc</SUB><SUP>*++</SUP>";
    pdg2str[-4424] = "&Xi;<SUB>cc</SUB><SUP>*--</SUP>";
    pdg2str[4432] = "&Omega;<SUB>cc</SUB><SUP>+</SUP>";
    pdg2str[-4432] = "&Omega;<SUB>cc</SUB><SUP>-</SUP>";
    pdg2str[4434] = "&Omega;<SUB>cc</SUB><SUP>*+</SUP>";
    pdg2str[-4434] = "&Omega;<SUB>cc</SUB><SUP>*-</SUP>";
    pdg2str[4444] = "&Omega;<SUB>ccc</SUB><SUP>++</SUP>";
    pdg2str[-4444] = "&Omega;<SUB>ccc</SUB><SUP>--</SUP>";
    //bottom baryons
    pdg2str[5122] = "&Lambda;<SUB>b</SUB><SUP>0</SUP>";
    pdg2str[-5122] = "&Lambda;<SUB>b</SUB><SUP>0</SUP>";
    pdg2str[5112] = "&Sigma;<SUB>b</SUB><SUP>-</SUP>";
    pdg2str[-5112] = "&Sigma;<SUB>b</SUB><SUP>-</SUP>";
    pdg2str[5212] = "&Sigma;<SUB>b</SUB><SUP>0</SUP>";
    pdg2str[-5212] = "&Sigma;<SUB>b</SUB><SUP>0</SUP>";
    pdg2str[5222] = "&Sigma;<SUB>b</SUB><SUP>+</SUP>";
    pdg2str[-5222] = "&Sigma;<SUB>b</SUB><SUP>+</SUP>";
    pdg2str[5114] = "&Sigma;<SUB>b</SUB><SUP>*-</SUP>";
    pdg2str[-5114] = "&Sigma;<SUB>b</SUB><SUP>*-</SUP>";
    pdg2str[5214] = "&Sigma;<SUB>b</SUB><SUP>*0</SUP>";
    pdg2str[-5214] = "&Sigma;<SUB>b</SUB><SUP>*0</SUP>";
    pdg2str[5224] = "&Sigma;<SUB>b</SUB><SUP>*+</SUP>";
    pdg2str[-5224] = "&Sigma;<SUB>b</SUB><SUP>*+</SUP>";
    pdg2str[5132] = "&Xi;<SUB>b</SUB><SUP>-</SUP>";
    pdg2str[-5132] = "&Xi;<SUB>b</SUB><SUP>-</SUP>";
    pdg2str[5232] = "&Xi;<SUB>b</SUB><SUP>0</SUP>";
    pdg2str[-5232] = "&Xi;<SUB>b</SUB><SUP>0</SUP>";
    pdg2str[5312] = "&Xi;<SUB>b</SUB><SUP>\'-</SUP>";
    pdg2str[-5312] = "&Xi;<SUB>b</SUB><SUP>\'-</SUP>";
    pdg2str[5322] = "&Xi;<SUB>b</SUB><SUP>\'0</SUP>";
    pdg2str[-5322] = "&Xi;<SUB>b</SUB><SUP>\'0</SUP>";
    pdg2str[5314] = "&Xi;<SUB>b</SUB><SUP>*-</SUP>";
    pdg2str[-5314] = "&Xi;<SUB>b</SUB><SUP>*-</SUP>";
    pdg2str[5324] = "&Xi;<SUB>b</SUB><SUP>*0</SUP>";
    pdg2str[-5324] = "&Xi;<SUB>b</SUB><SUP>*0</SUP>";
    pdg2str[5332] = "&Omega;<SUB>b</SUB><SUP>-</SUP>";
    pdg2str[-5332] = "&Omega;<SUB>b</SUB><SUP>-</SUP>";
    pdg2str[5334] = "&Omega;<SUB>b</SUB><SUP>*-</SUP>";
    pdg2str[-5334] = "&Omega;<SUB>b</SUB><SUP>*-</SUP>";
    pdg2str[5142] = "&Xi;<SUB>bc</SUB><SUP>0</SUP>";
    pdg2str[-5142] = "&Xi;<SUB>bc</SUB><SUP>0</SUP>";
    pdg2str[5242] = "&Xi;<SUB>bc</SUB><SUP>+</SUP>";
    pdg2str[-5242] = "&Xi;<SUB>bc</SUB><SUP>-</SUP>";
    pdg2str[5412] = "&Xi;<SUB>bc</SUB><SUP>\'0</SUP>";
    pdg2str[-5412] = "&Xi;<SUB>bc</SUB><SUP>\'0</SUP>";
    pdg2str[5422] = "&Xi;<SUB>bc</SUB><SUP>\'+</SUP>";
    pdg2str[-5422] = "&Xi;<SUB>bc</SUB><SUP>\'-</SUP>";
    pdg2str[5414] = "&Xi;<SUB>bc</SUB><SUP>*0</SUP>";
    pdg2str[-5414] = "&Xi;<SUB>bc</SUB><SUP>*0</SUP>";
    pdg2str[5424] = "&Xi;<SUB>bc</SUB><SUP>*+</SUP>";
    pdg2str[-5424] = "&Xi;<SUB>bc</SUB><SUP>*-</SUP>";
    pdg2str[5342] = "&Omega;<SUB>bc</SUB><SUP>0</SUP>";
    pdg2str[-5342] = "&Omega;<SUB>bc</SUB><SUP>0</SUP>";
    pdg2str[5432] = "&Omega;<SUB>bc</SUB><SUP>\'0</SUP>";
    pdg2str[-5432] = "&Omega;<SUB>bc</SUB><SUP>\'0</SUP>";
    pdg2str[5434] = "&Omega;<SUB>bc</SUB><SUP>*0</SUP>";
    pdg2str[-5434] = "&Omega;<SUB>bc</SUB><SUP>*0</SUP>";
    pdg2str[5442] = "&Omega;<SUB>bcc</SUB><SUP>+</SUP>";
    pdg2str[-5442] = "&Omega;<SUB>bcc</SUB><SUP>-</SUP>";
    pdg2str[5444] = "&Omega;<SUB>bcc</SUB><SUP>*+</SUP>";
    pdg2str[-5444] = "&Omega;<SUB>bcc</SUB><SUP>*-</SUP>";
    pdg2str[5512] = "&Xi;<SUB>bb</SUB><SUP>-</SUP>";
    pdg2str[-5512] = "&Xi;<SUB>bb</SUB><SUP>+</SUP>";
    pdg2str[5522] = "&Xi;<SUB>bb</SUB><SUP>0</SUP>";
    pdg2str[-5522] = "&Xi;<SUB>bb</SUB><SUP>0</SUP>";
    pdg2str[5514] = "&Xi;<SUB>bb</SUB><SUP>*-</SUP>";
    pdg2str[-5514] = "&Xi;<SUB>bb</SUB><SUP>*+</SUP>";
    pdg2str[5524] = "&Xi;<SUB>bb</SUB><SUP>*0</SUP>";
    pdg2str[-5524] = "&Xi;<SUB>bb</SUB><SUP>*0</SUP>";
    pdg2str[5532] = "&Omega;<SUB>bb</SUB><SUP>-</SUP>";
    pdg2str[-5532] = "&Omega;<SUB>bb</SUB><SUP>+</SUP>";
    pdg2str[5534] = "&Omega;<SUB>bb</SUB><SUP>*-</SUP>";
    pdg2str[-5534] = "&Omega;<SUB>bb</SUB><SUP>*+</SUP>";
    pdg2str[5542] = "&Omega;<SUB>bbc</SUB><SUP>0</SUP>";
    pdg2str[-5542] = "&Omega;<SUB>bbc</SUB><SUP>0</SUP>";
    pdg2str[5544] = "&Omega;<SUB>bbc</SUB><SUP>*0</SUP>";
    pdg2str[-5544] = "&Omega;<SUB>bbc</SUB><SUP>*0</SUP>";
    pdg2str[5554] = "&Omega;<SUB>bbb</SUB><SUP>-</SUP>";
    pdg2str[-5554] = "&Omega;<SUB>bbb</SUB><SUP>+</SUP>";

    return pdg2str;
}
