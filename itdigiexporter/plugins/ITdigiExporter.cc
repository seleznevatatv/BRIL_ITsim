// -*- C++ -*-
//
// Package:    BRIL_ITsim/ITdigiExporter
// Class:      ITdigiExporter
//
/**\class ITdigiExporter ITdigiExporter.cc BRIL_ITsim/ITdigiExporter/plugins/ITdigiExporter.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Georg Auzinger
//         Created:  Thu, 17 Jan 2019 13:12:52 GMT
// Edits: Alexander Ruede
//
//

// system include files
// system include files
#include <fstream>
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
//#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
//#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h" // digi

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <TTree.h>
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

struct ITDigiEvent {
    // Struct to store address and ADC information about all Digis in an event

    long int event;
    std::vector<bool> barrel;               // Whether the Digi is from TBPX (true) or TFPX/TEPX (false)
    std::vector<unsigned int> diskladder;   // The ladder (barrel) or disk number (forward)
    std::vector<unsigned int> ringlayer;    // The layer (barrel) or ring number (forward)
    std::vector<unsigned int> module;       // The module number
    std::vector<unsigned int> row;          // The pixel row number on the module
    std::vector<unsigned int> column;       // The pixel column number on the module
    std::vector<unsigned int> adc;          // ADC value for given pixel
    std::vector<std::vector<int>> clusters;          // ADC value for given pixel 16MSBs: x, 16 LSBs: y

    void clear(){
        this->event = -1;
        this->barrel.clear();
        this->diskladder.clear();
        this->ringlayer.clear();
        this->module.clear();
        this->row.clear();
        this->column.clear();
        this->adc.clear();
        this->clusters.clear();
    }
};

class ITdigiExporter : public edm::one::EDAnalyzer<edm::one::SharedResources> {
    public:
        explicit ITdigiExporter(const edm::ParameterSet&);
        ~ITdigiExporter();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        // ----------member data ---------------------------
        edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster>> m_tokenClusters;
        edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> m_tokenDigis; // digi

        // the pointers to geometry, topology and clusters
        // these are members so all functions can access them without passing as argument
        const TrackerTopology* tTopo = NULL;
        const TrackerGeometry* tkGeom = NULL;
        const edmNew::DetSetVector<SiPixelCluster>* clusters = NULL;
        const edm::DetSetVector<PixelDigi>* digis = NULL;  //defining pointer to digis - COB 26.02.19

        //from config file
        uint32_t m_disk;
        uint32_t m_ring;
        bool m_hittype;
        // uint32_t m_module;

        //event counter
        uint32_t m_nevents;

        // File service to write ROOT output

        edm::Service<TFileService> m_fileservice;
        ITDigiEvent m_event;
        TTree * m_tree;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ITdigiExporter::ITdigiExporter(const edm::ParameterSet& iConfig) :
    m_tokenClusters(consumes<edmNew::DetSetVector<SiPixelCluster>>(iConfig.getParameter<edm::InputTag>("clusters")))
    , m_tokenDigis(consumes<edm::DetSetVector<PixelDigi>>(iConfig.getParameter<edm::InputTag>("digis"))) //adding digis variable - COB 26.02.19
    , m_disk(iConfig.getUntrackedParameter<uint32_t>("disk"))
    , m_ring(iConfig.getUntrackedParameter<uint32_t>("ring"))
    , m_hittype(iConfig.getUntrackedParameter<bool>("hits"))
{
    // TTree version 0: Save everything into one TTree
    // TODO: Organize data from different rings, detector parts
    // either into same TTree or different TTrees
    usesResource("TFileService");
    m_tree = m_fileservice->make<TTree>("Digis","Digis");
    m_tree->Branch("event",  &m_event.event , "event/i");
    m_tree->Branch("barrel",  &m_event.barrel);
    m_tree->Branch("diskladder",   &m_event.diskladder);
    m_tree->Branch("ringlayer",   &m_event.ringlayer);
    m_tree->Branch("module", &m_event.module);
    if(m_hittype)
    {
        m_tree->Branch("row",    &m_event.row);
        m_tree->Branch("column", &m_event.column);
        m_tree->Branch("adc",    &m_event.adc);
    }
    else
        m_tree->Branch("clusters",   &m_event.clusters);

    //now do what ever initialization is needed
    m_nevents = 0;
}

ITdigiExporter::~ITdigiExporter()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void ITdigiExporter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    // Empty event container
    this->m_event.clear();
    this->m_event.event = iEvent.id().event();

    //get the digis - COB 26.02.19
    edm::Handle<edm::DetSetVector<PixelDigi>> tdigis;
    iEvent.getByToken(m_tokenDigis, tdigis);

    //get the clusters
    edm::Handle<edmNew::DetSetVector<SiPixelCluster>> tclusters;
    iEvent.getByToken(m_tokenClusters, tclusters);

    // Get the geometry
    edm::ESHandle<TrackerGeometry> tgeomHandle;
    iSetup.get<TrackerDigiGeometryRecord>().get("idealForDigi", tgeomHandle);

    // Get the topology
    edm::ESHandle<TrackerTopology> tTopoHandle;
    iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);

    //get the pointers to geometry and topology
    tTopo = tTopoHandle.product();
    //const TrackerGeometry* tkGeom = &(*tgeomHandle);
    tkGeom = tgeomHandle.product();
    clusters = tclusters.product();
    digis = tdigis.product();  //pointer to digis - COB 26.02.19

    // Module loop
    for (typename edm::DetSetVector<PixelDigi>::const_iterator DSVit = digis->begin(); DSVit != digis->end(); DSVit++) {

        //get the detid
        unsigned int rawid(DSVit->detId());
        DetId detId(rawid);


        // Determine whether it's barrel or endcap
        TrackerGeometry::ModuleType mType = tkGeom->getDetectorType(detId);
        bool barrel;
        if (mType == TrackerGeometry::ModuleType::Ph2PXF && detId.subdetId() == PixelSubdetector::PixelEndcap) {
            barrel = false;
        } else if (mType == TrackerGeometry::ModuleType::Ph2PXB  && detId.subdetId() == PixelSubdetector::PixelBarrel) {
            barrel = true;
            continue;

        }// else {
        //continue;
        //}

        //obtaining location of module
        unsigned int side, diskladder, ringlayer, module;
        diskladder = tTopo->pxfDisk(detId);
        ringlayer = tTopo->pxfBlade(detId);
        module = tTopo->pxfModule(detId);
        side = tTopo->pxfSide(detId);

        // Loop over the digis in each module
        if(side ==2)
        {
            //only consider +z
            if(m_hittype)
            {
                for (edm::DetSet<PixelDigi>::const_iterator digit = DSVit->begin(); digit != DSVit->end(); digit++) {
                    this->m_event.barrel.push_back(barrel);
                    this->m_event.ringlayer.push_back(ringlayer);
                    this->m_event.diskladder.push_back(diskladder);
                    this->m_event.module.push_back(module);
                    this->m_event.row.push_back(digit->row());
                    this->m_event.column.push_back(digit->column());
                    this->m_event.adc.push_back(digit->adc());
                }
            }
            else{
                //now find the DetSet for SiPixelClusters based on DetId
                edmNew::DetSetVector<SiPixelCluster>::const_iterator theit = clusters->find(detId);
                if (theit != clusters->end()) {
                    //now iterate the DetSet/Clusters for this DetID
                    for (edmNew::DetSet<SiPixelCluster>::const_iterator cluit = theit->begin(); cluit != theit->end(); cluit++) {

                        std::vector<int> tmpClu;
                        //for each cluster, get the size and push the x and y into a vector of pair
                        int size = cluit->size();
                        for (int i = 0; i < size; i++) {
                            SiPixelCluster::Pixel pix = cluit->pixel(i);
                            uint16_t x = pix.x;
                            uint16_t y = pix.y;
                            int tmp = x << 16 | y;
                            tmpClu.push_back(tmp);
                        }
                        //now push back the vector of pair in the clusters vector
                        this->m_event.barrel.push_back(barrel);
                        this->m_event.ringlayer.push_back(ringlayer);
                        this->m_event.diskladder.push_back(diskladder);
                        this->m_event.module.push_back(module);
                        this->m_event.clusters.push_back(tmpClu);
                    }
                }
                else
                    continue;

            }
        }
        else continue;
    }
    this->m_tree->Fill();
    m_nevents++;
}

// ------------ method called once each job just before starting event loop  ------------
void ITdigiExporter::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void ITdigiExporter::endJob()
{

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ITdigiExporter::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);

    //Specify that only 'tracks' is allowed
    //To use, remove the default given above and uncomment below
    //ParameterSetDescription desc;
    //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
    //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ITdigiExporter);
