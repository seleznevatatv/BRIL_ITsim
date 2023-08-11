// -*- C++ -*-
//
// Package:    BRIL_ITsim/ITClusterPixelExporter
// Class:      ITClusterPixelExporter
//
/**\class ITClusterPixelExporter ITClusterPixelExporter.cc BRIL_ITsim/ITClusterPixelExporter/plugins/ITClusterPixelExporter.cc

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
#include <map>

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
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"

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

struct ITModuleEvent{
    //to be written to TTree, one ITModuleEvent per TTree entry, so one entry in the TTree per module and event
    int event; // event number, one entry
    uint32_t side;  // side, one entry
    uint32_t disk;  // disk, one entry
    uint32_t layer; // front or backside of a single disk, one entry;
    uint32_t ring;  // ring, one entry
    uint32_t module;// module ID, one entry
    std::vector<uint32_t> row;//vector of row coordinates; all digis for this module and event
    std::vector<uint32_t> column;//same as above
    std::vector<uint32_t> adc;//same as above
    //std::vector<std::vector<uint64_t>> clusters;//vector of clusters, each cluster is one entry in the outer vector, each cluster represented by a vector of uint64_t with row in 16 MSB, then col, then ADC, then 0xDEAD in LSB
    std::vector<std::vector<int>> clusters;//vector of clusters, each cluster is one entry in the outer vector, each cluster represented by a vector of int with row in 16 MSB, then col
    std::vector<std::vector<uint32_t>> cluSimTrackId;//simTrackIds for each cluster: a vector of uin32_t for each cluster, each pixel in the clusters has one 32 bit simtrackid
    std::vector<uint32_t> cluMultiplicity;//the multiplicity of each cluster, vector for clusters
    std::vector<std::vector<int>> clusters_charge;

    void clear(){
        this->event = -1;
        this->side = 0;
        this->disk = 0;
        this->layer = 0;
        this->ring = 0;
        this->module = 0;
        this->row.clear();
        this->column.clear();
        this->adc.clear();
        this->clusters.clear();
        this->cluSimTrackId.clear();
        this->cluMultiplicity.clear();
        this->clusters_charge.clear();
    }

    void print()
    {
        std::cout << "Event " << this->event << " Side " << this->side << " Disk " << this->disk << " Layer " << this->layer << " Ring " << this->ring << "Module " << this->module << std::endl;
        //if(this->row.size() == this->column.size() == this->adc.size())
        //{
        std::cout << "Hits: ";
        for (size_t i = 0; i < this->row.size(); i++)
        {
            std::cout << "(" << this->row.at(i) <<","<<this->column.at(i)<<","<<this->adc.at(i)<<"), ";
        }
        std::cout << std::endl << std::endl;
        //}
        //else std::cout << "Error, vectors are not of same size!" <<std::endl;

        int clucount = 0;
        for(auto cluit : this->clusters)
        {
            std::cout << "cluster " << clucount << " : ";
            for(auto hitit : cluit) std::cout << "(" << (hitit>>16) <<","<<(hitit & 0xFFFF) << ") ";
            std::cout << std::endl;
            clucount++;
        }
    }

    void print2()
    {
        int clusterpixelcount = 0;
        for(auto cluster : this->clusters)
            clusterpixelcount += cluster.size();

        std::cout << "Digi size: " << this->row.size() << " | " << this->column.size() << " ---  Cluster size " << clusterpixelcount << std::endl;
    }

    void setModule(int event, unsigned int side, unsigned int disk, unsigned int layer, unsigned int ring, unsigned int module)
    {
        this->event = event;
        this->side = side;
        this->disk = disk;
        this->layer = layer;
        this->ring = ring;
        this->module = module;
    }

    void fillDigis(int row, int column, int adc)
    {
        this->row.push_back(static_cast<uint32_t>(row));
        this->column.push_back(static_cast<uint32_t>(column));
        this->adc.push_back(static_cast<uint32_t>(adc));
    }

    void fillCluster(std::vector<int> tmp_clu)
    {
        this->clusters.push_back(tmp_clu);
    }

    void fillClusterCharge(std::vector<int> tmp_clu)
    {
        this->clusters_charge.push_back(tmp_clu);
    }

    void fillSimLinks(std::vector<uint32_t> SimTrackIds, size_t cluMultiplicity)
    {
        this->cluSimTrackId.push_back(SimTrackIds);
        this->cluMultiplicity.push_back(static_cast<uint32_t>(cluMultiplicity));
    }
};


class ITClusterPixelExporter : public edm::one::EDAnalyzer<edm::one::SharedResources> {
    public:
        explicit ITClusterPixelExporter(const edm::ParameterSet&);
        ~ITClusterPixelExporter();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        // ----------member data ---------------------------
        edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> tgeomHandle;
        edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> tTopoHandle;
        edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster>> m_tokenClusters;
        edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> m_tokenDigis; // digi
        edm::EDGetTokenT<edm::DetSetVector<PixelDigiSimLink>> m_tokenSimLinks;

        // the pointers to geometry, topology and clusters
        // these are members so all functions can access them without passing as argument
        const TrackerTopology* tTopo = NULL;
        const TrackerGeometry* tkGeom = NULL;
        const edmNew::DetSetVector<SiPixelCluster>* clusters = NULL;
        const edmNew::DetSetVector<SiPixelCluster>* clusters_charge = NULL;
        const edm::DetSetVector<PixelDigi>* digis = NULL;  //defining pointer to digis - COB 26.02.19
        const edm::DetSetVector<PixelDigiSimLink>* simlinks = NULL;

        //from config file
        int m_Eventoffset;

        //event counter
        uint32_t m_nevents;

        //internal debug flag
        bool m_debug;
        int modulecounter;
        int toomanycounter;
        int toofewcounter;

        // File service to write ROOT output

        edm::Service<TFileService> m_fileservice;
        ITModuleEvent m_event;
        TTree * m_tree;
};

//
// constants, enums and typedefs
//
std::vector<std::tuple<uint32_t,uint32_t>>   backside_panel_ids({{348655e3, 348680e3},
                        {349175e3, 349205e3},
                        {349703e3, 349723e3},
                        {350227e3, 350248e3},
                        {357042e3, 357063e3},
                        {357567e3, 357588e3},
                        {358091e3, 358112e3},
                        {358616e3, 358637e3}});
//
// static data member definitions
//

//
// constructors and destructor
//
ITClusterPixelExporter::ITClusterPixelExporter(const edm::ParameterSet& iConfig) :
    tgeomHandle(esConsumes())
    , tTopoHandle(esConsumes())
    , m_tokenClusters(consumes<edmNew::DetSetVector<SiPixelCluster>>(iConfig.getParameter<edm::InputTag>("clusters")))
    // , m_tokenClusters(consumes<edmNew::DetSetVector<SiPixelCluster>>(iConfig.getParameter<edm::InputTag>("clusters_charge")))
    , m_tokenDigis(consumes<edm::DetSetVector<PixelDigi>>(iConfig.getParameter<edm::InputTag>("digis"))) //adding digis variable - COB 26.02.19
    , m_tokenSimLinks(consumes<edm::DetSetVector<PixelDigiSimLink>>(iConfig.getParameter<edm::InputTag>("simlinks")))
    , m_Eventoffset(iConfig.getParameter<int>("eventNoOffset"))
{
    usesResource("TFileService");
    m_tree = m_fileservice->make<TTree>("Digis","Digis");
    m_tree->Branch("event",  &m_event.event , "event/i");
    m_tree->Branch("side",  &m_event.side);
    m_tree->Branch("disk",  &m_event.disk);
    m_tree->Branch("layer",  &m_event.layer);
    m_tree->Branch("ring",  &m_event.ring);

    m_tree->Branch("module",  &m_event.module);
    m_tree->Branch("row",    &m_event.row);
    m_tree->Branch("column", &m_event.column);
    m_tree->Branch("adc",    &m_event.adc);
    m_tree->Branch("clusters",   &m_event.clusters);
    //m_tree->Branch("cluSimTrackId",   &m_event.cluSimTrackId);
    m_tree->Branch("cluMultiplicity",   &m_event.cluMultiplicity);
    m_tree->Branch("clusters_charge",   &m_event.clusters_charge);

    //now do what ever initialization is needed
    m_nevents = 0;
    m_debug = true;
    modulecounter = 0;
    toomanycounter = 0;
    toofewcounter = 0;
}

ITClusterPixelExporter::~ITClusterPixelExporter()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void ITClusterPixelExporter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    // Empty event container
    // this is useful here just as a safety mesure but actually I need to do it in the beginning of each module
    this->m_event.clear();

    //get the digis - COB 26.02.19
    edm::Handle<edm::DetSetVector<PixelDigi>> tdigis;
    iEvent.getByToken(m_tokenDigis, tdigis);

    //get the clusters
    edm::Handle<edmNew::DetSetVector<SiPixelCluster>> tclusters;
    iEvent.getByToken(m_tokenClusters, tclusters);

    //get the simlinks
    edm::Handle<edm::DetSetVector<PixelDigiSimLink>> tsimlinks;
    iEvent.getByToken(m_tokenSimLinks, tsimlinks);

    // Get the geometry
    //edm::ESHandle<TrackerGeometry> tgeomHandle;
    //iSetup.get<TrackerDigiGeometryRecord>().get("idealForDigi", tgeomHandle);
    iSetup.getHandle(tgeomHandle);

    // Get the topology
    //edm::ESHandle<TrackerTopology> tTopoHandle;
    //iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);
    iSetup.getHandle(tTopoHandle);

    //get the pointers to geometry and topology
    //tTopo = tTopoHandle.product();
    //const TrackerGeometry* tkGeom = &(*tgeomHandle);
    //tkGeom = tgeomHandle.product();
    tTopo = &iSetup.getData(tTopoHandle);
    tkGeom = &iSetup.getData(tgeomHandle);
    clusters = tclusters.product();
    clusters_charge = tclusters.product();
    digis = tdigis.product();  //pointer to digis - COB 26.02.19
    simlinks = tsimlinks.product();

    int event  = iEvent.id().event();
    event+= m_Eventoffset;


    // Module loop
    for (typename edm::DetSetVector<PixelDigi>::const_iterator DSVit = digis->begin(); DSVit != digis->end(); DSVit++) {

        //get the detid
        unsigned int rawid(DSVit->detId());
        DetId detId(rawid);


        // Determine whether it's barrel or endcap
        TrackerGeometry::ModuleType mType = tkGeom->getDetectorType(detId);
        if (mType == TrackerGeometry::ModuleType::Ph2PXF && detId.subdetId() == PixelSubdetector::PixelEndcap)
        {
            // it's an endcap module
            // so now let's get the exact location

            unsigned int disk = tTopo->pxfDisk(detId);
            if(disk > 8)
            {
                bool thisevent = false;
                //it's TEPX
                unsigned int side = tTopo->pxfSide(detId);
                unsigned int ring = tTopo->pxfBlade(detId);
                unsigned int module = tTopo->pxfModule(detId);
                unsigned int layer = 1;

                for (auto [min_id, max_id] : backside_panel_ids)
                    layer = (detId > min_id && detId < max_id) ? 2 : layer;
                //std::cout << module << std::endl;
                if (ring == 1 && layer ==2) {
                    module += 10;
                } else if (ring == 2 && layer == 2) {
                     module += 14;
                } else if (ring == 3 && layer == 2) {
                     module += 18;
                } else if (ring == 4 && layer == 2) {
                     module += 22;
                } else if (ring == 5 && layer == 2) {
                    module += 24;
                }
                //std::cout << module << std::endl;

                //since every module is it's own entry in the tree I need to clear for every module, otherwise the vector explodes
                this->m_event.clear();
                // now fill the important variables for location
                m_event.setModule(event, side, disk, layer, ring, module);

                //introduce some counters for the size of the digi collection for this module and the total number of pixels in the clusters
                int nDigis=0;
                int nClusterPixels=0;
                int nClusters = 0;
                modulecounter++;

                //if((event == 48 || event == 54) && (side ==2) && (disk == 11) && (ring == 4) && (module == 32))
                //thisevent = true;

                std::map<std::pair<int,int>, uint32_t> adc_data_pixel_map;

                //loop over the digis for this module
                for (edm::DetSet<PixelDigi>::const_iterator digit = DSVit->begin(); digit != DSVit->end(); digit++)
                {
                    // m_event.fillDigis(digit->row(), digit->column(), digit->adc());
                    adc_data_pixel_map[std::make_pair(digit->column(), digit->row())] = digit->adc();
                }
                //append the det set vector iterator size
                // nDigis += DSVit->size();

                //now find the DetSet for SiPixelClusters based on DetId
                edmNew::DetSetVector<SiPixelCluster>::const_iterator theit = clusters->find(detId);
                if (theit != clusters->end())
                {
                    //increment the cluster counter by the size of the DetSetVector size
                    nClusters += theit->size();

                    //now get the simlink detset
                    edm::DetSetVector<PixelDigiSimLink>::const_iterator simLinkDSViter = simlinks->find(detId);
                    if(simLinkDSViter != simlinks->end())
                    {
                        //now iterate the DetSet/Clusters for this DetID
                        for (edmNew::DetSet<SiPixelCluster>::const_iterator cluit = theit->begin(); cluit != theit->end(); cluit++)
                        {
                            //temporary vector to hold the individual pixels of the cluster
                            //std::vector<uint64_t> tmpClu;
                            std::vector<int> tmpClu;
                            //temporary vector to hold the SimTrackIds of the pixels in the cluster
                            std::vector<uint32_t> tmpSimTrackIds;
                            //temporary std::set to count the number of unique SimTrackIds;
                            std::set<unsigned int> simTrackIds;

                            //for each cluster, get the size and push the x and y into a vector of pair
                            int size = cluit->size();
                            nClusterPixels += size;
                            //if(size ==1) std::cout << "Single Pixel Cluster" << std::endl;

                            for (int i = 0; i < size; i++) {
                                SiPixelCluster::Pixel pix = cluit->pixel(i);
                                uint16_t x = pix.x;
                                uint16_t y = pix.y;
                                uint64_t adc = pix.adc;

                                auto adc_val = adc_data_pixel_map.find(std::make_pair(x,y));

                                if(adc_val != adc_data_pixel_map.end()) {
                                    m_event.fillDigis((int)y, (int)x, adc_val->second);
                                
                                    std::cout << (int)x << " " << (int)y << " " << adc_val->second << std::endl;
                                }
                                //uint64_t dead = 0xDEAD;>
                                //int tmp = x << 48 | y << 32 | adc << 16 | dead;
                                int tmp = x << 16 | y;//
                                tmpClu.push_back(tmp);

                                //now deal with the simLinks
                                unsigned int clusterChannel = PixelDigi::pixelToChannel(pix.x, pix.y);

                                //loop over all SimLinks for this detId and check pixel by pixel wheather there is a sim link
                                bool found = false;
                                for (edm::DetSet<PixelDigiSimLink>::const_iterator it = simLinkDSViter->data.begin(); it != simLinkDSViter->data.end(); it++)
                                {
                                    if (clusterChannel == it->channel())
                                    {
                                        found = true;
                                        //indeed, this channel has a simLink
                                        simTrackIds.insert(it->SimTrackId());
                                        tmpSimTrackIds.push_back(it->SimTrackId());
                                    }
                                }//end of search loop
                                if(thisevent) std::cout << std::endl;
                                if(found == false) std::cout << "warning, nothing found in SimLink collection for pixel " << x << " " << y << " channel " << clusterChannel << std::endl;
                            }//end of loop over all pixels in cluster

                            //now push back the vector with the pixels in the cluster to the cluster vector
                            m_event.fillCluster(tmpClu);
                            m_event.fillClusterCharge(tmpClu);

                            //now push back the SimTrackIds for this cluster
                            m_event.fillSimLinks(tmpSimTrackIds, simTrackIds.size());

                            //quick sanity check
                            //if(tmpClu.size() != tmpSimTrackIds.size()) std::cout << "Warning - clustersize: "<< tmpClu.size() << " sim link vector size " << tmpSimTrackIds.size() << " # of unique SimTrackIDs " << simTrackIds.size() << std::endl;
                            //else std::cout << "Same size!" << std::endl;
                        }//end of cluster iteration



                    }//end of detID found in Sim links
                    else
                    {
                        std::cout << "No simlinks for this DetID" << std::endl;
                        continue;
                    }
                }//end of no clusters for this DetID condition
                else
                {
                    std::cout << "No clusters for this DetID" << std::endl;
                    continue;
                }

                nDigis += this->m_event.adc.size();

                if(thisevent)m_event.print();
                this->m_tree->Fill();

                //print some debug information
                if(nDigis != nClusterPixels)
                {
                    std::cout << "DEBUG: Evnet "<<event << " Side " << side << " Disk " << disk << "layer" << layer << " Ring " << ring << " Module " << module << " --- Digis: "<<nDigis << " | Pixels in Clusters: " << nClusterPixels << " --- nClusters: " << nClusters << std::endl;
                    this->m_event.print2();
                    if(nClusterPixels > nDigis)
                    {
                        toomanycounter++;
                    }
                    else
                        toofewcounter++;

                }

            }//end of disk > 8 condition
        }//end of is endcap condition
    }//end of Module loop

    m_nevents++;
}

// ------------ method called once each job just before starting event loop  ------------
void ITClusterPixelExporter::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void ITClusterPixelExporter::endJob()
{
    //this->m_tree->Write();
    std::cout << "Processed "<< m_nevents << " events with " << modulecounter << " modules (total) of which " << toomanycounter << " having more pixels in clusters than in the digis and " << toofewcounter << " having more digis than pixels in clusters" << std::endl;

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ITClusterPixelExporter::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
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
DEFINE_FWK_MODULE(ITClusterPixelExporter);
