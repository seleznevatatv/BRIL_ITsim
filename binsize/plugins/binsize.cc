// -*- C++ -*-
//
// Package:    BRIL_ITsim/binsize
// Class:      binsize
//
/**\class binsize binsize.cc BRIL_ITsim/binsize/plugins/binsize.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Georg Auzinger
//         Created:  Tue, 31 Mar 2020 13:31:06 GMT
//
//


// system include files
#include <algorithm>
#include <memory>
#include <assert.h>
#include <math.h>

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
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
//#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <TH2F.h>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.



class binsize : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
    public:
        explicit binsize(const edm::ParameterSet&);
        ~binsize();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;
        int getQuarterRing(DetId theID);
        int getBin(int ring, int quarter);

        // ----------member data ---------------------------
        edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster>> m_tokenClusters;
        //edm::EDGetTokenT<edm::DetSetVector<PixelDigiSimLink>> m_tokenSimLinks;
        edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> m_tokenDigis;

        // the pointers to geometry, topology and clusters
        // these are members so all functions can access them without passing as argument
        const TrackerTopology* tTopo = NULL;
        const TrackerGeometry* tkGeom = NULL;
        const edmNew::DetSetVector<SiPixelCluster>* clusters = NULL;
        //const edm::DetSetVector<PixelDigiSimLink>* simlinks = NULL;
        const edm::DetSetVector<PixelDigi>* digis = NULL;  //defining pointer to digis - COB 26.02.19

        //max bins of Counting histogram
        uint32_t m_maxBinTEPX;
        uint32_t m_maxBinD4R1;
        //event counter
        uint32_t m_nevents;
        // bunches per orbit
        uint32_t m_bxperorbit;
        // trigger frequency for TEPX and D4R1
        // in kHz
        uint32_t m_trgfrqTEPX;
        uint32_t m_trgfrqD4R1;
        // number of events to integrate for both parts of the detector
        uint32_t m_nEventsTEPX;
        uint32_t m_nEventsD4R1;
        //counter to hold the integrated number of clusters per ring, resetting every N events
        //disk/ring/quarter
        uint32_t TEPXcounter[8][5][4];
        //same for D4R1, resetting every M events
        //side/quarter
        uint32_t D4R1counter[2][4];
        //Histogram
        TH2F* m_diskHistosTEPX[8];
        TH2F* m_diskHistosD4R1[2];
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
binsize::binsize(const edm::ParameterSet& iConfig)
    : //m_tokenClusters(consumes<edmNew::DetSetVector<SiPixelCluster>> ("clusters"))
        m_tokenClusters(consumes<edmNew::DetSetVector<SiPixelCluster>>(iConfig.getParameter<edm::InputTag>("clusters")))
        //, m_tokenSimLinks(consumes<edm::DetSetVector<PixelDigiSimLink>>(iConfig.getParameter<edm::InputTag>("simlinks")))
        //, m_tokenDigis(consumes<edm::DetSetVector<PixelDigi>>(iConfig.getParameter<edm::InputTag>("digis"))) //adding digis variable - COB 26.02.19
        , m_maxBinTEPX(iConfig.getUntrackedParameter<uint32_t>("maxBinTEPX"))
        , m_maxBinD4R1(iConfig.getUntrackedParameter<uint32_t>("maxBinD4R1"))
        , m_bxperorbit(iConfig.getUntrackedParameter<uint32_t>("bxperorbit"))
        , m_trgfrqTEPX(iConfig.getUntrackedParameter<uint32_t>("trgfrqTEPX"))
        , m_trgfrqD4R1(iConfig.getUntrackedParameter<uint32_t>("trgfrqD4R1"))

{
    //now do what ever initialization is needed
    m_nevents = 0;
    std::memset(TEPXcounter,0, sizeof(TEPXcounter));
    std::memset(D4R1counter,0, sizeof(D4R1counter));

    //now calculate the number of events to integrate for the two parts of the detector
    //assume the colliding BX per orbit from the config and orbit frequency
    float flhc = 11.2233;
    float eventRate = flhc * (float)m_bxperorbit;
    //this is the number of orbits needed to sample every collding BX
    //calculated as eventRate/trigger rate
    //alternatively colliding_bx / (trigger rate/orbit frequency)
    float orbits_to_sample_all_TEPX = ceil(eventRate/m_trgfrqTEPX);
    float orbits_to_sample_all_D4R1 = ceil(eventRate/m_trgfrqD4R1);
    float orbits_in_NB4 = 16384.; //2^14
    m_nEventsTEPX = (uint32_t)floor(orbits_in_NB4/orbits_to_sample_all_TEPX);
    m_nEventsD4R1 = (uint32_t)floor(orbits_in_NB4/orbits_to_sample_all_D4R1);

    std::cout << "**************************************************************************************************************************" << std::endl;
    std::cout << "SUMMARY OF SIMULATION PARAMETERS" << std::endl;
    std::cout << "Assuming " << m_bxperorbit << " colliding bunches per Orbit and fLHC = " << flhc << " kHz" << std::endl;
    std::cout << "Assuming fTRG_TEPX = " << m_trgfrqTEPX << " kHz and fTRG_D4R1 = " << m_trgfrqD4R1 << " kHz" << std::endl;
    std::cout << "**************************************************************************************************************************" << std::endl;
    std::cout << "Resulting in an Event rate of " << eventRate/1000. << " MHz " << std::endl;
    std::cout << "Number of Orbits to sample every colliding bunch in the Orbit( " << m_bxperorbit << " )" << std::endl;
    std::cout << "TEPX: " << orbits_to_sample_all_TEPX << " | D4R1: " << orbits_to_sample_all_D4R1 << std::endl;
    std::cout << "With 16384 orbits in a NB4 this allows to sample every collding bunch N times: " << std::endl;
    std::cout << "TEPX: " << m_nEventsTEPX << " | D4R1: " << m_nEventsD4R1 << std::endl;
    std::cout << "**************************************************************************************************************************" << std::endl;


}


binsize::~binsize()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
//helper to quickly return quarter Ring as function of module location
//returns -1 is wrong disk
int binsize::getQuarterRing(DetId theID)
{
    //module type => need phase 2 pixel forward module, in endcap
    TrackerGeometry::ModuleType mType = tkGeom->getDetectorType(theID);
    if (mType == TrackerGeometry::ModuleType::Ph2PXF && theID.subdetId() == PixelSubdetector::PixelEndcap)
    {
        //obtaining side, disk, ring of module
        //unsigned int side = (tTopo->pxfSide(detId));
        unsigned int disk = (tTopo->pxfDisk(theID));
        unsigned int ring = (tTopo->pxfBlade(theID));
        unsigned int module = (tTopo->pxfModule (theID));

        if (disk > 8) { // TEPX modules
            switch(ring) {
                case 1:
                    //20 modules in TEPX Ring 1
                    if(module <=5) return 0;
                    else if (module >5 && module <=10) return 1;
                    else if (module > 10 && module <=15) return 2;
                    else return 3;
                    break;
                case 2:
                    //28 modules in TEPX Ring 2
                    if(module <=7) return 0;
                    else if (module >7 && module <=14) return 1;
                    else if (module > 14 && module <=21) return 2;
                    else return 3;
                    break;
                case 3:
                    //36 modules in TEPX Ring 2
                    if(module <=9) return 0;
                    else if (module >9 && module <=18) return 1;
                    else if (module > 18 && module <=27) return 2;
                    else return 3;
                    break;
                case 4:
                    //44 modules in TEPX Ring 2
                    if(module <=11) return 0;
                    else if (module >11 && module <=22) return 1;
                    else if (module > 22 && module <=33) return 2;
                    else return 3;
                    break;
                case 5:
                    //48 modules in TEPX Ring 2
                    if(module <=12) return 0;
                    else if (module >12 && module <=24) return 1;
                    else if (module > 24 && module <=36) return 2;
                    else return 3;
                    break;
            }
        }
        else return -1;
    }
    return -1;
}

//helper to get histogram bin for a given ring+quarter combination
int binsize::getBin(int ring, int quarter){
    //ring goes from 0 to 4
    //quarter goes from 0 to 3
    //bins go from 1 to 20
    //if(ring>0&& ring <=5) {
    //ring-=1;
    return 4*ring+quarter+1;
    //}
    //else
    //return -1;
}

// ------------ method called for each event  ------------
    void
binsize::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    //get the digis - COB 26.02.19
    //edm::Handle<edm::DetSetVector<PixelDigi>> tdigis;
    //iEvent.getByToken(m_tokenDigis, tdigis);

    //get the clusters
    edm::Handle<edmNew::DetSetVector<SiPixelCluster>> tclusters;
    iEvent.getByToken(m_tokenClusters, tclusters);

    //get the simlinks
    //edm::Handle<edm::DetSetVector<PixelDigiSimLink>> tsimlinks;
    //iEvent.getByToken(m_tokenSimLinks, tsimlinks);

    // Get the geometry
    edm::ESHandle<TrackerGeometry> tgeomHandle;
    iSetup.get<TrackerDigiGeometryRecord>().get("idealForDigi", tgeomHandle);

    // Get the topology
    edm::ESHandle<TrackerTopology> tTopoHandle;
    iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);

    //get the pointers to geometry, topology and clusters
    tTopo = tTopoHandle.product();
    //const TrackerGeometry* tkGeom = &(*geomHandle);
    tkGeom = tgeomHandle.product();
    clusters = tclusters.product();
    //simlinks = tsimlinks.product();
    //digis = tdigis.product();  //pointer to digis - COB 26.02.19

    //change accordingly for hits or clusters
    //for (typename edm::DetSetVector<PixelDigi>::const_iterator DSVit = digis->begin(); DSVit != digis->end(); DSVit++) {
    for (typename edmNew::DetSetVector<SiPixelCluster>::const_iterator DSVit = clusters->begin(); DSVit != clusters->end(); DSVit++) {

        //get the detid
        unsigned int rawid(DSVit->detId());
        DetId detId(rawid);

        //debugging...
        //std::cout << "DetId " << std::hex << "0x" << rawid << std::dec << " " << detId.det() << " " << detId.subdetId() << " "
        //          << ((rawid >> 23) & 0x3) << " " << ((rawid >> 18) & 0xF) << " " << ((rawid >> 12) & 0x3F) << " "
        //          << ((rawid >> 2) & 0xFF) << std::endl;
        //std::cout << "side " << tTopo->pxfSide(detId) << " layer " << tTopo->pxfDisk(detId) << " ring " << tTopo->pxfBlade(detId)
        //          << " module " << tTopo->pxfModule(detId) << std::endl;

        //module type => need phase 2 pixel forward module, in endcap
        TrackerGeometry::ModuleType mType = tkGeom->getDetectorType(detId);
        if (mType != TrackerGeometry::ModuleType::Ph2PXF && detId.subdetId() != PixelSubdetector::PixelEndcap)
            continue;

        //obtaining side, disk, ring of module
        unsigned int side = (tTopo->pxfSide(detId));
        unsigned int disk = (tTopo->pxfDisk(detId));
        unsigned int ring = (tTopo->pxfBlade(detId));
        //unsigned int module = (tTopo->pxfModule (detId));
        int quarterRing =  this->getQuarterRing(detId);

        if (disk > 8) { // TEPX modules

            //index in histogram map
            int hist_id = -1;
            unsigned int ring_id = ring - 1;
            if (side == 1) {
                //this is a TEPX- hit on side1
                hist_id = disk - 9;
            } else if (side == 2) {
                //this is a TEPX+ hit on side 2
                hist_id = 4 + disk - 9;
            }

            //std::cout << "DEBUG " << hist_id << " " << ring_id << " " << quarterRing << std::endl;
            TEPXcounter[hist_id][ring_id][quarterRing]+=DSVit->size();
            if(disk ==12 && ring ==1)
                D4R1counter[side-1][quarterRing]+=DSVit->size();

        }//end of disk>8 = TEPX
    }//end of DSVIterator = module loop
    m_nevents++;

    //now if the integration interval is reached, fill the histogramm with the counter value and reset the counter
    if(m_nevents % m_nEventsTEPX == 0)
    {
        //TODO: here fill the histogram with the counter value for the respective bin and then reset the counter array
        for(unsigned disk = 0; disk < 8; disk++)
        {
            for(unsigned ring = 0; ring < 5; ring++)
            {
                for(unsigned quarter = 0; quarter  < 4; quarter++)
                {
                    m_diskHistosTEPX[disk]->Fill(this->getBin(ring, quarter), TEPXcounter[disk][ring][quarter]);
                    //std::cout << "Disk: " << disk << " Ring: " << ring << " Quarter: " << quarter << " counts in " << m_nEventsTEPX << " = " << TEPXcounter[disk][ring][quarter] << std::endl;
                }
            }

        }
        std::cout << "TEPX: Integrated " << m_nEventsTEPX << " Events - filling Histogram and resetting counters" << std::endl;
        std::memset(TEPXcounter, 0, sizeof(TEPXcounter));
    }
    if(m_nevents % m_nEventsD4R1 == 0)
    {
        //TODO: here fill the histogram with the counter value for the respective bin and then reset the counter array
        for(unsigned side = 0; side < 2; side++)
        {
            for(unsigned quarter = 0; quarter < 4; quarter++)
            {
                m_diskHistosD4R1[side]->Fill(quarter+1, D4R1counter[side][quarter]);
                //std::cout << "Side: " << side << " Quarter: " << quarter << " counts in " << m_nEventsD4R1 << " = " << D4R1counter[side][quarter] << std::endl;
            }
        }
        std::cout << "D4R1: Integrated " << m_nEventsD4R1 << " Events - filling Histogram and resetting counters" << std::endl;
        std::memset(D4R1counter, 0, sizeof(D4R1counter));
    }
}


// ------------ method called once each job just before starting event loop  ------------
    void
binsize::beginJob()
{
    edm::Service<TFileService> fs;

    fs->file().cd("/");
    TFileDirectory td = fs->mkdir("TEPX");

    //now lets create the histograms
    for (unsigned int i = 0; i < 8; i++) {
        int disk = (i < 4) ? i - 4 : i - 3;
        std::stringstream histoname;
        histoname << "Number of clusters for quarter Rings Disk " << disk << ";Ring+Quarter;# of Clusters per NB4";
        std::stringstream histotitle;
        histotitle << "Number of clusters for Disk " << disk;
        //name, name, nbinX, Xlow, Xhigh, nbinY, Ylow, Yhigh
        m_diskHistosTEPX[i] = td.make<TH2F>(histotitle.str().c_str(), histoname.str().c_str(), 20, .5, 20.5, m_maxBinTEPX/100, 0, m_maxBinTEPX);
        m_diskHistosTEPX[i]->GetXaxis()->SetBinLabel(1,"Ring 1, Quarter 1");
        m_diskHistosTEPX[i]->GetXaxis()->SetBinLabel(2,"Ring 1, Quarter 2");
        m_diskHistosTEPX[i]->GetXaxis()->SetBinLabel(3,"Ring 1, Quarter 3");
        m_diskHistosTEPX[i]->GetXaxis()->SetBinLabel(4,"Ring 1, Quarter 4");
        m_diskHistosTEPX[i]->GetXaxis()->SetBinLabel(5,"Ring 2, Quarter 1");
        m_diskHistosTEPX[i]->GetXaxis()->SetBinLabel(6,"Ring 2, Quarter 2");
        m_diskHistosTEPX[i]->GetXaxis()->SetBinLabel(7,"Ring 2, Quarter 3");
        m_diskHistosTEPX[i]->GetXaxis()->SetBinLabel(8,"Ring 2, Quarter 4");
        m_diskHistosTEPX[i]->GetXaxis()->SetBinLabel(9,"Ring 3, Quarter 1");
        m_diskHistosTEPX[i]->GetXaxis()->SetBinLabel(10,"Ring 3, Quarter 2");
        m_diskHistosTEPX[i]->GetXaxis()->SetBinLabel(11,"Ring 3, Quarter 3");
        m_diskHistosTEPX[i]->GetXaxis()->SetBinLabel(12,"Ring 3, Quarter 4");
        m_diskHistosTEPX[i]->GetXaxis()->SetBinLabel(13,"Ring 4, Quarter 1");
        m_diskHistosTEPX[i]->GetXaxis()->SetBinLabel(14,"Ring 4, Quarter 2");
        m_diskHistosTEPX[i]->GetXaxis()->SetBinLabel(15,"Ring 4, Quarter 3");
        m_diskHistosTEPX[i]->GetXaxis()->SetBinLabel(16,"Ring 4, Quarter 4");
        m_diskHistosTEPX[i]->GetXaxis()->SetBinLabel(17,"Ring 5, Quarter 1");
        m_diskHistosTEPX[i]->GetXaxis()->SetBinLabel(18,"Ring 5, Quarter 2");
        m_diskHistosTEPX[i]->GetXaxis()->SetBinLabel(19,"Ring 5, Quarter 3");
        m_diskHistosTEPX[i]->GetXaxis()->SetBinLabel(20,"Ring 5, Quarter 4");
    }
    m_diskHistosD4R1[0] = td.make<TH2F>("Number of Clusters per Quarter Ring D4R1-z","Number of Clusters per Quarter Ring D4R1 -z; Quarter; # of Clusters per NB4",4, .5,4.5, m_maxBinD4R1/100, 0, m_maxBinD4R1);
    m_diskHistosD4R1[1] = td.make<TH2F>("Number of Clusters per Quarter Ring D4R1+z","Number of Clusters per Quarter Ring D4R1 +z; Quarter; # of Clusters per NB4",4, .5,4.5, m_maxBinD4R1/100, 0, m_maxBinD4R1);
}

// ------------ method called once each job just after ending the event loop  ------------
    void
binsize::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
binsize::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(binsize);
