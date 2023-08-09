// -*- C++ -*-
//
// Package:    BRIL_ITsim/ITdigiExporter
// Class:      ITdigiExporter
//
/**\class ITdigiExporter ITdigiExporter.cc
BRIL_ITsim/ITdigiExporter/plugins/ITdigiExporter.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Georg Auzinger
//         Created:  Thu, 17 Jan 2032 13:12:52 GMT
// Edits: Max Bensink
//
//

// #define DEBUG

// system include files
// system include files
#include <fstream>
#include <memory>
#include <algorithm>
#include <ranges>

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/ESInputTag.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CondFormats/DataRecord/interface/TrackerDetToDTCELinkCablingMapRcd.h"
#include "CondFormats/SiPhase2TrackerObjects/interface/DTCELinkId.h"
#include "CondFormats/SiPhase2TrackerObjects/interface/TrackerDetToDTCELinkCablingMap.h"

#include <TTree.h>

#include "ClusterModel.h"
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.
constexpr int32_t row_max = 1354;
constexpr int32_t column_max = 434;
constexpr int32_t row_cutoff = 10;
constexpr int32_t column_cutoff = 2;
constexpr int32_t actual_row_max = row_max - row_cutoff;
constexpr int32_t actual_column_max = column_max - column_cutoff;

constexpr int32_t second_col_chip_start = actual_column_max / 2;
constexpr int32_t second_row_chip_start = actual_row_max / 2;

typedef std::pair<int, int> vertex_t;

struct ITDigiModule
{
  // side of the collider (1 == left and 2 == right)
  uint8_t side;
  // disk per side (max 4)
  uint8_t disk;
  // module per ring (max 25)
  uint8_t module;
  // layer per disk (1 == front and 2 == back)
  uint8_t layer;
  // ring per disk (max 5)
  uint8_t ring;
  // number of clusters per module
  uint32_t hlt_nclusters;
  // number of clusters per module
  uint32_t bril_nclusters;
  // connected elink
  uint8_t elink;
  uint8_t lpgbt;
  uint16_t dtc;

  // id of module
  uint64_t detid;

  std::vector<unsigned short> row;    // The pixel row number on the module
  std::vector<unsigned short> column; // The pixel column number on the module
  std::vector<unsigned char> adc;     // ADC value for given pixel

  std::vector<unsigned char> bril_hit_stage; // through which stage did this pixel go

  std::vector<unsigned short> bril_cluster_row;
  std::vector<unsigned short> bril_cluster_column;
  std::vector<unsigned short> bril_cluster_size_x;
  std::vector<unsigned short> bril_cluster_size_y;
  std::vector<unsigned short> bril_cluster_nbits;

  std::vector<unsigned short> hlt_nbits;
  std::vector<unsigned short> hlt_size_x;
  std::vector<unsigned short> hlt_size_y;
  std::vector<unsigned short> hlt_cluster_row;
  std::vector<unsigned short> hlt_cluster_column;

  uint32_t col_not_adjacent;
  uint32_t col_already_used;
  uint32_t col_not_found;
  uint32_t col_rhs_must_be_bigger;
  uint32_t row_processor_buffer_not_empty;
  uint32_t row_already_used;
  uint32_t row_not_found;
};

std::vector<std::tuple<uint32_t, uint32_t>>
    backside_panel_ids({{348655e3, 348680e3},
                        {349175e3, 349205e3},
                        {349703e3, 349723e3},
                        {350227e3, 350248e3},
                        {357042e3, 357063e3},
                        {357567e3, 357588e3},
                        {358091e3, 358112e3},
                        {358616e3, 358637e3}});

struct ITDigiEvent
{
  // Struct to store address and ADC information about all Digis in an event

  long int event;
  unsigned int nmodule;

  TTree *tree = NULL;
  std::vector<ITDigiModule> modules; // The module number

  std::vector<unsigned int> branch_detid;
  std::vector<unsigned char> branch_module;
  std::vector<unsigned char> branch_disk;
  std::vector<unsigned char> branch_layer;
  std::vector<unsigned char> branch_ring;

  std::vector<uint8_t> branch_side;
  std::vector<uint8_t> branch_elink;
  std::vector<uint8_t> branch_lpgbt;
  std::vector<uint16_t> branch_dtc;

  std::vector<std::vector<unsigned short>> branch_row;    // The pixel row number on the module
  std::vector<std::vector<unsigned short>> branch_column; // The pixel column number on the module
  std::vector<std::vector<unsigned char>> branch_adc;     // ADC value for given pixel

  std::vector<std::vector<unsigned char>> branch_bril_hit_stage; // through which stage did this pixel go

  std::vector<unsigned int> branch_bril_nclusters;
  std::vector<std::vector<unsigned short>> branch_bril_cluster_row;
  std::vector<std::vector<unsigned short>> branch_bril_cluster_column;
  std::vector<std::vector<unsigned short>> branch_bril_cluster_size_x;
  std::vector<std::vector<unsigned short>> branch_bril_cluster_size_y;
  std::vector<std::vector<unsigned short>> branch_bril_cluster_nbits;

  std::vector<unsigned int> branch_hlt_nclusters;
  std::vector<std::vector<unsigned short>> branch_hlt_cluster_row;
  std::vector<std::vector<unsigned short>> branch_hlt_cluster_column;
  std::vector<std::vector<unsigned short>> branch_hlt_cluster_size_x;
  std::vector<std::vector<unsigned short>> branch_hlt_cluster_size_y;
  std::vector<std::vector<unsigned short>> branch_hlt_cluster_nbits;

  std::vector<unsigned int> branch_col_not_adjacent;
  std::vector<unsigned int> branch_col_already_used;
  std::vector<unsigned int> branch_col_not_found;
  std::vector<unsigned int> branch_col_rhs_must_be_bigger;
  std::vector<unsigned int> branch_row_processor_buffer_not_empty;
  std::vector<unsigned int> branch_row_already_used;
  std::vector<unsigned int> branch_row_not_found;

  void clear()
  {
    this->event = -1;
    this->nmodule = 0;
    this->modules.clear();

    this->branch_side.clear();
    this->branch_detid.clear();
    this->branch_module.clear();
    this->branch_disk.clear();
    this->branch_layer.clear();
    this->branch_ring.clear();

    this->branch_elink.clear();
    this->branch_lpgbt.clear();
    this->branch_dtc.clear();

    this->branch_row.clear();
    this->branch_column.clear();
    this->branch_adc.clear();

    this->branch_hlt_nclusters.clear();
    this->branch_bril_nclusters.clear();

    this->branch_bril_hit_stage.clear();

    this->branch_bril_cluster_nbits.clear();
    this->branch_bril_cluster_column.clear();
    this->branch_bril_cluster_row.clear();
    this->branch_bril_cluster_size_x.clear();
    this->branch_bril_cluster_size_y.clear();

    this->branch_hlt_cluster_nbits.clear();
    this->branch_hlt_cluster_column.clear();
    this->branch_hlt_cluster_row.clear();
    this->branch_hlt_cluster_size_x.clear();
    this->branch_hlt_cluster_size_y.clear();

    this->branch_col_not_adjacent.clear();
    this->branch_col_already_used.clear();
    this->branch_col_not_found.clear();
    this->branch_col_rhs_must_be_bigger.clear();
    this->branch_row_processor_buffer_not_empty.clear();
    this->branch_row_already_used.clear();
    this->branch_row_not_found.clear();
  }
  void attach_tree(TTree *tree)
  {
    this->tree = tree;
    tree->Branch("event", &this->event, "event/I");
    tree->Branch("nmodule", &this->nmodule, "nmodule/I");

    tree->Branch("detid", &this->branch_detid);
    tree->Branch("module", &this->branch_module);
    tree->Branch("disk", &this->branch_disk);
    tree->Branch("layer", &this->branch_layer);
    tree->Branch("ring", &this->branch_ring);

    tree->Branch("side", &this->branch_side);
    tree->Branch("elink", &this->branch_elink);
    tree->Branch("lpgbt", &this->branch_lpgbt);
    tree->Branch("dtc", &this->branch_dtc);

    tree->Branch("row", &this->branch_row);
    tree->Branch("column", &this->branch_column);
    tree->Branch("adc", &this->branch_adc);

    tree->Branch("hlt_nClusters", &this->branch_hlt_nclusters);

    tree->Branch("hlt_cluster_row", &this->branch_hlt_cluster_row);
    tree->Branch("hlt_cluster_column", &this->branch_hlt_cluster_column);
    tree->Branch("hlt_cluster_size_x", &this->branch_hlt_cluster_size_x);
    tree->Branch("hlt_cluster_size_y", &this->branch_hlt_cluster_size_y);
    tree->Branch("hlt_cluster_nBits", &this->branch_hlt_cluster_nbits);

    tree->Branch("bril_nClusters", &this->branch_bril_nclusters);
    tree->Branch("bril_hit_stage", &this->branch_bril_hit_stage);

    tree->Branch("bril_cluster_row", &this->branch_bril_cluster_row);
    tree->Branch("bril_cluster_column", &this->branch_bril_cluster_column);
    tree->Branch("bril_cluster_size_x", &this->branch_bril_cluster_size_x);
    tree->Branch("bril_cluster_size_y", &this->branch_bril_cluster_size_y);
    tree->Branch("bril_cluster_nBits", &this->branch_bril_cluster_nbits);

    tree->Branch("col_not_adjacent", &this->branch_col_not_adjacent);
    tree->Branch("col_already_used", &this->branch_col_already_used);
    tree->Branch("col_not_found", &this->branch_col_not_found);
    tree->Branch("col_rhs_must_be_bigger", &this->branch_col_rhs_must_be_bigger);
    tree->Branch("row_przocessor_buffer_not_empty", &this->branch_row_processor_buffer_not_empty);
    tree->Branch("row_already_used", &this->branch_row_already_used);
    tree->Branch("row_not_found", &this->branch_row_not_found);
  }
  void serialize()
  {
    this->nmodule = this->modules.size();
    for (auto module : this->modules)
    {
      this->branch_detid.push_back(module.detid);
      this->branch_module.push_back(module.module);
      this->branch_disk.push_back(module.disk);
      this->branch_layer.push_back(module.layer);
      this->branch_ring.push_back(module.ring);
      this->branch_side.push_back(module.side);

      this->branch_hlt_nclusters.push_back(module.hlt_nclusters);
      this->branch_bril_nclusters.push_back(module.bril_nclusters);

      this->branch_elink.push_back(module.elink);
      this->branch_lpgbt.push_back(module.lpgbt);
      this->branch_dtc.push_back(module.dtc);

      this->branch_row.push_back(module.row);
      this->branch_column.push_back(module.column);
      this->branch_adc.push_back(module.adc);

      this->branch_hlt_cluster_column.push_back(module.hlt_cluster_column);
      this->branch_hlt_cluster_row.push_back(module.hlt_cluster_row);
      this->branch_hlt_cluster_size_x.push_back(module.hlt_size_x);
      this->branch_hlt_cluster_size_y.push_back(module.hlt_size_y);
      this->branch_hlt_cluster_nbits.push_back(module.hlt_nbits);

      this->branch_bril_hit_stage.push_back(module.bril_hit_stage);
      this->branch_bril_cluster_row.push_back(module.bril_cluster_row);
      this->branch_bril_cluster_column.push_back(module.bril_cluster_column);
      this->branch_bril_cluster_size_x.push_back(module.bril_cluster_size_x);
      this->branch_bril_cluster_size_y.push_back(module.bril_cluster_size_y);
      this->branch_bril_cluster_nbits.push_back(module.bril_cluster_nbits);

      this->branch_col_not_adjacent.push_back(module.col_not_adjacent);
      this->branch_col_already_used.push_back(module.col_already_used);
      this->branch_col_not_found.push_back(module.col_not_found);
      this->branch_col_rhs_must_be_bigger.push_back(module.col_rhs_must_be_bigger);
      this->branch_row_processor_buffer_not_empty.push_back(module.row_processor_buffer_not_empty);
      this->branch_row_already_used.push_back(module.row_already_used);
      this->branch_row_not_found.push_back(module.row_not_found);
    }

    this->tree->Fill();
  }
};

class ITdigiExporter : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit ITdigiExporter(const edm::ParameterSet &);
  ~ITdigiExporter();

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event &, const edm::EventSetup &) override;
  virtual void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> m_tokenDigis;            // digi
  edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster>> m_tokenClusters; // clusters

  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> m_tokenGeometry; // geometry
  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> m_tokenTopo;            // geometry
  edm::ESGetToken<TrackerDetToDTCELinkCablingMap, TrackerDetToDTCELinkCablingMapRcd> m_tokenMap;

  // the pointers to geometry, topology and clusters
  // these are members so all functions can access them without passing as
  // argument
  const TrackerTopology *tTopo = NULL;
  const TrackerGeometry *tkGeom = NULL;
  const edm::DetSetVector<PixelDigi> *digis = NULL;
  const edmNew::DetSetVector<SiPixelCluster> *clusters = NULL;

  // File service to write ROOT output
  edm::Service<TFileService> m_fileservice;
  ITDigiEvent m_event;
  TTree *m_tree;
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

ITdigiExporter::ITdigiExporter(const edm::ParameterSet &iConfig) : m_tokenDigis(consumes<edm::DetSetVector<PixelDigi>>(iConfig.getParameter<edm::InputTag>("digis"))),
                                                                   m_tokenClusters(consumes<edmNew::DetSetVector<SiPixelCluster>>(iConfig.getParameter<edm::InputTag>("clusters"))),
                                                                   m_tokenGeometry(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>(edm::ESInputTag("", "idealForDigi"))),
                                                                   m_tokenTopo(esConsumes<TrackerTopology, TrackerTopologyRcd>()),
                                                                   m_tokenMap(esConsumes<TrackerDetToDTCELinkCablingMap, TrackerDetToDTCELinkCablingMapRcd>())
{
  // TTree version 0: Save everything into one TTree
  // TODO: Organize data from different rings, detector parts
  // either into same TTree or different TTrees
  usesResource("TFileService");
  m_tree = m_fileservice->make<TTree>("Digis", "Digis");
  this->m_event.attach_tree(this->m_tree);
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
void ITdigiExporter::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  // Empty event container
  this->m_event.clear();
  this->m_event.event = iEvent.id().event();

  // get the digis - COB 26.02.19
  edm::Handle<edm::DetSetVector<PixelDigi>> tdigis;
  iEvent.getByToken(m_tokenDigis, tdigis);

  edm::Handle<edmNew::DetSetVector<SiPixelCluster>> tcluster;
  iEvent.getByToken(m_tokenClusters, tcluster);

  // Get the geometry
  edm::ESHandle<TrackerGeometry> tgeomHandle = iSetup.getHandle(m_tokenGeometry);

  // Get the topology
  edm::ESHandle<TrackerTopology> tTopoHandle = iSetup.getHandle(m_tokenTopo);

  // Get the cabling map
  edm::ESHandle<TrackerDetToDTCELinkCablingMap> cablingMapHandle = iSetup.getHandle(m_tokenMap);

  TrackerDetToDTCELinkCablingMap const *cablingMap = cablingMapHandle.product();

  // get the pointers to geometry and topology
  tTopo = tTopoHandle.product();
  const TrackerGeometry *tkGeom = &(*tgeomHandle);
  tkGeom = tgeomHandle.product();
  digis = tdigis.product(); // pointer to digis - COB 26.02.19
  clusters = tcluster.product();

  // Loop over modules
  for (typename edm::DetSetVector<PixelDigi>::const_iterator DSV_it = digis->begin(); DSV_it != digis->end(); DSV_it++)
  {

    ITDigiModule imodule;

    // get the detid
    unsigned int rawid = DSV_it->detId();

    DetId detId(rawid);

    imodule.detid = rawid;

    // Determine whether it's barrel or endcap
    TrackerGeometry::ModuleType mType = tkGeom->getDetectorType(detId);

    if (!(mType == TrackerGeometry::ModuleType::Ph2PXF || mType == TrackerGeometry::ModuleType::Ph2PXF3D) && detId.subdetId() == PixelSubdetector::PixelEndcap)
      continue;

    // obtaining location of module
    if (tTopo->pxfDisk(detId) <= 8)
      continue;

    imodule.disk = tTopo->pxfDisk(detId);
    imodule.ring = tTopo->pxfBlade(detId);
    imodule.module = tTopo->pxfModule(detId);
    imodule.side = tTopo->pxfSide(detId);

    imodule.layer = 1;

    for (auto [min_id, max_id] : backside_panel_ids)
      imodule.layer = (detId > min_id && detId < max_id) ? 2 : imodule.layer;

    // Find the geometry of the module associated to the digi
    const GeomDetUnit *geomDetUnit(tkGeom->idToDetUnit(detId));

    if (!geomDetUnit)
      continue;

    // Loop over the digis in each module
    auto elink_ids = cablingMap->detIdToDTCELinkId(detId);

    imodule.elink = (*elink_ids.first).second.elink_id();
    imodule.lpgbt = (*elink_ids.first).second.gbtlink_id();
    imodule.dtc = (*elink_ids.first).second.dtc_id();

    imodule.hlt_nclusters = 0;

    std::map<vertex_t, int> pixel_mapping;

    size_t pixel_index = 0;

    // if (!DSV_it->empty())
    // {
    //   for (edm::DetSet<PixelDigi>::const_iterator digit = DSV_it->begin(); digit != DSV_it->end(); digit++)
    //   {
    //     int32_t r = digit->row(), c = digit->column();

    //     if (r < row_cutoff / 2 && r >= (row_max - row_cutoff / 2))
    //       continue;
    //     if (c < column_cutoff / 2 && c >= (column_max - column_cutoff / 2))
    //       continue;

    //     r -= row_cutoff / 2;
    //     c -= column_cutoff / 2;

    //     if (r < 1 || r >= (actual_row_max - 1))
    //       continue;
    //     if (c < 1 || c >= (actual_column_max - 1))
    //       continue;
    //     if (r >= (second_row_chip_start - 1) && r < (second_row_chip_start + 1))
    //       continue;
    //     if (c >= (second_col_chip_start - 1) && c < (second_col_chip_start + 1))
    //       continue;

    //     imodule.row.push_back(r);
    //     imodule.column.push_back(c);

    //     pixel_mapping[std::make_pair(c, r)] = i++;

    //     imodule.adc.push_back(digit->adc());
    //   }
    // }

    if (!clusters->empty())
    {
      auto check_id = [&rawid](const auto &item) -> bool
      {
        return item.detId() == rawid;
      };

      typename edmNew::DetSetVector<SiPixelCluster>::const_iterator cluster_it = std::find_if(clusters->begin(), clusters->end(), check_id);

      if (cluster_it == clusters->end())
        continue;

      imodule.hlt_nclusters = cluster_it->size();

      for (edmNew::DetSet<SiPixelCluster>::const_iterator cluit = cluster_it->begin(); cluit != cluster_it->end(); cluit++)
      {
        std::vector<unsigned short> pixels_x, pixels_y;

        const size_t n_pixels = cluit->size();

        
        for (size_t i = 0; i < n_pixels; i++)
        {
          uint16_t c = cluit->pixel(i).x, r = cluit->pixel(i).y;

          pixels_x.push_back(c);
          pixels_y.push_back(r);

          imodule.row.push_back(c);
          imodule.column.push_back(r);
          pixel_mapping[std::make_pair((int)c, (int)r)] = pixel_index++;

        }

        for (const auto & adc_value : cluit->pixelADC())
        {
          imodule.adc.push_back(adc_value);
        }

        int16_t min_x = pixels_x.size() > 0 ? (int16_t)*std::min_element(pixels_x.begin(), pixels_x.end()) : 0;
        int16_t min_y = pixels_y.size() > 0 ? (int16_t)*std::min_element(pixels_y.begin(), pixels_y.end()) : 0;

        auto &r = min_x, &c = min_y;

        // if (r < row_cutoff / 2 && r >= (row_max - row_cutoff / 2))
        //   continue;
        // if (c < column_cutoff / 2 && c >= (column_max - column_cutoff / 2))
        //   continue;

        // r -= row_cutoff / 2;
        // c -= column_cutoff / 2;

        // if (r < 1 || r >= (actual_row_max - 1))
        //   continue;
        // if (c < 1 || c >= (actual_column_max - 1))
        //   continue;
        // if (r >= (second_row_chip_start - 1) && r < (second_row_chip_start + 1))
        //   continue;
        // if (c >= (second_col_chip_start - 1) && c < (second_col_chip_start + 1))
        //   continue;

        // if (c < 0)
        //   imodule.hlt_cluster_column.push_back(0);
        // else
        imodule.hlt_cluster_column.push_back((uint16_t)c);

        // if (r < 0)
        //   imodule.hlt_cluster_row.push_back(0);
        // else
        imodule.hlt_cluster_row.push_back((uint16_t)r);

        imodule.hlt_size_x.push_back(cluit->sizeY());
        imodule.hlt_size_y.push_back(cluit->sizeX());
        imodule.hlt_nbits.push_back(cluit->size());
      }
    }
    else
    {
#ifdef DEBUG
      std::cout << "no clusters!" << std::endl;
#endif
      imodule.hlt_nclusters = 0;
      imodule.hlt_cluster_column.clear();
      imodule.hlt_cluster_row.clear();
      imodule.hlt_nbits.clear();
      imodule.hlt_size_x.clear();
      imodule.hlt_size_y.clear();
    }
#ifdef DEBUG
    std::cout << "retrieved HLT clusteres" << std::endl;
#endif

#ifdef DEBUG
    std::cout << "retrieved pixel data" << std::endl;
#endif

    std::map<vertex_t, QuarterCore> qcores_map;

    for (const auto &pix : pixel_mapping)
    {
      const auto &hit = pix.first;

      int32_t col_in_qcore = hit.first % SIZE_QCORE_HORIZONTAL;
      int32_t row_in_qcore = hit.second % SIZE_QCORE_VERTICAL;

      int32_t qcol = (hit.first - col_in_qcore) / SIZE_QCORE_HORIZONTAL;
      int32_t qrow = (hit.second - row_in_qcore) / SIZE_QCORE_VERTICAL;

      vertex_t qpos = std::make_pair(qcol, qrow);

      if (qcores_map.find(qpos) == qcores_map.end())
        qcores_map[qpos] = QuarterCore(qcol, qrow);

      qcores_map[qpos].add_hit(col_in_qcore, row_in_qcore);
    }

#ifdef DEBUG
    std::cout << "determined used qcores" << std::endl;
#endif

    std::vector<QuarterCore> qcores;

    for (const auto &item : qcores_map)
    {
      qcores.push_back(item.second);
    }

    auto qcores_metadata = determine_metadeta(qcores);

#ifdef DEBUG
    std::cout << "metadata is determined" << std::endl;
#endif

    std::vector<ProcessedQuarterCoreModel> qcore_processed;

    for (auto qcore : qcores_metadata)
      qcore_processed.push_back(ProcessedQuarterCoreModel(qcore));

#ifdef DEBUG
    std::cout << "Qcores are processed" << std::endl;
#endif

    DistributorModel distributor(qcore_processed);

    auto [single_clusters, row_buffer] = distributor.run();

#ifdef DEBUG
    std::cout << "Distributor is done" << std::endl;
#endif

    RowMergerModel row_merger(row_buffer);

    auto [row_clusters, col_buffer] = row_merger.run();

#ifdef DEBUG
    std::cout << "row merger is done" << std::endl;
#endif

    ColMergerModel col_merger(col_buffer);

    auto col_clusters = col_merger.run();

#ifdef DEBUG
    std::cout << "col merger is done" << std::endl;
#endif

    imodule.bril_hit_stage = std::vector<unsigned char>(imodule.adc.size(), 0);

#ifdef DEBUG
    std::cout << "bril hits array is allocated" << std::endl;
#endif

    auto fill_data = [&imodule, &pixel_mapping](const std::vector<ClusterModel> &data, char stage) -> void
    {
      for (auto cluster : data)
      {
        imodule.bril_cluster_row.push_back(cluster.row);
        imodule.bril_cluster_column.push_back(cluster.col);

        std::vector<unsigned short> cluster_x, cluster_y;

        for (const auto &hit : cluster.hit_map)
        {
          const auto [x, y] = hit;

          vertex_t global_coordinate = std::make_pair(x + SIZE_QCORE_HORIZONTAL * cluster.col, y + SIZE_QCORE_VERTICAL * cluster.row);

          if (pixel_mapping.find(global_coordinate) != pixel_mapping.end())
            imodule.bril_hit_stage[pixel_mapping[global_coordinate]] = stage;
        }

        imodule.bril_cluster_size_x.push_back(*std::max_element(cluster_x.begin(), cluster_x.end()) - *std::min_element(cluster_x.begin(), cluster_x.end()));
        imodule.bril_cluster_size_y.push_back(*std::max_element(cluster_y.begin(), cluster_y.end()) - *std::min_element(cluster_y.begin(), cluster_y.end()));

        imodule.bril_cluster_nbits.push_back(cluster.hit_map.size());
      }
    };

    fill_data(single_clusters, 1);
#ifdef DEBUG
    std::cout << "single_clusters is added to data" << std::endl;
#endif
    fill_data(row_clusters, 2);
#ifdef DEBUG
    std::cout << "row_clusters is added to data" << std::endl;
#endif
    fill_data(col_clusters, 3);
#ifdef DEBUG
    std::cout << "col_clusters is added to data" << std::endl;
#endif
    imodule.bril_nclusters = single_clusters.size() + row_clusters.size() + col_clusters.size();

    imodule.col_not_adjacent = col_merger.get_not_adjacent();
    imodule.col_already_used = col_merger.get_already_used();
    imodule.col_not_found = col_merger.get_not_found();
    imodule.col_rhs_must_be_bigger = col_merger.get_rhs_must_be_bigger();
    imodule.row_processor_buffer_not_empty = row_merger.get_processor_buffer_not_empty();
    imodule.row_already_used = row_merger.get_already_used();
    imodule.row_not_found = row_merger.get_not_found();
#ifdef DEBUG
    std::cout << "module cluster metadata is added to data" << std::endl;
#endif
    this->m_event.modules.push_back(imodule);
  }

  this->m_event.serialize();
}

// ------------ method called once each job just before starting event loop
// ------------
void ITdigiExporter::beginJob() {}

// ------------ method called once each job just after ending the event loop
// ------------
void ITdigiExporter::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the
// module  ------------
void ITdigiExporter::fillDescriptions(
    edm::ConfigurationDescriptions &descriptions)
{
  // The following says we do not know what parameters are allowed so do no
  // validation
  //  Please change this to state exactly what you do use, even if it is no
  //  parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  // Specify that only 'tracks' is allowed
  // To use, remove the default given above and uncomment below
  // ParameterSetDescription desc;
  // desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  // descriptions.addDefault(desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(ITdigiExporter);
