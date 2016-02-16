#ifndef LARLITE_SHOWERANA_CXX
#define LARLITE_SHOWERANA_CXX

#include "ShowerAna.h"

namespace larlite {
  
  // Called at beginning of event loop
  bool ShowerAna::initialize() {

    // Not sure if events work in same way for reco2D
    _evtN = 0;

    // Initialise TTree with all relevant branches
    _t_ch = new TTree("ch_tree","");
    _t_ch->Branch("evt",&_evtN,"evt/i");
    _t_ch->Branch("hitNo",&_hitNo,"hitNo/i");
    _t_ch->Branch("hitNoU",&_hitNoU,"hitNoU/i");
    _t_ch->Branch("hitNoV",&_hitNoV,"hitNoV/i");
    _t_ch->Branch("hitNoY",&_hitNoY,"hitNoY/i");
    _t_ch->Branch("TDCstd",&_TDCstd,"TDCstd/F");
    _t_ch->Branch("TDCstdU",&_TDCstdU,"TDCstdU/F");
    _t_ch->Branch("TDCstdV",&_TDCstdV,"TDCstdV/F");
    _t_ch->Branch("TDCstdY",&_TDCstdY,"TDCstdY/F");
    _t_ch->Branch("TDCiqr",&_TDCiqr,"TDCiqr/F");
    _t_ch->Branch("TDCiqrU",&_TDCiqrU,"TDCiqrU/F");
    _t_ch->Branch("TDCiqrV",&_TDCiqrV,"TDCiqrV/F");
    _t_ch->Branch("TDCiqrY",&_TDCiqrY,"TDCiqrY/F");
    _t_ch->Branch("ADCamp",&_ADCamp,"ADCamp/F");
    _t_ch->Branch("ADCampU",&_ADCampU,"ADCampU/F");
    _t_ch->Branch("ADCampV",&_ADCampV,"ADCampV/F");
    _t_ch->Branch("ADCampY",&_ADCampY,"ADCampY/F");
    _t_ch->Branch("WFint",&_WFint,"WFint/F");
    _t_ch->Branch("WFintU",&_WFintU,"WFintU/F");
    _t_ch->Branch("WFintV",&_WFintV,"WFintV/F");
    _t_ch->Branch("WFintY",&_WFintY,"WFintY/F");
    _t_ch->Branch("showerNo",&_showerNo,"showerNo/i");
    _t_ch->Branch("energy",&_Energy,"Energy/F");
    _t_ch->SetDirectory(0);

    return true;
  }
  
  // Called for every event
  bool ShowerAna::analyze(storage_manager* storage) {
   
    bool showerFlag = 1;

    // Get MCShower information
    auto ev_mcs = storage->get_data<event_mcshower>("mcreco");
    // Check if it exists / are there any showers
    if ( (!ev_mcs) || (!ev_mcs->size())){
    	print (msg::kERROR,__FUNCTION__,"MCShower data product not found!");
        // Probably no point doing any analysis if no showers are present
        // might want a control
        //return false;
        // If you still want to do stuff remove line above
    	showerFlag = 0;
    }

    // Re set shower count
    _showerNo = 0;

    if (showerFlag){
      // Loop over number of showers? not sure this is how MCShower works
      for (size_t j=0; j < ev_mcs->size(); ++j){
        auto const& shower = (*ev_mcs).at(j);
        std::cout << j << std::endl;
        // Get energy, not sure if this is right energy
        _Energy = shower.Start().E();
        // Count number of showers per event
        _showerNo += 1;
      }
    }


    // Use storage to get larlite::hit object (gaushit, cccluster, pandoraCosmicKHitRemoval)
    auto hitdata = storage->get_data<event_hit>("gaushit");
    // Display error if hit data not present
    if ( (!hitdata) || (!hitdata->size())){
    	print (msg::kERROR,__FUNCTION__,"Hit data product not found!");
    	return false;
    }

    // Initialise hit counters for event
    _isHit = 0;
    uHit = 0;
    vHit = 0;
    yHit = 0;

    // Initialize integration counter for event
    intADC = 0;
    UintADC = 0;
    VintADC = 0;
    YintADC = 0;

    // Clear TDC vectors
    TDCvec.clear();
    UTDCvec.clear();
    VTDCvec.clear();
    YTDCvec.clear();

    // Clear ADC vectors
    ADCvec.clear();
    UADCvec.clear();
    VADCvec.clear();
    YADCvec.clear();

    // Get channel ID
    for (size_t i=0; i < hitdata->size(); ++i){
      auto const& hit = (*hitdata).at(i);
      int chnum = hit.Channel();
      _isHit += 1;
      if(chnum<2400) uHit += 1;
      if(chnum>=2400&&chnum<4800) vHit += 1;
      if(chnum>=4800&&chnum<8256) yHit += 1;
    }
    // Count up total number of hits per event
    // Add total number of hits per event to ttree
    _hitNo = _isHit;
    _hitNoU = uHit;
    _hitNoV = vHit;
    _hitNoY = yHit;

    _t_ch->Fill();

    _evtN += 1;
  
    return true;
  }

  // Called at end of event loop
  bool ShowerAna::finalize() {

    // Save histograms to .root file
   if(_fout){
     _fout->cd();
     std::cout << "writing ch tree" << std::endl;
     _t_ch->Write();
   }
  
    return true;
  }

}
#endif
