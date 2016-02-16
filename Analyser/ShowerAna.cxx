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
    _t_ch->SetDirectory(0);

    return true;
  }
  
  // Called for every event
  bool ShowerAna::analyze(storage_manager* storage) {
  
    // Use storage to get larlite::hit object (one per event?)
    auto hitdata = storage->get_data<event_hit>("gaushit");
    // Display error if hit data not present
    if ( (!hitdata) || (!hitdata->size())){
    	print (msg::kERROR,__FUNCTION__,"Hit data product not found!");
    	return false;
    }
    _isHit = 0;
    // Get channel ID
    for (size_t i=0; i < hitdata->size(); ++i){
      auto const& hit = (*hitdata).at(i);
      int chnum = hit.Channel();
      std::cout<<chnum<<std::endl;
      _isHit += 1;
    }
    // Count up total number of hits per event
    // Add total number of hits per event to ttree
    _hitNo = _isHit;
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
