#ifndef LARLITE_ROIFINDER_CXX
#define LARLITE_ROIFINDER_CXX

#include "ROIfinder.h"

namespace larlite {

  bool ROIfinder::initialize() {

    // Not sure if events work in same way for reco2D
    _evtN = 0;
    // Initialise TTree with all relevant branches
    _t_ch = new TTree("ch_tree","");
    _t_ch->Branch("ADCvec",&ADCvec);
    _t_ch->Branch("TDCvec",&TDCvec);
    _t_ch->Branch("ROIvec",&ROIvec);
    _t_ch->SetDirectory(0);

    return true;
  }
  
  bool ROIfinder::analyze(storage_manager* storage) {
  
    // Get data
    auto hitdata = storage->get_data<event_hit>("gaushit");
    if ( (!hitdata) || (!hitdata->size())){
    	print (msg::kERROR,__FUNCTION__,"Hit data product not found!");
    	return false;
    }

    TDCvec.clear();
    ADCvec.clear();
    ROIvec.clear();
    // Loop over once
    for (size_t i=0; i < hitdata->size(); ++i){
      auto const& hit = (*hitdata).at(i);
      float amp = hit.PeakAmplitude();
      float tdc = hit.PeakTime();
      // Get time of all hits > 15 ADCs
      if(amp>15){ROIvec.push_back(tdc);}
    }
    // Group all hits within 300 ticks of each other
    std::sort(ROIvec.begin(),ROIvec.end());
    std::vector<float> start, end;
    bool is_start = false;
    start.push_back(ROIvec[0]);
    for (size_t j=0; j < ROIvec.size(); ++j){
      if (j<ROIvec.size()-2){
      if (ROIvec[j+1]-ROIvec[j]>300&&ROIvec[j+2]-ROIvec[j+1]>50){
        if(is_start){
          start.push_back(ROIvec[j]);
          is_start = false;
        }
        else {
          end.push_back(ROIvec[j]);
          is_start = true;
        }
      }
      }
      else end.push_back(ROIvec[j]);
    }
    // If <= 3 hits within 300 ticks then ignore and call last point start/end of ROI
    // add start & end of ROIs to a vector
    
    ROIvec.clear();

    // Loop over data as many times as ROIs
    for (size_t k=0; k < start.size(); ++k){
      for (size_t l=0; l < hitdata->size(); ++l){
        auto const& hit = (*hitdata).at(l);
        float amp = hit.PeakAmplitude();
        float tdc = hit.PeakTime();
        if (tdc>start[k]&&tdc<end[k]){
          ADCvec.push_back(amp);
          TDCvec.push_back(tdc);
          ROIvec.push_back(k);
        }
      }
    }
    // Calculate usual variables and push on to vector
    // Vectors will be as long as the number of ROIs
    // Fill the ttree with the vectors
    // Probably won't use in this way, do filtering here and just output decision
    // For each ROI see if variable(s) match shower or non shower 
    // Need to work out how to do this ^ (decision trees?)
    //
    _t_ch->Fill();

    _evtN += 1;
  
    return true;
  }

  bool ROIfinder::finalize() {

   if(_fout){
     _fout->cd();
     std::cout << "writing ch tree" << std::endl;
     _t_ch->Write();
   }
  
    return true;
  }

}
#endif
