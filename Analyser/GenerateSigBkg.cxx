#ifndef LARLITE_GENERATESIGBKG_CXX
#define LARLITE_GENERATESIGBKG_CXX

#include "GenerateSigBkg.h"

// Function to calculate standard deviation from vector of TDC times
float gTDCstd(std::vector<float> const & TDCvec) {

    float stdTDC(0);
    float meanTDC(0);

    for(size_t f=0; f<TDCvec.size(); ++f) {
        meanTDC += (TDCvec[f]);
    }
    meanTDC /= ((float)TDCvec.size());

    for(size_t k=0; k<TDCvec.size(); ++k)
      stdTDC += (TDCvec[k]-meanTDC)*(TDCvec[k]-meanTDC);
    stdTDC = sqrt( stdTDC / ((float)TDCvec.size()));
    return stdTDC;

}

// Function to calculate interquartile range from vector of TDC times
float gTDCiqr(std::vector<float> TDCvec, int hitNo) {

    std::sort(TDCvec.begin(),TDCvec.end());
    int quart = hitNo/4;
    float iqrTDC = (TDCvec[3*quart] - TDCvec[quart]);
    return iqrTDC;

}

float gMaxVal(std::vector<float> ADCvec){

    // Sort in reverse order
    std::sort(ADCvec.rbegin(),ADCvec.rend());
    float avmax = 0;
    if(ADCvec.size()>20){
      for(int i = 0; i < 20; ++i){
        avmax += ADCvec[i];
      }
      return avmax/20;
    } else return avmax;

}

namespace larlite {

  bool GenerateSigBkg::initialize() {

    // Initialise TTree for Y variables and 960x346 boxes
    sig_tree = new TTree("sig_tree","");
    sig_tree->Branch("evt",&_evtN,"evt/i");
    sig_tree->Branch("hitNo",&_hitNoY,"hitNo/i");
    sig_tree->Branch("TDCstd",&_TDCstdY,"TDCstd/F");
    sig_tree->Branch("TDCiqr",&_TDCiqrY,"TDCiqr/F");
    sig_tree->Branch("ADCamp",&_ADCampY,"ADCamp/F");
    sig_tree->Branch("WFint",&_WFintY,"WFint/F");
    sig_tree->Branch("Meanamp",&_MeanampY,"Meanamp/F");
    sig_tree->Branch("Meanint",&_MeanintY,"Meanint/F");
    sig_tree->Branch("LowDen",&_LowDenY,"LowDen/i");
    sig_tree->Branch("HiDen",&_HiDenY,"HiDen/i");
    sig_tree->Branch("MeanRMS",&_MeanRMSY,"MeanRMS/F");
    sig_tree->Branch("MeanMult",&_MeanMultY,"MeanMult/F");
    sig_tree->Branch("Wirestd",&_WirestdY,"Wirestd/F");
    sig_tree->Branch("Wireiqr",&_WireiqrY,"Wireiqr/F");
    sig_tree->Branch("NumMult",&_NumMultY,"NumMult/i");
    sig_tree->SetDirectory(0);

    // Initialise background TTree for Y variables and 960x346 boxes
    bkg_tree = sig_tree->CloneTree();
    bkg_tree->SetDirectory(0);
    bkg_tree->SetName("bkg_tree");

    return true;
  }
  
  bool GenerateSigBkg::analyze(storage_manager* storage) {
  
    bool showerFlag = 1;

    // Get MCShower information
    auto ev_mcs = storage->get_data<event_mcshower>("mcreco");
    // Check if it exists / are there any showers
    if ( (!ev_mcs) || (!ev_mcs->size())){
    	print (msg::kERROR,__FUNCTION__,"MCShower data product not found!");
        // If you still want to do stuff remove line above
    	showerFlag = 0;
    }

    // Re set shower count
    _showerNo = 0;
    _detShowerNo = 0;

    if (showerFlag){
      // Loop over number of showers? not sure this is how MCShower works
      for (size_t j=0; j < ev_mcs->size(); ++j){
        auto const& shower = (*ev_mcs).at(j);
        // Get energy, not sure if this is right energy
        float Energy = shower.Start().E();
        float startX = shower.Start().X()/100;
        float startY = shower.Start().Y()/100;
        float startZ = shower.Start().Z()/100;
        // Count number of showers per event
        if((startX>0&&startX<2.5604)&&(startY>-1.2&&startY<1.2)&&(startZ>0&&startZ<10.5)&&Energy>100){
          _detShowerNo += 1;
        }
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
    _hitNoY = 0;
    // Initialize integration counters for event
    _WFintY = 0;
    // Clear TDC vectors
    YTDCvec.clear();
    // Clear ADC vectors
    YADCvec.clear();
    // Clear Channel number vectors
    YChvec.clear();
    // Initialize TDC standard deviation counters
    _TDCstdY = 0;
    // Initialize TDC interquartile range counters
    _TDCiqrY = 0;
    // Initialize max (or total?) ADC amplitude counters
    _ADCampY = 0;
    // Initialize mean adc amplitude counters
    _MeanampY = 0;
    // Initailize mean integrated waveform counters
    _MeanintY = 0;
    // Initialize hit density between 5-35 ADCs
    _LowDenY = 0;
    // Initialize hit density > 35
    _HiDenY = 0;
    // Initialize mean RMS counters
    _MeanRMSY = 0;
    // Initialize mean multiplicity counters
    _MeanMultY = 0;
    // Initialize Wire standard deviation counter
    _WirestdY = 0;
    // Initialize Wire interquartile range counter
    _WireiqrY = 0;
    _NumMultY = 0;

    // Loop over hits in each event
    for (size_t i=0; i < hitdata->size(); ++i){
      auto const& hit = (*hitdata).at(i);
      // Get the channel number
      int chnum = hit.Channel();
      float amp = hit.PeakAmplitude();
      float wint = hit.Integral();
      float tdc = hit.PeakTime();
      float rms = hit.RMS();
      float mult = hit.Multiplicity();
      if(amp>10&&tdc>3600&&tdc<9600){
      // Do same for Y plane CHECK THIS
      if(chnum>=4800&&chnum<8256){
        _hitNoY += 1;
        _MeanampY += amp;
        _WFintY += wint;
        _MeanRMSY += rms;
        _MeanMultY += mult;
        YTDCvec.push_back(tdc);
        YADCvec.push_back(amp);
        YChvec.push_back(chnum);
        if(amp<=26)_LowDenY += 1;
        if(amp>26)_HiDenY += 1;
      }
      }
    }

    _ADCampY = gMaxVal(YADCvec);
    _MeanampY = _MeanampY/_hitNoY;
    _MeanintY = _WFintY/_hitNoY;
    _MeanRMSY = _MeanRMSY/_hitNoY;
    _MeanMultY = _MeanMultY/_hitNoY;
    // Y plane
    if(_hitNoY!=0){
      _TDCstdY = gTDCstd(YTDCvec);
      _TDCiqrY = gTDCiqr(YTDCvec,_hitNoY);
      _WirestdY = gTDCstd(YChvec);
      _WireiqrY = gTDCiqr(YChvec,_hitNoY);
    }

    int current_num_hits, next_num_hits;

    for(int k(4800); k<8256; k++){
      current_num_hits = 0;
      next_num_hits = 0;
      for(int l(0); l<YChvec.size(); l++){
        if(YChvec[l]==k) current_num_hits++;
        if(YChvec[l]==k+1) next_num_hits++;
      }
      if(std::abs(current_num_hits-next_num_hits)>3) _NumMultY++;
    }

    // Fill TTree
    if(_hitNoY==0||_ADCampY>2000) return false; 

    if(_detShowerNo==0) bkg_tree->Fill();
    else sig_tree->Fill();

    _evtN += 1;
  
    return true;
  }

  bool GenerateSigBkg::finalize() {

   // Save TTrees to .root file
   if(_fout){
     _fout->cd();
     std::cout << "writing ch tree" << std::endl;
     sig_tree->Write();
     bkg_tree->Write();
   }
  
    return true;
  }

}
#endif
