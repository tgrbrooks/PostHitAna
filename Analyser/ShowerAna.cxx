#ifndef LARLITE_SHOWERANA_CXX
#define LARLITE_SHOWERANA_CXX

#include "ShowerAna.h"

// Function to calculate standard deviation from vector of TDC times
float TDCstd(std::vector<float> const & TDCvec) {

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
float TDCiqr(std::vector<float> TDCvec, int hitNo) {

    std::sort(TDCvec.begin(),TDCvec.end());
    int quart = hitNo/4;
    float iqrTDC = (TDCvec[3*quart] - TDCvec[quart]);
    return iqrTDC;

}

float MaxVal(std::vector<float> ADCvec){

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
    _t_ch->Branch("Evec",&Evec);
    _t_ch->Branch("ShowerStartEnd",&ShowerStartEnd);

    _t_ch->Branch("Meanamp",&_Meanamp,"Meanamp/F");
    _t_ch->Branch("MeanampU",&_MeanampU,"MeanampU/F");
    _t_ch->Branch("MeanampV",&_MeanampV,"MeanampV/F");
    _t_ch->Branch("MeanampY",&_MeanampY,"MeanampY/F");

    _t_ch->Branch("Meanint",&_Meanint,"Meanint/F");
    _t_ch->Branch("MeanintU",&_MeanintU,"MeanintU/F");
    _t_ch->Branch("MeanintV",&_MeanintV,"MeanintV/F");
    _t_ch->Branch("MeanintY",&_MeanintY,"MeanintY/F");

    _t_ch->Branch("LowDen",&_LowDen,"LowDen/F");
    _t_ch->Branch("LowDenU",&_LowDenU,"LowDenU/F");
    _t_ch->Branch("LowDenV",&_LowDenV,"LowDenV/F");
    _t_ch->Branch("LowDenY",&_LowDenY,"LowDenY/F");

    _t_ch->Branch("HiDen",&_HiDen,"HiDen/F");
    _t_ch->Branch("HiDenU",&_HiDenU,"HiDenU/F");
    _t_ch->Branch("HiDenV",&_HiDenV,"HiDenV/F");
    _t_ch->Branch("HiDenY",&_HiDenY,"HiDenY/F");

    _t_ch->Branch("MeanRMS",&_MeanRMS,"MeanRMS/F");
    _t_ch->Branch("MeanRMSU",&_MeanRMSU,"MeanRMSU/F");
    _t_ch->Branch("MeanRMSV",&_MeanRMSV,"MeanRMSV/F");
    _t_ch->Branch("MeanRMSY",&_MeanRMSY,"MeanRMSY/F");

    _t_ch->Branch("MeanMult",&_MeanMult,"MeanMult/F");
    _t_ch->Branch("MeanMultU",&_MeanMultU,"MeanMultU/F");
    _t_ch->Branch("MeanMultV",&_MeanMultV,"MeanMultV/F");
    _t_ch->Branch("MeanMultY",&_MeanMultY,"MeanMultY/F");

    _t_ch->Branch("Wirestd",&_Wirestd,"Wirestd/F");
    _t_ch->Branch("WirestdU",&_WirestdU,"WirestdU/F");
    _t_ch->Branch("WirestdV",&_WirestdV,"WirestdV/F");
    _t_ch->Branch("WirestdY",&_WirestdY,"WirestdY/F");

    _t_ch->Branch("Wireiqr",&_Wireiqr,"Wireiqr/F");
    _t_ch->Branch("WireiqrU",&_WireiqrU,"WireiqrU/F");
    _t_ch->Branch("WireiqrV",&_WireiqrV,"WireiqrV/F");
    _t_ch->Branch("WireiqrY",&_WireiqrY,"WireiqrY/F");
   
    _t_ch->SetDirectory(0);

    return true;
  }
  
  // Called for every event
  bool ShowerAna::analyze(storage_manager* storage) {
   
    //float DriftVel = 1600;

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

    ShowerStartEnd.clear();
    Evec.clear();

    if (showerFlag){
      // Loop over number of showers? not sure this is how MCShower works
      for (size_t j=0; j < ev_mcs->size(); ++j){
        auto const& shower = (*ev_mcs).at(j);
        std::cout << j << std::endl;
        // Get energy, not sure if this is right energy
        Evec.push_back(shower.Start().E());
        float startX = shower.Start().X()/100;
        float endX = shower.End().X()/100;
        if((startX>0&&startX<2.5604)&&(endX>0&&endX<2.5604)){
          ShowerStartEnd.push_back(std::make_pair(startX, endX));
        }
        else ShowerStartEnd.push_back(std::make_pair(0, 0));
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
    _hitNo = 0;_hitNoU = 0;_hitNoV = 0;_hitNoY = 0;
    // Initialize integration counters for event
    _WFint = 0;_WFintU = 0;_WFintV = 0;_WFintY = 0;
    // Clear TDC vectors
    TDCvec.clear();UTDCvec.clear();VTDCvec.clear();YTDCvec.clear();
    // Clear ADC vectors
    ADCvec.clear();UADCvec.clear();VADCvec.clear();YADCvec.clear();
    // Clear Channel number vectors
    Chvec.clear();UChvec.clear();VChvec.clear();YChvec.clear();
    // Initialize TDC standard deviation counters
    _TDCstd = 0;_TDCstdU = 0;_TDCstdV = 0;_TDCstdY = 0;
    // Initialize TDC interquartile range counters
    _TDCiqr = 0;_TDCiqrU = 0;_TDCiqrV = 0;_TDCiqrY = 0;
    // Initialize max (or total?) ADC amplitude counters
    _ADCamp = 0;_ADCampU = 0;_ADCampV = 0;_ADCampY = 0;
    // Initialize mean adc amplitude counters
    _Meanamp = 0;_MeanampU = 0;_MeanampV = 0;_MeanampY = 0;
    // Initailize mean integrated waveform counters
    _Meanint = 0;_MeanintU = 0;_MeanintV = 0;_MeanintY = 0;
    // Initialize hit density between 5-35 ADCs
    _LowDen = 0;_LowDenU = 0;_LowDenV = 0;_LowDenY = 0;
    // Initialize hit density > 35
    _HiDen = 0;_HiDenU = 0;_HiDenV = 0;_HiDenY = 0;
    // Initialize mean RMS counters
    _MeanRMS = 0;_MeanRMSU = 0;_MeanRMSV = 0;_MeanRMSY = 0;
    // Initialize mean multiplicity counters
    _MeanMult = 0;_MeanMultU = 0;_MeanMultV = 0;_MeanMultY = 0;
    // Initialize Wire standard deviation counter
    _Wirestd = 0;_WirestdU = 0;_WirestdV = 0;_WirestdY = 0;
    // Initialize Wire interquartile range counter
    _Wireiqr = 0;_WireiqrU = 0;_WireiqrV = 0;_WireiqrY = 0;

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

      // Count number of hits
      _hitNo += 1;
      // Add up the hit amplitudes
      _Meanamp += amp;
      // Add up the integrated waveforms
      _WFint += wint;
      // Add up the peak RMS'
      _MeanRMS += rms;
      // Add up the multiplicities
      _MeanMult += mult;
      if(amp>5 && amp<=35) _LowDen += 1;
      if(amp>35) _HiDen += 1;
      // Record times, amplitudes and channels of each hit
      TDCvec.push_back(tdc);
      ADCvec.push_back(amp);
      Chvec.push_back(chnum);

      // Do same for U plane CHECK THIS
      if(chnum<2400){ 
        _hitNoU += 1;
        _MeanampU += amp;
        _WFintU += wint;
        _MeanRMSU += rms;
        _MeanMultU += mult;
        UTDCvec.push_back(tdc);
        UADCvec.push_back(amp);
        UChvec.push_back(chnum);
        if(amp>5 && amp<=35)_LowDenU += 1;
        if(amp>35)_HiDenU += 1;
      }

      // Do same for V plane CHECK THIS
      if(chnum>=2400&&chnum<4800){
        _hitNoV += 1;
        _MeanampV += amp;
        _WFintV += wint;
        _MeanRMSV += rms;
        _MeanMultV += mult;
        VTDCvec.push_back(tdc);
        VADCvec.push_back(amp);
        VChvec.push_back(chnum);
        if(amp>5 && amp<=35)_LowDenV += 1;
        if(amp>35)_HiDenV += 1;
      }

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
        if(amp>5 && amp<=35)_LowDenY += 1;
        if(amp>35)_HiDenY += 1;
      }

    }

    _LowDen = _LowDen/_hitNo;
    _LowDenU = _LowDenU/_hitNoU;
    _LowDenV = _LowDenV/_hitNoV;
    _LowDenY = _LowDenY/_hitNoY;

    _HiDen = _HiDen/_hitNo;
    _HiDenU = _HiDenU/_hitNoU;
    _HiDenV = _HiDenV/_hitNoV;
    _HiDenY = _HiDenY/_hitNoY;

    _ADCamp = MaxVal(ADCvec);
    _ADCampU = MaxVal(UADCvec);
    _ADCampV = MaxVal(VADCvec);
    _ADCampY = MaxVal(YADCvec);

    // Average the ADC amplitudes
    _Meanamp = _Meanamp/_hitNo;
    _MeanampU = _MeanampU/_hitNoU;
    _MeanampV = _MeanampV/_hitNoV;
    _MeanampY = _MeanampY/_hitNoY;

    // Average the integrated waveforms
    _Meanint = _WFint/_hitNo;
    _MeanintU = _WFintU/_hitNoU;
    _MeanintV = _WFintV/_hitNoV;
    _MeanintY = _WFintY/_hitNoY;

    _MeanRMS = _MeanRMS/_hitNo;
    _MeanRMSU = _MeanRMSU/_hitNoU;
    _MeanRMSV = _MeanRMSV/_hitNoV;
    _MeanRMSY = _MeanRMSY/_hitNoY;

    _MeanMult = _MeanMult/_hitNo;
    _MeanMultU = _MeanMultU/_hitNoU;
    _MeanMultV = _MeanMultV/_hitNoV;
    _MeanMultY = _MeanMultY/_hitNoY;

    // Calculate time spread
    if(_hitNo!=0){
      _TDCstd = TDCstd(TDCvec);
      _TDCiqr = TDCiqr(TDCvec,_hitNo);
      _Wirestd = TDCstd(Chvec);
      _Wireiqr = TDCiqr(Chvec,_hitNo);
    }
    // U plane
    if(_hitNoU!=0){
      _TDCstdU = TDCstd(UTDCvec);
      _TDCiqrU = TDCiqr(UTDCvec,_hitNoU);
      _WirestdU = TDCstd(UChvec);
      _WireiqrU = TDCiqr(UChvec,_hitNoU);
    }
    // V plane
    if(_hitNoV!=0){
      _TDCstdV = TDCstd(VTDCvec);
      _TDCiqrV = TDCiqr(VTDCvec,_hitNoV);
      _WirestdV = TDCstd(VChvec);
      _WireiqrV = TDCiqr(VChvec,_hitNoV);
    }
    // Y plane
    if(_hitNoY!=0){
      _TDCstdY = TDCstd(YTDCvec);
      _TDCiqrY = TDCiqr(YTDCvec,_hitNoY);
      _WirestdY = TDCstd(YChvec);
      _WireiqrY = TDCiqr(YChvec,_hitNoY);
    }

    // Fill TTree
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
