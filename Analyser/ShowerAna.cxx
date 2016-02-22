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
int TDCiqr(std::vector<float> TDCvec, int hitNo) {

    std::sort(TDCvec.begin(),TDCvec.end());
    int quart =(int) hitNo/4;
    int iqrTDC = (TDCvec[3*quart] - TDCvec[quart]);
    return iqrTDC;

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
    _t_ch->Branch("TDCvec",&TDCvec);
    _t_ch->Branch("UTDCvec",&UTDCvec);
    _t_ch->Branch("VTDCvec",&VTDCvec);
    _t_ch->Branch("YTDCvec",&YTDCvec);
    _t_ch->Branch("ADCvec",&ADCvec);
    _t_ch->Branch("UADCvec",&UADCvec);
    _t_ch->Branch("VADCvec",&VADCvec);
    _t_ch->Branch("YADCvec",&YADCvec);
    _t_ch->Branch("UChvec",&UChvec);
    _t_ch->Branch("VChvec",&VChvec);
    _t_ch->Branch("YChvec",&YChvec);
    _t_ch->Branch("ShowerStartEnd",&ShowerStartEnd);
/*    _t_ch->Branch("Meanamp",&_Meanamp,"Meanamp/F");
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
    _t_ch->Branch("StartEnd",&StartEnd);
    _t_ch->Branch("StartEndU",&StartEndU);
    _t_ch->Branch("StartEndV",&StartEndV);
    _t_ch->Branch("StartEndY",&StartEndY);
    _t_ch->Branch("ChStartEnd",&ChStartEnd);
    _t_ch->Branch("ChStartEndU",&ChStartEndU);
    _t_ch->Branch("ChStartEndV",&ChStartEndV);
    _t_ch->Branch("ChStartEndY",&ChStartEndY);
*/
    _t_ch->SetDirectory(0);

    return true;
  }
  
  // Called for every event
  bool ShowerAna::analyze(storage_manager* storage) {
   
    float DriftVel = 1600;

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
        float startTime = shower.Start().T();
        float endTime = shower.End().T();
        float startX = shower.Start().X()/100;
        float endX = shower.End().X()/100;
        if((startX>0&&startX<2.5604)||(endX>0&&endX<2.5604)){
          ShowerStartEnd.push_back(std::make_pair(startTime + (startX/(DriftVel*0.0000005)), endTime + (endX/(DriftVel*0.0000005))));
          std::cout<<"Start T "<<startTime<<" End T "<<endTime<<std::endl;
          std::cout<<"Start X "<<startX<<" End X "<<endX<<std::endl;
          std::cout<<"Start "<<startTime + startX*DriftVel<<" End "<<endTime + endX*DriftVel<<std::endl;
        }
        else ShowerStartEnd.push_back(std::make_pair(0, 0));
        // Count number of showers per event
        _showerNo += 1;
      }
    }


    // Use storage to get larlite::hit object (gaushit, cccluster, pandoraCosmicKHitRemoval)
    auto hitdata = storage->get_data<event_hit>("pandoraCosmicKHitRemoval");
    // Display error if hit data not present
    if ( (!hitdata) || (!hitdata->size())){
    	print (msg::kERROR,__FUNCTION__,"Hit data product not found!");
    	return false;
    }

    // Initialise hit counters for event
    _hitNo = 0;
    _hitNoU = 0;
    _hitNoV = 0;
    _hitNoY = 0;

    // Initialize integration counters for event
    _WFint = 0;
    _WFintU = 0;
    _WFintV = 0;
    _WFintY = 0;

    // Clear TDC vectors
    TDCvec.clear();
    UTDCvec.clear();
    VTDCvec.clear();
    YTDCvec.clear();

    ADCvec.clear();
    UADCvec.clear();
    VADCvec.clear();
    YADCvec.clear();

    UChvec.clear();
    VChvec.clear();
    YChvec.clear();

    // Initialize TDC standard deviation counters
    _TDCstd = 0;
    _TDCstdU = 0;
    _TDCstdV = 0;
    _TDCstdY = 0;

    // Initialize TDC interquartile range counters
    _TDCiqr = 0;
    _TDCiqrU = 0;
    _TDCiqrV = 0;
    _TDCiqrY = 0;

    // Initialize ADC amplitude counters
    _ADCamp = 0;
    _ADCampU = 0;
    _ADCampV = 0;
    _ADCampY = 0;

    // Loop over hits in each event
    for (size_t i=0; i < hitdata->size(); ++i){
      auto const& hit = (*hitdata).at(i);
      // Get the channel number
      int chnum = hit.Channel();

      // Count number of hits
      _hitNo += 1;
      // Add up the hit amplitudes
      _ADCamp += hit.PeakAmplitude();
      // Add up the integrated waveforms
      _WFint += hit.Integral();
      // Record times of each hit
      TDCvec.push_back(hit.PeakTime());
      ADCvec.push_back(hit.PeakAmplitude());

      // Do same for U plane CHECK THIS
      if(chnum<2400){ 
        _hitNoU += 1;
        _ADCampU += hit.PeakAmplitude();
        _WFintU += hit.Integral();
        UTDCvec.push_back(hit.PeakTime());
        UADCvec.push_back(hit.PeakAmplitude());
        UChvec.push_back(chnum);
      }

      // Do same for V plane CHECK THIS
      if(chnum>=2400&&chnum<4800){
        _hitNoV += 1;
        _ADCampV += hit.PeakAmplitude();
        _WFintV += hit.Integral();
        VTDCvec.push_back(hit.PeakTime());
        VADCvec.push_back(hit.PeakAmplitude());
        VChvec.push_back(chnum-2400);
      }

      // Do same for Y plane CHECK THIS
      if(chnum>=4800&&chnum<8256){
        _hitNoY += 1;
        _ADCampY += hit.PeakAmplitude();
        _WFintY += hit.Integral();
        YTDCvec.push_back(hit.PeakTime());
        YADCvec.push_back(hit.PeakAmplitude());
        YChvec.push_back(chnum-4800);
      }

    }

    // Average the ADC amplitudes
    _ADCamp = _ADCamp/_hitNo;
    _ADCampU = _ADCampU/_hitNoU;
    _ADCampV = _ADCampV/_hitNoV;
    _ADCampY = _ADCampY/_hitNoY;

    // Calculate time spread
    if(_hitNo!=0){
      _TDCstd = TDCstd(TDCvec);
      _TDCiqr = TDCiqr(TDCvec,_hitNo);
    }
    // U plane
    if(_hitNoU!=0){
      _TDCstdU = TDCstd(UTDCvec);
      _TDCiqrU = TDCiqr(UTDCvec,_hitNoU);
    }
    // V plane
    if(_hitNoV!=0){
      _TDCstdV = TDCstd(VTDCvec);
      _TDCiqrV = TDCiqr(VTDCvec,_hitNoV);
    }
    // Y plane
    if(_hitNoY!=0){
      _TDCstdY = TDCstd(YTDCvec);
      _TDCiqrY = TDCiqr(YTDCvec,_hitNoY);
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
