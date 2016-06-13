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
    sig_tree->Branch("detShowerNo",&_detShowerNo,"detShowerNo/i");
    sig_tree->Branch("Type",&_Type,"Type/i");
    sig_tree->Branch("hitNoU",&_hitNoU,"hitNoU/i");
    sig_tree->Branch("hitNoV",&_hitNoV,"hitNoV/i");
    sig_tree->Branch("hitNoY",&_hitNoY,"hitNoY/i");
    sig_tree->Branch("TDCstdU",&_TDCstdU,"TDCstdU/F");
    sig_tree->Branch("TDCstdV",&_TDCstdV,"TDCstdV/F");
    sig_tree->Branch("TDCstdY",&_TDCstdY,"TDCstdY/F");
    sig_tree->Branch("TDCiqrU",&_TDCiqrU,"TDCiqrU/F");
    sig_tree->Branch("TDCiqrV",&_TDCiqrV,"TDCiqrV/F");
    sig_tree->Branch("TDCiqrY",&_TDCiqrY,"TDCiqrY/F");
    sig_tree->Branch("ADCampU",&_ADCampU,"ADCampU/F");
    sig_tree->Branch("ADCampV",&_ADCampV,"ADCampV/F");
    sig_tree->Branch("ADCampY",&_ADCampY,"ADCampY/F");
    sig_tree->Branch("WFintU",&_WFintU,"WFintU/F");
    sig_tree->Branch("WFintV",&_WFintV,"WFintV/F");
    sig_tree->Branch("WFintY",&_WFintY,"WFintY/F");
    sig_tree->Branch("MeanampU",&_MeanampU,"MeanampU/F");
    sig_tree->Branch("MeanampV",&_MeanampV,"MeanampV/F");
    sig_tree->Branch("MeanampY",&_MeanampY,"MeanampY/F");
    sig_tree->Branch("MeanintU",&_MeanintU,"MeanintU/F");
    sig_tree->Branch("MeanintV",&_MeanintV,"MeanintV/F");
    sig_tree->Branch("MeanintY",&_MeanintY,"MeanintY/F");
    sig_tree->Branch("LowDenU",&_LowDenU,"LowDenU/i");
    sig_tree->Branch("LowDenV",&_LowDenV,"LowDenV/i");
    sig_tree->Branch("LowDenY",&_LowDenY,"LowDenY/i");
    sig_tree->Branch("HiDenU",&_HiDenU,"HiDenU/i");
    sig_tree->Branch("HiDenV",&_HiDenV,"HiDenV/i");
    sig_tree->Branch("HiDenY",&_HiDenY,"HiDenY/i");
    sig_tree->Branch("MeanRMSU",&_MeanRMSU,"MeanRMSU/F");
    sig_tree->Branch("MeanRMSV",&_MeanRMSV,"MeanRMSV/F");
    sig_tree->Branch("MeanRMSY",&_MeanRMSY,"MeanRMSY/F");
    sig_tree->Branch("MeanMultU",&_MeanMultU,"MeanMultU/F");
    sig_tree->Branch("MeanMultV",&_MeanMultV,"MeanMultV/F");
    sig_tree->Branch("MeanMultY",&_MeanMultY,"MeanMultY/F");
    sig_tree->Branch("WirestdU",&_WirestdU,"WirestdU/F");
    sig_tree->Branch("WirestdV",&_WirestdV,"WirestdV/F");
    sig_tree->Branch("WirestdY",&_WirestdY,"WirestdY/F");
    sig_tree->Branch("WireiqrU",&_WireiqrU,"WireiqrU/F");
    sig_tree->Branch("WireiqrV",&_WireiqrV,"WireiqrV/F");
    sig_tree->Branch("WireiqrY",&_WireiqrY,"WireiqrY/F");
    sig_tree->Branch("NumMultU",&_NumMultU,"NumMultU/i");
    sig_tree->Branch("NumMultV",&_NumMultV,"NumMultV/i");
    sig_tree->Branch("NumMultY",&_NumMultY,"NumMultY/i");
    sig_tree->SetDirectory(0);

    // Initialise background TTree for Y variables and 960x346 boxes
    bkg_tree = sig_tree->CloneTree();
    bkg_tree->SetDirectory(0);
    bkg_tree->SetName("bkg_tree");

    return true;
  }
  
  bool GenerateSigBkg::analyze(storage_manager* storage) {


    bool truthflag = true;

    auto ev_mct = storage->get_data<event_mctruth>("generator");

    // Display error if rawdigit data not present
    if ( (!ev_mct) || (!ev_mct->size())){
    	print (msg::kERROR,__FUNCTION__,"MCTruth data product not found!");
    	truthflag = false;
    }
    if(truthflag){
        auto const& mct = (*ev_mct).at(0);
        // Get neutrino inf
        auto const& neut = (mct.GetNeutrino());
        // Get charged current or neutral current
        int ccnc = neut.CCNC();
	// Get interaction type
	int intMode = neut.Mode();
	if(intMode==0&&ccnc==0){_Type = 0;}//CCQE
	else if(intMode==0&&ccnc==1){_Type = 1;}//NCQE
	else if(intMode==1&&ccnc==0){_Type = 2;}//CCRE
	else if(intMode==1&&ccnc==1){_Type = 3;}//NCRE
	else if(intMode==2&&ccnc==0){_Type = 4;}//CCDIS
	else if(intMode==2&&ccnc==1){_Type = 5;}//NCDIS
	else if(intMode==3&&ccnc==0){_Type = 6;}//CCCO
	else if(intMode==3&&ccnc==1){_Type = 7;}//NCCO
        else{_Type = 8;}
    } else _Type = 8;
  
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
        float endX = shower.End().X()/100;
        float endY = shower.End().Y()/100;
        float endZ = shower.End().Z()/100;
        // Count number of showers per event
        if(((startX>0&&startX<2.5604)&&(startY>-1.2&&startY<1.2)&&(startZ>0&&startZ<10.5)&&Energy>100&&Energy<500)
           ||((endX>0&&endX<2.5604)&&(endY>-1.2&&endY<1.2)&&(endZ>0&&endZ<10.5)&&Energy>100&&Energy<500)){
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
    _hitNoU = 0; _hitNoV = 0; _hitNoY = 0;
    // Initialize integration counters for event
    _WFintU = 0; _WFintV = 0; _WFintY = 0;
    // Clear TDC vectors
    UTDCvec.clear(); VTDCvec.clear(); YTDCvec.clear();
    // Clear ADC vectors
    UADCvec.clear(); VADCvec.clear(); YADCvec.clear();
    // Clear Channel number vectors
    UChvec.clear(); VChvec.clear(); YChvec.clear();
    // Initialize TDC standard deviation counters
    _TDCstdU = 0; _TDCstdV = 0; _TDCstdY = 0;
    // Initialize TDC interquartile range counters
    _TDCiqrU = 0; _TDCiqrV = 0; _TDCiqrY = 0;
    // Initialize max (or total?) ADC amplitude counters
    _ADCampU = 0; _ADCampV = 0; _ADCampY = 0;
    // Initialize mean adc amplitude counters
    _MeanampU = 0; _MeanampV = 0; _MeanampY = 0;
    // Initailize mean integrated waveform counters
    _MeanintU = 0; _MeanintV = 0; _MeanintY = 0;
    // Initialize hit density between 5-35 ADCs
    _LowDenU = 0; _LowDenV = 0; _LowDenY = 0;
    // Initialize hit density > 35
    _HiDenU = 0; _HiDenV = 0; _HiDenY = 0;
    // Initialize mean RMS counters
    _MeanRMSU = 0; _MeanRMSV = 0; _MeanRMSY = 0;
    // Initialize mean multiplicity counters
    _MeanMultU = 0; _MeanMultV = 0; _MeanMultY = 0;
    // Initialize Wire standard deviation counter
    _WirestdU = 0; _WirestdV = 0; _WirestdY = 0;
    // Initialize Wire interquartile range counter
    _WireiqrU = 0; _WireiqrV = 0; _WireiqrY = 0;
    _NumMultU = 0; _NumMultV = 0; _NumMultY = 0;

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
        if(amp<=26)_LowDenU += 1;
        if(amp>26)_HiDenU += 1;
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
        if(amp<=26)_LowDenV += 1;
        if(amp>26)_HiDenV += 1;
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
        if(amp<=26)_LowDenY += 1;
        if(amp>26)_HiDenY += 1;
      }
      }
    }

    _ADCampU = gMaxVal(UADCvec);
    _MeanampU = _MeanampU/_hitNoU;
    _MeanintU = _WFintU/_hitNoU;
    _MeanRMSU = _MeanRMSU/_hitNoU;
    _MeanMultU = _MeanMultU/_hitNoU;

    _ADCampV = gMaxVal(VADCvec);
    _MeanampV = _MeanampV/_hitNoV;
    _MeanintV = _WFintV/_hitNoV;
    _MeanRMSV = _MeanRMSV/_hitNoV;
    _MeanMultV = _MeanMultV/_hitNoV;

    _ADCampY = gMaxVal(YADCvec);
    _MeanampY = _MeanampY/_hitNoY;
    _MeanintY = _WFintY/_hitNoY;
    _MeanRMSY = _MeanRMSY/_hitNoY;
    _MeanMultY = _MeanMultY/_hitNoY;

    // U plane
    if(_hitNoU!=0){
      _TDCstdU = gTDCstd(UTDCvec);
      _TDCiqrU = gTDCiqr(UTDCvec,_hitNoU);
      _WirestdU = gTDCstd(UChvec);
      _WireiqrU = gTDCiqr(UChvec,_hitNoU);
    }

    if(_hitNoV!=0){
      _TDCstdV = gTDCstd(VTDCvec);
      _TDCiqrV = gTDCiqr(VTDCvec,_hitNoV);
      _WirestdV = gTDCstd(VChvec);
      _WireiqrV = gTDCiqr(VChvec,_hitNoV);
    }

    if(_hitNoY!=0){
      _TDCstdY = gTDCstd(YTDCvec);
      _TDCiqrY = gTDCiqr(YTDCvec,_hitNoY);
      _WirestdY = gTDCstd(YChvec);
      _WireiqrY = gTDCiqr(YChvec,_hitNoY);
    }

    int current_num_hits, next_num_hits;

    for(int k(0); k<2400; k++){
      current_num_hits = 0;
      next_num_hits = 0;
      for(int l(0); l<UChvec.size(); l++){
        if(UChvec[l]==k) current_num_hits++;
        if(UChvec[l]==k+1) next_num_hits++;
      }
      if(std::abs(current_num_hits-next_num_hits)>3) _NumMultU++;
    }

    for(int k(2400); k<4800; k++){
      current_num_hits = 0;
      next_num_hits = 0;
      for(int l(0); l<VChvec.size(); l++){
        if(VChvec[l]==k) current_num_hits++;
        if(VChvec[l]==k+1) next_num_hits++;
      }
      if(std::abs(current_num_hits-next_num_hits)>3) _NumMultV++;
    }

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
    if(_hitNoU==0||_ADCampU>2000) return false; 
    if(_hitNoV==0||_ADCampV>2000) return false;
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
