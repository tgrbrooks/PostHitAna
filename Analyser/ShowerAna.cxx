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
    _t_ch->Branch("detShowerNo",&_detShowerNo,"detShowerNo/i");
    _t_ch->Branch("Evec",&Evec);
    _t_ch->Branch("Energy",&_Energy);
    _t_ch->Branch("ShowerStartEnd",&ShowerStartEnd);

    _t_ch->Branch("Meanamp",&_Meanamp,"Meanamp/F");
    _t_ch->Branch("MeanampU",&_MeanampU,"MeanampU/F");
    _t_ch->Branch("MeanampV",&_MeanampV,"MeanampV/F");
    _t_ch->Branch("MeanampY",&_MeanampY,"MeanampY/F");

    _t_ch->Branch("Meanint",&_Meanint,"Meanint/F");
    _t_ch->Branch("MeanintU",&_MeanintU,"MeanintU/F");
    _t_ch->Branch("MeanintV",&_MeanintV,"MeanintV/F");
    _t_ch->Branch("MeanintY",&_MeanintY,"MeanintY/F");

    _t_ch->Branch("LowDen",&_LowDen,"LowDen/i");
    _t_ch->Branch("LowDenU",&_LowDenU,"LowDenU/i");
    _t_ch->Branch("LowDenV",&_LowDenV,"LowDenV/i");
    _t_ch->Branch("LowDenY",&_LowDenY,"LowDenY/i");

    _t_ch->Branch("HiDen",&_HiDen,"HiDen/i");
    _t_ch->Branch("HiDenU",&_HiDenU,"HiDenU/i");
    _t_ch->Branch("HiDenV",&_HiDenV,"HiDenV/i");
    _t_ch->Branch("HiDenY",&_HiDenY,"HiDenY/i");

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
/*
    _t_ch->Branch("ShowerX",&_ShowerX,"ShowerX/F");
    _t_ch->Branch("ShowerY",&_ShowerY,"ShowerY/F");
    _t_ch->Branch("ShowerZ",&_ShowerZ,"ShowerZ/F");
    _t_ch->Branch("ShowerST",&_ShowerST,"ShowerST/F");
    _t_ch->Branch("ShowerET",&_ShowerET,"ShowerET/F");
    _t_ch->Branch("ShowerTheta",&_ShowerTheta,"ShowerTheta/F");
    _t_ch->Branch("ShowerPhi",&_ShowerPhi,"ShowerPhi/F");

    _t_ch->Branch("Type",&_Type,"Type/i");
    _t_ch->Branch("Length",&_Length,"Length/F");
    _t_ch->Branch("LengthU",&_LengthU,"LengthU/F");
    _t_ch->Branch("LengthV",&_LengthV,"LengthV/F");
    _t_ch->Branch("LengthY",&_LengthY,"LengthY/F");
    _t_ch->Branch("wLengthU",&_wLengthU,"wLengthU/F");
    _t_ch->Branch("wLengthV",&_wLengthV,"wLengthV/F");
    _t_ch->Branch("wLengthY",&_wLengthY,"wLengthY/F");
 
    _t_ch->Branch("ADCvec",&ADCvec);
    _t_ch->Branch("UADCvec",&UADCvec);
    _t_ch->Branch("VADCvec",&VADCvec);
    _t_ch->Branch("YADCvec",&YADCvec);
    _t_ch->Branch("TDCvec",&TDCvec);
    _t_ch->Branch("UTDCvec",&UTDCvec);
    _t_ch->Branch("VTDCvec",&VTDCvec);
    _t_ch->Branch("YTDCvec",&YTDCvec);  
    _t_ch->Branch("UChvec",&UChvec);
    _t_ch->Branch("VChvec",&VChvec);
    _t_ch->Branch("YChvec",&YChvec);
*/
    _t_ch->SetDirectory(0);

    return true;
  }
  
  // Called for every event
  bool ShowerAna::analyze(storage_manager* storage) {
   
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
	// Store neutrino energy
	_Energy = neut.Nu().Trajectory()[0].E();
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
    _detShowerNo = 0;
    _ShowerX = 0;
    _ShowerY = 0;
    _ShowerZ = 0;
    _ShowerST = 0;
    _ShowerET = 0;
    _ShowerTheta = 0;
    _ShowerPhi = 0;

    ShowerStartEnd.clear();
    Evec.clear();
    _Energy = 0;

    if (showerFlag){
      // Loop over number of showers? not sure this is how MCShower works
      for (size_t j=0; j < ev_mcs->size(); ++j){
        auto const& shower = (*ev_mcs).at(j);
        // Get energy, not sure if this is right energy
        _Energy = shower.Start().E();
        Evec.push_back(shower.Start().E());
        _ShowerX = shower.Start().X();
        _ShowerY = shower.Start().Y();
        _ShowerZ = shower.Start().Z();
        _ShowerST = shower.Start().T();
        _ShowerET = shower.End().T();
        float sX = shower.End().X()-_ShowerX;
        float sY = shower.End().Y()-_ShowerY;
        float sZ = shower.End().Z()-_ShowerZ;
        _ShowerPhi = atan(sY/sX);
        _ShowerTheta = acos(sZ/sqrt(sX*sX+sY*sY+sZ*sZ));
        float startX = _ShowerX/100;
        float startY = _ShowerY/100;
        float startZ = _ShowerZ/100;
        float endX = shower.End().X()/100;
        if((startX>0&&startX<2.5604)&&(endX>0&&endX<2.5604)){
          ShowerStartEnd.push_back(std::make_pair(startX, endX));
        }
        else ShowerStartEnd.push_back(std::make_pair(0, 0));
        // Count number of showers per event
        if((startX>0&&startX<2.5604)&&(startY>-1.2&&startY<1.2)&&(startZ>0&&startZ<10.5)&&_Energy>100){
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
/*    _Length = 0;_LengthU = 0;_LengthV = 0;_LengthY = 0; 
    _wLengthU = 0; _wLengthV = 0; _wLengthY = 0;
*/
    // Calculate the center of gravity
    float ampsumu = 0; float ampsumv = 0; float ampsumy = 0;
    // Loop over hits in each event
    for (size_t i=0; i < hitdata->size(); ++i){
      auto const& hit = (*hitdata).at(i);
      int chnum = hit.Channel();
      if(chnum<2400) ampsumu += hit.PeakAmplitude();
      if(chnum>=2400&&chnum<4800) ampsumv += hit.PeakAmplitude();
      if(chnum>=4800&&chnum<8256) ampsumy += hit.PeakAmplitude();
    }
    float tdcsumu = 0; float tdcsumv = 0; float tdcsumy = 0;
    float chsumu = 0; float chsumv = 0; float chsumy = 0;
    // Loop over hits in each event
    for (size_t i=0; i < hitdata->size(); ++i){
      auto const& hit = (*hitdata).at(i);
      // Get the channel number
      int chnum = hit.Channel();
      float amp = hit.PeakAmplitude();
      float tdc = hit.PeakTime();
      if(chnum<2400) {
        tdcsumu += tdc*amp/ampsumu;
        chsumu += chnum*amp/ampsumu;
      }
      if(chnum>=2400&&chnum<4800) {
        tdcsumv += tdc*amp/ampsumv;
        chsumv += chnum*amp/ampsumv;
      }
      if(chnum>=4800&&chnum<8256) {
        tdcsumy += tdc*amp/ampsumy;
        chsumy += chnum*amp/ampsumy;
      }
    }

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
      if(amp>10){
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
      if(amp<=20) _LowDen += 1;
      if(amp>20) _HiDen += 1;
      // Record times, amplitudes and channels of each hit
      TDCvec.push_back(tdc);
      ADCvec.push_back(amp);
      Chvec.push_back(chnum);

      // Do same for U plane CHECK THIS
      if(chnum<2400&&tdc>tdcsumu-320&&tdc<tdcsumu+320&&chnum>chsumu-80&&chnum<chsumu+80){ 
        _hitNoU += 1;
        _MeanampU += amp;
        _WFintU += wint;
        _MeanRMSU += rms;
        _MeanMultU += mult;
        UTDCvec.push_back(tdc);
        UADCvec.push_back(amp);
        UChvec.push_back(chnum);
        if(amp<=20)_LowDenU += 1;
        if(amp>20)_HiDenU += 1;
      }

      // Do same for V plane CHECK THIS
      if(chnum>=2400&&chnum<4800&&tdc>tdcsumv-320&&tdc<tdcsumv+320&&chnum>chsumv-80&&chnum<chsumv+80){
        _hitNoV += 1;
        _MeanampV += amp;
        _WFintV += wint;
        _MeanRMSV += rms;
        _MeanMultV += mult;
        VTDCvec.push_back(tdc);
        VADCvec.push_back(amp);
        VChvec.push_back(chnum);
        if(amp<=20)_LowDenV += 1;
        if(amp>20)_HiDenV += 1;
      }

      // Do same for Y plane CHECK THIS
      if(chnum>=4800&&chnum<8256&&tdc>tdcsumy-320&&tdc<tdcsumy+320&&chnum>chsumy-115&&chnum<chsumy+115){
        _hitNoY += 1;
        _MeanampY += amp;
        _WFintY += wint;
        _MeanRMSY += rms;
        _MeanMultY += mult;
        YTDCvec.push_back(tdc);
        YADCvec.push_back(amp);
        YChvec.push_back(chnum);
        if(amp<=20)_LowDenY += 1;
        if(amp>20)_HiDenY += 1;
      }
      }
    }

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
/*
    std::sort(TDCvec.begin(),TDCvec.end());
    std::sort(UTDCvec.begin(),UTDCvec.end());
    std::sort(VTDCvec.begin(),VTDCvec.end());
    std::sort(YTDCvec.begin(),YTDCvec.end());

    _Length = TDCvec[TDCvec.size()-5]-TDCvec[5];
    _LengthU = UTDCvec[UTDCvec.size()-5]-UTDCvec[5];
    _LengthV = VTDCvec[VTDCvec.size()-5]-VTDCvec[5];
    _LengthY = YTDCvec[YTDCvec.size()-5]-YTDCvec[5];

    std::sort(UChvec.begin(),UChvec.end());
    std::sort(VChvec.begin(),VChvec.end());
    std::sort(YChvec.begin(),YChvec.end());

    _wLengthU = UChvec[UChvec.size()-5]-UChvec[5];
    _wLengthV = VChvec[VChvec.size()-5]-VChvec[5];
    _wLengthY = YChvec[YChvec.size()-5]-YChvec[5];
*/
    // Fill TTree
    if(_hitNoY==0||_hitNoU==0||_hitNoV==0) return false; 
///////////////////////////////////////////////////////////////////////////////////////////
    if(_showerNo==0) _t_ch->Fill();
////////////////////////////////////////////////////////////////////////////////////////////
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
