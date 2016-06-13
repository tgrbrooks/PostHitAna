/**
 * \file GenerateSigBkg.h
 *
 * \ingroup Analyser
 * 
 * \brief Class def header for a class GenerateSigBkg
 *
 * @author tom
 */

/** \addtogroup Analyser

    @{*/

#ifndef LARLITE_GENERATESIGBKG_H
#define LARLITE_GENERATESIGBKG_H

#include "Analysis/ana_base.h"
#include <iostream>
#include <vector>
#include <utility>
//#include <cmath>
#include "DataFormat/hit.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mctruth.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <string>
#include <fstream>
#include <math.h>

namespace larlite {
  /**
     \class GenerateSigBkg
     User custom analysis class made by SHELL_USER_NAME
   */
  class GenerateSigBkg : public ana_base{
  
  public:

    /// Default constructor
    GenerateSigBkg(){ _name="GenerateSigBkg"; _fout=0;}

    /// Default destructor
    virtual ~GenerateSigBkg(){}

    /** IMPLEMENT in GenerateSigBkg.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in GenerateSigBkg.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in GenerateSigBkg.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    int _evtN;

    int _showerNo, _detShowerNo;

    int _Type;

    int _hitNoU, _hitNoV, _hitNoY;

    float _TDCstdU, _TDCstdV, _TDCstdY;

    float _TDCiqrU, _TDCiqrV, _TDCiqrY;

    float _ADCampU, _ADCampV, _ADCampY;

    float _MeanampU, _MeanampV, _MeanampY;

    float _WFintU, _WFintV, _WFintY;

    float _MeanintU, _MeanintV, _MeanintY;

    int _LowDenU, _LowDenV, _LowDenY;

    int _HiDenU, _HiDenV, _HiDenY;

    float _MeanRMSU, _MeanRMSV, _MeanRMSY;

    float _MeanMultU, _MeanMultV, _MeanMultY;

    float _WirestdU, _WirestdV, _WirestdY;

    float _WireiqrU, _WireiqrV, _WireiqrY;

    int _NumMultU, _NumMultV, _NumMultY;

    std::vector<float> TDCvec, UTDCvec, VTDCvec, YTDCvec;
   
    std::vector<float> ADCvec, UADCvec, VADCvec, YADCvec;

    std::vector<float> Chvec, UChvec, VChvec, YChvec;

    TTree* sig_tree; TTree* bkg_tree; 
    
  };
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
