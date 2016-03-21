/**
 * \file ShowerAna.h
 *
 * \ingroup Analyser
 * 
 * \brief Class def header for a class ShowerAna
 *
 * @author tbrooks
 */

/** \addtogroup Analyser

    @{*/

#ifndef LARLITE_SHOWERANA_H
#define LARLITE_SHOWERANA_H

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
     \class ShowerAna
     User custom analysis class made by SHELL_USER_NAME
   */
  class ShowerAna : public ana_base{
  
  public:

    /// Default constructor
    ShowerAna(){ _name="ShowerAna"; _fout=0;}

    /// Default destructor
    virtual ~ShowerAna(){}

    /** IMPLEMENT in ShowerAna.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in ShowerAna.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in ShowerAna.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    int _evtN;

    int _showerNo, _detShowerNo;

    int _hitNo, _hitNoU, _hitNoV, _hitNoY, _Type;

    float _TDCstd, _TDCstdU, _TDCstdV, _TDCstdY;

    float _TDCiqr, _TDCiqrU, _TDCiqrV, _TDCiqrY;

    float _ADCamp, _ADCampU, _ADCampV, _ADCampY;

    float _Meanamp, _MeanampU, _MeanampV, _MeanampY;

    float _WFint, _WFintU, _WFintV, _WFintY;

    float _Meanint, _MeanintU, _MeanintV, _MeanintY;

    int _LowDen, _LowDenU, _LowDenV, _LowDenY;

    int _HiDen, _HiDenU, _HiDenV, _HiDenY;

    float _MeanRMS, _MeanRMSU, _MeanRMSV, _MeanRMSY;

    float _MeanMult, _MeanMultU, _MeanMultV, _MeanMultY;

    float _Wirestd, _WirestdU, _WirestdV, _WirestdY;

    float _Wireiqr, _WireiqrU, _WireiqrV, _WireiqrY;

    float _ShowerX, _ShowerY, _ShowerZ, _ShowerTheta, _ShowerPhi;

    float _Length, _LengthU, _LengthV, _LengthY;

    float _wLengthU, _wLengthV, _wLengthY;

    std::vector<float> TDCvec, UTDCvec, VTDCvec, YTDCvec;
   
    std::vector<float> ADCvec, UADCvec, VADCvec, YADCvec;

    std::vector<float> Chvec, UChvec, VChvec, YChvec;

    std::vector<std::pair<float,float>> ShowerStartEnd;

    float _Energy;

    std::vector<float> Evec;

    TTree* _t_ch;
    
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
