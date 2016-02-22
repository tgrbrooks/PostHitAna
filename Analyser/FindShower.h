/**
 * \file FindShower.h
 *
 * \ingroup Analyser
 * 
 * \brief Class def header for a class FindShower
 *
 * @author tbrooks
 */

/** \addtogroup Analyser

    @{*/

#ifndef LARLITE_FINDSHOWER_H
#define LARLITE_FINDSHOWER_H

#include "Analysis/ana_base.h"
#include <iostream>
#include <vector>
#include <utility>
#include <cmath>
#include "DataFormat/hit.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <string>
#include <fstream>

namespace larlite {
  /**
     \class FindShower
     User custom analysis class made by SHELL_USER_NAME
   */
  class FindShower : public ana_base{
  
  public:

    /// Default constructor
    FindShower(){ _name="FindShower"; _fout=0;}

    /// Default destructor
    virtual ~FindShower(){}

    /** IMPLEMENT in FindShower.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in FindShower.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in FindShower.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    int _evtN;
    std::vector<float> TDCvec;
    std::vector<float> Chvec;
    // Vector of Tracks
    std::vector<std::vector<std::pair<float,float>>> Tracks;
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
