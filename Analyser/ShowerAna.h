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
#include <cmath>
#include "DataFormat/hit.h"
#include "DataFormat/mctruth.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string>
#include <fstream>

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

    int _hitNo;

    // Numbers of hits
    int _isHit;   

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
