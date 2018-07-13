#ifndef CORR_SELECTOR_FACTORY_H
#define CORR_SELECTOR_FACTORY_H

#include "enums.h"
#include "mom_project.h"
#include "symm_data_multiplicities.h"


PiPiCorrelatorSelector* getSelector(const PiPiCorrSelector selector, const std::vector<threeMomentum> &pimom, const threeMomentum &p_tot, const std::string &data_dir,
				    const PiPiProjector proj_src, const PiPiProjector proj_snk, const PiPiMomAllowed allow){
  
  if(selector == PiPiCorrSelector::Basic){  
    return new PiPiCorrelatorBasicSelector(proj_src, proj_snk, allow, p_tot);
  }else if(selector == PiPiCorrSelector::SymmetrySubset){  
    PiPiSymmetrySubset corrs_avail(data_dir);
    return new PiPiMomSelectSymmetrySubset( corrs_avail.createSelector(pimom, p_tot, proj_src, proj_snk, allow) );
  }else{
    error_exit(std::cout << "getSelector unknown selector " << selector << std::endl);
  }
}

PiPiCorrelatorSelector* getSelector(const PiPiCorrSelector selector, const std::vector<threeMomentum> &pimom, const threeMomentum &p_tot, const std::string &data_dir,
				    PiPiProject* proj_src, PiPiProject* proj_snk, PiPiMomAllow* allow){
  
  if(selector == PiPiCorrSelector::Basic){  
    return new PiPiCorrelatorBasicSelector(proj_src, proj_snk, allow);
  }else if(selector == PiPiCorrSelector::SymmetrySubset){  
    PiPiSymmetrySubset corrs_avail(data_dir);
    return new PiPiMomSelectSymmetrySubset( corrs_avail.createSelector(pimom, p_tot, *proj_src, *proj_snk, *allow) );
  }else{
    error_exit(std::cout << "getSelector unknown selector " << selector << std::endl);
  }
}


#endif
