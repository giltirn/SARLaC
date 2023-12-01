#ifndef _SARLAC_PRINTMANAGER_H_
#define _SARLAC_PRINTMANAGER_H_

#include<utils/macros.h>
#include<distribution/distribution_print/printers.h>
SARLAC_START_NAMESPACE

//A class that stores a singleton copy of the printer for a given type. The current printer is used in the stream operators for the distributions. The printer can be overridden at arbitrary time
template<typename DistributionType>
struct distributionPrint{

  //Call with no arguments to return the current printer for this type, or call with arguments to change the printer (with optional deletion of old one)
  static distributionPrinter<DistributionType>* printer(distributionPrinter<DistributionType>* change = NULL, bool delete_old = true){
    static distributionPrinter<DistributionType>* p = NULL;
    static bool initialized = false;
    if(!initialized){ p = new basicDistributionPrinter<DistributionType>; initialized = true; }

    if(change != NULL){ if(p!=NULL && delete_old) delete p; p = change; }
    return p;
  }

};

SARLAC_END_NAMESPACE
#endif
