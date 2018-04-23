#ifndef PIPI_MOM_PROJECT_FACTORY_H_
#define PIPI_MOM_PROJECT_FACTORY_H_

PiPiProject* getProjector(const PiPiProjector p){
  switch(p){
  case A1:
    return (PiPiProject*)(new PiPiProjectA1);
  case Avg4:
    return (PiPiProject*)(new PiPiProjectAvg4);
  case Avg2:
    return (PiPiProject*)(new PiPiProjectAvg2);
  case Solo:
    return (PiPiProject*)(new PiPiProjectSolo);
  default:
    error_exit(std::cout << "getProjector unknown projector " << p << std::endl);
  }
}

PiPiMomAllow* getMomPairFilter(const PiPiMomAllowed p){
  switch(p){
  case All:
    return (PiPiMomAllow*)(new PiPiMomAllowAll);
  case Orig64:
    return (PiPiMomAllow*)(new PiPiMomAllowOrig64);
  case ParityAxisPermSymmReduced:
    return (PiPiMomAllow*)(new PiPiMomAllowParityAxisPermSymmReduced);
  case AuxDiagSymmReduced:
    return (PiPiMomAllow*)(new PiPiMomAllowAuxDiagSymmReduced);
  case AuxDiagParityAxisPermSymmReduced:
    return (PiPiMomAllow*)(new PiPiMomAllowAuxDiagParityAxisPermSymmReduced);
  default:
    error_exit(std::cout << "getMomPairFilter unknown filter " << p << std::endl);
  }
}

#endif
