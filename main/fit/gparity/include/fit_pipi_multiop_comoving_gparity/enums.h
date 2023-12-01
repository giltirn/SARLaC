#ifndef _FIT_PIPI_COMOVING_GPARITY_ENUMS_H
#define _FIT_PIPI_COMOVING_GPARITY_ENUMS_H

GENERATE_ENUM_AND_PARSER(Operator, (PiPiComoveGnd)(PiPiComoveExc1)(PiPiComoveExc2) );
void write(HDF5writer &wr, const Operator op, const std::string &nm){
  int p = (int)op;
  _force_external_lookup<int>::fwrite(wr, p, nm);
}
void read(HDF5reader &rd, Operator &op, const std::string &nm){
  int p;
  _force_external_lookup<int>::fread(rd, p, nm);
  op = (Operator)p;
}

#endif
