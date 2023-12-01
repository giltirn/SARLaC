#ifndef _FIT_PIPI_GND_EXC_SIGMA_GPARITY_ENUMS_H_
#define _FIT_PIPI_GND_EXC_SIGMA_GPARITY_ENUMS_H_

GENERATE_ENUM_AND_PARSER(Operator, (PiPiGnd)(PiPiExc)(Sigma) );
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
