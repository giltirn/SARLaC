#ifndef _FIT_PIPI_COMOVING_RAW_DATA_H
#define _FIT_PIPI_COMOVING_RAW_DATA_H

class RawData{
  NumericSquareMatrix<bubbleDataAllMomenta> pipi_bubble;
  NumericSquareMatrix<rawCorrelationFunction> correlators;
  std::set<std::pair<Operator,Operator> > contains;
public:
  //Real self-contractions
  inline bubbleDataAllMomenta & PiPiBubble(const Operator srcop, const Operator snkop){  return pipi_bubble((int)srcop, (int)snkop);  }
  inline const bubbleDataAllMomenta & PiPiBubble(const Operator srcop, const Operator snkop) const{  return pipi_bubble((int)srcop, (int)snkop);  }

  //Correlator data
  inline rawCorrelationFunction & correlator(const Operator srcop, const Operator snkop){ return correlators((int)srcop, (int)snkop); }
  inline const rawCorrelationFunction & correlator(const Operator srcop, const Operator snkop) const{ return correlators((int)srcop, (int)snkop); }
  
  RawData(): pipi_bubble(3), correlators(3){}

  inline bool doOp(const Operator op, const std::vector<Operator> &incl_ops) const{ 
    return std::find(incl_ops.begin(), incl_ops.end(), op) != incl_ops.end();
  }

  bool haveData(const Operator opa, const Operator opb) const{ 
    return contains.find({opa,opb}) != contains.end();
  }

  void write(HDF5writer &wr, const std::string &nm) const{
    wr.enter(nm);
    CPSfit::write(wr, pipi_bubble, "pipi_bubble");
    CPSfit::write(wr, correlators, "correlators");
    CPSfit::write(wr, contains, "contains");
    wr.leave();
  }
  void read(HDF5reader &rd, const std::string &nm){
    rd.enter(nm);
    CPSfit::read(rd, pipi_bubble, "pipi_bubble");
    CPSfit::read(rd, correlators, "correlators");
    CPSfit::read(rd, contains, "contains");
    rd.leave();
  }

  void read(const int Lt, const std::string &data_dir, const int traj_start, const int traj_inc, const int traj_lessthan,
	    const std::string &pipi_fig_file_fmt, const std::string &pipi_bubble_file_fmt, const int tsep_pipi, const int tstep_pipi,
	    const threeMomentum &p_tot, const std::vector<Operator> &incl_ops, bool filemap_allow_ptot_parity = false){

    static const std::vector< std::pair<Operator,Operator> > op_pairs = { {Operator::PiPiComoveGnd, Operator::PiPiComoveGnd},
									  {Operator::PiPiComoveGnd, Operator::PiPiComoveExc1},
									  {Operator::PiPiComoveGnd, Operator::PiPiComoveExc2},
									  {Operator::PiPiComoveExc1, Operator::PiPiComoveExc1},
									  {Operator::PiPiComoveExc1, Operator::PiPiComoveExc2},
									  {Operator::PiPiComoveExc2, Operator::PiPiComoveExc2} };

    static const std::map<Operator, PiPiProjector> op_proj = { {Operator::PiPiComoveGnd, PiPiProjector::MovingSwaveGround},
							       {Operator::PiPiComoveExc1, PiPiProjector::MovingSwaveExc1},
							       {Operator::PiPiComoveExc2, PiPiProjector::MovingSwaveExc2} };
    
    for(int i=0; i<op_pairs.size(); i++){
      auto o1 = op_pairs[i].first; auto o2 = op_pairs[i].second;

      if(doOp(o1, incl_ops) && doOp(o2, incl_ops)){
	readPiPi2pt(correlator(o1,o2), PiPiBubble(o1,o2), data_dir, 
		    pipi_fig_file_fmt, pipi_bubble_file_fmt, tsep_pipi, tstep_pipi, Lt, 
		    traj_start, traj_inc, traj_lessthan, 
		    p_tot, op_proj.find(o1)->second, op_proj.find(o2)->second,
		    0, filemap_allow_ptot_parity);
	contains.insert({o1,o2});
      }
    }
  }
};

inline void read(HDF5reader &rd, RawData &data, const std::string &nm){
  data.read(rd,nm);
}
inline void write(HDF5writer &wr, const RawData &data, const std::string &nm){
  data.write(wr,nm);
}


#endif
