#ifndef _FIT_PIPI_KTOPIPI_SIM_GPARITY_FITFUNC_H_
#define _FIT_PIPI_KTOPIPI_SIM_GPARITY_FITFUNC_H_

//Fit in the 7 basis
struct TwoPointThreePointSimFitParams{
  double AK;
  double mK;
  double Apipi;
  double Epipi;
  double Cpipi; //constant term
  std::vector<double> M;

  enum {nprelim = 5};
  
  TwoPointThreePointSimFitParams(): AK(1), mK(0.5), Apipi(0.16), Epipi(0.36), Cpipi(0), M(7){
    M[0] = -0.003; M[1] = -0.003; M[2] = 0.005; M[3] = -0.01; M[4] = -0.02; M[5] = 0.02; M[6] = 0.07;
  }
  TwoPointThreePointSimFitParams(const double _AK, const double _mK, const double _Apipi, const double _Epipi, const double _Cpipi, const double _M): AK(_AK), mK(_mK), Apipi(_Apipi), Epipi(_Epipi), Cpipi(_Cpipi), M(7,_M){}
  TwoPointThreePointSimFitParams(const TwoPointThreePointSimFitParams &r) = default;
  TwoPointThreePointSimFitParams(TwoPointThreePointSimFitParams &&r) = default;

  TwoPointThreePointSimFitParams& operator=(const TwoPointThreePointSimFitParams &r) = default;
  TwoPointThreePointSimFitParams& operator=(TwoPointThreePointSimFitParams &&r) = default;
    
  inline int size() const{ return 12; }

  inline std::string print() const{
    std::ostringstream os; os << "AK=" << AK << " mK=" << mK << " Apipi=" << Apipi << " Epipi=" << Epipi << " Cpipi=" << Cpipi;
    for(int i=0;i<7;i++) os << " M[" << i << "]=" << M[i];
    return os.str();
  }
  inline void zero(){
    AK=mK=Apipi=Epipi=Cpipi=0.;
    for(int i=0;i<7;i++) M[i] = 0.;
  }
  inline double &operator()(const int i){
    if(i<5){
      switch(i){
      case 0:
	return AK;
      case 1:
	return mK;
      case 2:
	return Apipi;
      case 3:
	return Epipi;
      case 4:
	return Cpipi;
      }
    }else return M[i-5];
  }
  inline const double &operator()(const int i) const{ return const_cast<const double &>( const_cast<TwoPointThreePointSimFitParams*>(this)->operator()(i) ); }

  typedef TwoPointThreePointSimFitParams ET_tag;
  
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,TwoPointThreePointSimFitParams >::value, int>::type = 0>
  TwoPointThreePointSimFitParams(U&& expr): M(7){
    this->AK = expr[0];
    this->mK = expr[1];
    this->Apipi = expr[2];
    this->Epipi = expr[3];
    this->Cpipi = expr[4];
    for(int i=5;i<12;i++) this->M[i-5] = expr[i];
  }
};
template<>
struct CPSfit::getElem<TwoPointThreePointSimFitParams>{
  static inline auto elem(const TwoPointThreePointSimFitParams &v, const int i)->decltype(v(i)){ return v(i); }
  static inline int common_properties(const TwoPointThreePointSimFitParams &v){ return 0; }
};

template<>
struct CPSfit::printStats<jackknifeDistribution<TwoPointThreePointSimFitParams> >{
  inline static std::string centralValue(const jackknifeDistribution<TwoPointThreePointSimFitParams> &d){
    auto best = d.best();
    std::ostringstream os; os << "(" << best.AK << ", " << best.mK << ", " << best.Apipi << ", " << best.Epipi << ", " << best.Cpipi;
    for(int i=0;i<7;i++) os << ", " << best.M[i];
    os << ")";
    return os.str();
  }
  inline static std::string error(const jackknifeDistribution<TwoPointThreePointSimFitParams> &d){
    auto err = d.standardError();
    std::ostringstream os; os << "(" << err.AK << ", " << err.mK << ", " << err.Apipi << ", " << err.Epipi << ", " << err.Cpipi;
    for(int i=0;i<7;i++) os << ", " << err.M[i];
    os << ")";
    return os.str();
  }
};

GENERATE_PARSER(TwoPointThreePointSimFitParams, (double, AK)(double, mK)(double, Apipi)(double, Epipi)(double, Cpipi)(std::vector<double>, M) )


class FitTwoPointThreePointSim{
  const double Lt;
  const double tsep_pipi;
  const double Ascale;
  const double Cscale;

public:
  typedef TwoPointThreePointSimFitParams Params;
  typedef double ValueType;
  typedef Params ParameterType;
  typedef Params ValueDerivativeType; //derivative wrt parameters
  typedef amplitudeDataCoordSim GeneralizedCoordinate; //use idx = -1 for pipi data

  FitTwoPointThreePointSim(const int _Lt, const double _tsep_pipi, const double _Ascale=1e13, const double _Cscale=1e13): Lt(_Lt),tsep_pipi(_tsep_pipi),Ascale(_Ascale),Cscale(_Cscale){}
  
  ValueType value(const GeneralizedCoordinate &c, const ParameterType &p) const{
    if(c.idx == -1){
      return p.Apipi * Ascale * ( exp(-p.Epipi*c.t) + exp(-p.Epipi*(Lt-2.*tsep_pipi-c.t)) ) + p.Cpipi * Cscale;
    }else{    
      return p.AK*sqrt(p.Apipi*Ascale)*p.M[c.idx]*exp(-p.Epipi * c.tsep_k_pi)*exp( -(p.mK - p.Epipi)*c.t )/sqrt(2.);
    }
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &c, const ParameterType &p) const{
    ValueDerivativeType yderivs;  yderivs.zero();
    if(c.idx == -1){
      yderivs.Apipi = Ascale * ( exp(-p.Epipi*c.t) + exp(-p.Epipi*(Lt-2.*tsep_pipi-c.t)) );
      yderivs.Epipi = p.Apipi * Ascale * ( (-c.t)*exp(-p.Epipi*c.t) - (Lt-2.*tsep_pipi-c.t)*exp(-p.Epipi*(Lt-2.*tsep_pipi-c.t)) );
      yderivs.Cpipi = Cscale;    
    }else{
      double v = value(c,p);
      yderivs.AK = v/p.AK;
      yderivs.Apipi = 0.5*v/p.Apipi;
      yderivs.Epipi = -c.tsep_k_pi*v +c.t*v;
      yderivs.mK = -c.t*v;
      yderivs.M[c.idx] = v/p.M[c.idx];
    }
    return yderivs;
  }
    
  static inline int Nparams(){ return 12; }
};

















struct TwoExpTwoPointThreePointSimFitParams{
  double AK;
  double mK;
  double Apipi0;
  double Epipi0;
  double Apipi1;
  double Epipi1;  
  double Cpipi; //constant term
  std::vector<double> M;

  enum {nprelim = 7};
  
  TwoExpTwoPointThreePointSimFitParams(): AK(1), mK(0.5), Apipi0(0.16), Epipi0(0.36), Apipi1(0.1), Epipi1(1.0), Cpipi(0), M(7){
    M[0] = -0.003; M[1] = -0.003; M[2] = 0.005; M[3] = -0.01; M[4] = -0.02; M[5] = 0.02; M[6] = 0.07;
  }
  TwoExpTwoPointThreePointSimFitParams(const double _AK, const double _mK,
				       const double _Apipi0, const double _Epipi0,
				       const double _Apipi1, const double _Epipi1,
				       const double _Cpipi, const double _M): AK(_AK), mK(_mK),
									      Apipi0(_Apipi0), Epipi0(_Epipi0),
									      Apipi1(_Apipi1), Epipi1(_Epipi1),
									      Cpipi(_Cpipi), M(7,_M){}
  TwoExpTwoPointThreePointSimFitParams(const TwoExpTwoPointThreePointSimFitParams &r) = default;
  TwoExpTwoPointThreePointSimFitParams(TwoExpTwoPointThreePointSimFitParams &&r) = default;

  TwoExpTwoPointThreePointSimFitParams& operator=(const TwoExpTwoPointThreePointSimFitParams &r) = default;
  TwoExpTwoPointThreePointSimFitParams& operator=(TwoExpTwoPointThreePointSimFitParams &&r) = default;
    
  inline int size() const{ return 14; }

  inline std::string print() const{
    std::ostringstream os; os << "AK=" << AK << " mK=" << mK << " Apipi0=" << Apipi0 << " Epipi0=" << Epipi0 << " Apipi1=" << Apipi1 << " Epipi1=" << Epipi1 << " Cpipi=" << Cpipi;
    for(int i=0;i<7;i++) os << " M[" << i << "]=" << M[i];
    return os.str();
  }
  inline void zero(){
    AK=mK=Apipi0=Epipi0=Apipi1=Epipi1=Cpipi=0.;
    for(int i=0;i<7;i++) M[i] = 0.;
  }
  inline double &operator()(const int i){
    if(i<7){
      switch(i){
      case 0:
	return AK;
      case 1:
	return mK;
      case 2:
	return Apipi0;
      case 3:
	return Epipi0;
      case 4:
	return Apipi1;
      case 5:
	return Epipi1;	
      case 6:
	return Cpipi;
      }
    }else return M[i-7];
  }
  inline const double &operator()(const int i) const{ return const_cast<const double &>( const_cast<TwoExpTwoPointThreePointSimFitParams*>(this)->operator()(i) ); }

  typedef TwoExpTwoPointThreePointSimFitParams ET_tag;
  
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,TwoExpTwoPointThreePointSimFitParams >::value, int>::type = 0>
  TwoExpTwoPointThreePointSimFitParams(U&& expr): M(7){
    this->AK = expr[0];
    this->mK = expr[1];
    this->Apipi0 = expr[2];
    this->Epipi0 = expr[3];
    this->Apipi1 = expr[4];
    this->Epipi1 = expr[5];    
    this->Cpipi = expr[6];
    for(int i=7;i<14;i++) this->M[i-7] = expr[i];
  }
};
template<>
struct CPSfit::getElem<TwoExpTwoPointThreePointSimFitParams>{
  static inline auto elem(const TwoExpTwoPointThreePointSimFitParams &v, const int i)->decltype(v(i)){ return v(i); }
  static inline int common_properties(const TwoExpTwoPointThreePointSimFitParams &v){ return 0; }
};

template<>
struct CPSfit::printStats<jackknifeDistribution<TwoExpTwoPointThreePointSimFitParams> >{
  inline static std::string centralValue(const jackknifeDistribution<TwoExpTwoPointThreePointSimFitParams> &d){
    auto best = d.best();
    std::ostringstream os; os << "(" << best.AK << ", " << best.mK << ", " << best.Apipi0 << ", " << best.Epipi0 << ", " << best.Apipi1 << ", " << best.Epipi1 << ", " << best.Cpipi;
    for(int i=0;i<7;i++) os << ", " << best.M[i];
    os << ")";
    return os.str();
  }
  inline static std::string error(const jackknifeDistribution<TwoExpTwoPointThreePointSimFitParams> &d){
    auto err = d.standardError();
    std::ostringstream os; os << "(" << err.AK << ", " << err.mK << ", " << err.Apipi0 << ", " << err.Epipi0 << ", " << err.Apipi1 << ", " << err.Epipi1 << ", " << err.Cpipi;
    for(int i=0;i<7;i++) os << ", " << err.M[i];
    os << ")";
    return os.str();
  }
};

GENERATE_PARSER(TwoExpTwoPointThreePointSimFitParams, (double, AK)(double, mK)(double, Apipi0)(double, Epipi0)(double, Apipi1)(double, Epipi1)(double, Cpipi)(std::vector<double>, M) )


class FitTwoExpTwoPointThreePointSim{
  const double Lt;
  const double tsep_pipi;
  const double Ascale;
  const double Cscale;

public:
  typedef TwoExpTwoPointThreePointSimFitParams Params;
  typedef double ValueType;
  typedef Params ParameterType;
  typedef Params ValueDerivativeType; //derivative wrt parameters
  typedef amplitudeDataCoordSim GeneralizedCoordinate; //use idx = -1 for pipi data

  FitTwoExpTwoPointThreePointSim(const int _Lt, const double _tsep_pipi, const double _Ascale=1e13, const double _Cscale=1e13): Lt(_Lt),tsep_pipi(_tsep_pipi),Ascale(_Ascale),Cscale(_Cscale){}
  
  ValueType value(const GeneralizedCoordinate &c, const ParameterType &p) const{
    if(c.idx == -1){
      return p.Apipi0 * Ascale * ( exp(-p.Epipi0*c.t) + exp(-p.Epipi0*(Lt-2.*tsep_pipi-c.t)) )
	+ p.Apipi1 * Ascale * ( exp(-p.Epipi1*c.t) + exp(-p.Epipi1*(Lt-2.*tsep_pipi-c.t)) )
	+ p.Cpipi * Cscale;
    }else{    
      return p.AK*sqrt(p.Apipi0*Ascale)*p.M[c.idx]*exp(-p.Epipi0 * c.tsep_k_pi)*exp( -(p.mK - p.Epipi0)*c.t )/sqrt(2.);
    }
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &c, const ParameterType &p) const{
    ValueDerivativeType yderivs;  yderivs.zero();
    if(c.idx == -1){
      yderivs.Apipi0 = Ascale * ( exp(-p.Epipi0*c.t) + exp(-p.Epipi0*(Lt-2.*tsep_pipi-c.t)) );
      yderivs.Epipi0 = p.Apipi0 * Ascale * ( (-c.t)*exp(-p.Epipi0*c.t) - (Lt-2.*tsep_pipi-c.t)*exp(-p.Epipi0*(Lt-2.*tsep_pipi-c.t)) );
      
      yderivs.Apipi1 = Ascale * ( exp(-p.Epipi1*c.t) + exp(-p.Epipi1*(Lt-2.*tsep_pipi-c.t)) );
      yderivs.Epipi1 = p.Apipi1 * Ascale * ( (-c.t)*exp(-p.Epipi1*c.t) - (Lt-2.*tsep_pipi-c.t)*exp(-p.Epipi1*(Lt-2.*tsep_pipi-c.t)) );

      yderivs.Cpipi = Cscale;    
    }else{
      double v = value(c,p);
      yderivs.AK = v/p.AK;
      yderivs.Apipi0 = 0.5*v/p.Apipi0;
      yderivs.Epipi0 = -c.tsep_k_pi*v +c.t*v;
      yderivs.Apipi1 = 0;
      yderivs.Epipi1 = 0;
      yderivs.mK = -c.t*v;
      yderivs.M[c.idx] = v/p.M[c.idx];
    }
    return yderivs;
  }
    
  static inline int Nparams(){ return 14; }
};





#endif
