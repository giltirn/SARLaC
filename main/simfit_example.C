#include <fstream>
#include <minimizer.h>
#include <cost_function.h>
#include <fitfunc.h>
#include <fitfunc_mapping.h>
#include <random.h>
#include <plot.h>
#include <distribution.h>
#include <data_series.h>
#include <parser.h>
#include <common_defs.h>
#include <sstream>

#include <generic_ET.h>

//Coordinate, parameters and param derivatives for aggregate fit func
struct Coord{
  double t; //time
  int idx; //data type index
};
std::ostream & operator<<(std::ostream &os, const Coord &c){
  os << "(" << c.t << "," << c.idx << ")";
  return os;
}


template<typename DataSeriesType>
class SingleTypePlotAccessor{
  const DataSeriesType &series;
  std::vector<int> sidx;
public:
  SingleTypePlotAccessor(const DataSeriesType &_series, const int _didx):series(_series){
    for(int i=0;i<series.size();i++) if(series.coord(i).idx == _didx) sidx.push_back(i);
  }
  
  int size() const{ return sidx.size(); }

  inline double x(const int i) const{ return series.coord(sidx[i]).t ; }
  inline double dxp(const int i) const{ return 0; }
  inline double dxm(const int i) const{ return 0; }  

  inline double y(const int i) const{ return series.value(sidx[i]).mean(); }
  inline double dyp(const int i) const{ return series.value(sidx[i]).standardError(); }
  inline double dym(const int i) const{ return series.value(sidx[i]).standardError(); }  

  inline double upper(const int i) const{ return y(i)+dyp(i); }
  inline double lower(const int i) const{ return y(i)-dym(i); }
};



class Params{
public:
  double a;
  double b;
  double c;
  double d;

  inline double& operator()(const int i){
    switch(i){
    case 0:
      return a;
    case 1:
      return b;
    case 2:
      return c;
    case 3:
      return d;
    default:
      assert(0);
    }
  }
  inline const double& operator()(const int i) const{ return const_cast<const double&>( const_cast<Params*>(this)->operator()(i) ); }

  inline void zero(){ a=b=c=d=0.; }
  
  std::string print() const{
    std::ostringstream os; os << "(a="<<a<<", b="<<b<<", c="<<c<<", d="<<d <<")";
    return os.str();
  }
  inline int size() const{ return 4; }

  Params &operator+=(const Params &p){
    a += p.a;
    b += p.b;
    c += p.c;
    d += p.d;
    return *this;
  }
  Params &operator=(const double v){
    a = b = c =d = v;
  }

  Params(){}

  typedef Params ET_tag;
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,Params>::value, int>::type = 0>
  Params(U&& expr){
    for(int i=0;i<4;i++) (*this)(i) = expr[i];
  }
  
};

template<>
struct getElem<Params>{
  static inline auto elem(const Params &v, const int i)->decltype(v(i)){ return v(i); }
  static inline auto elem(Params &v, const int i)->decltype(v(i)){ return v(i); }
  static int common_properties(const Params &v){ return 0; }
};







std::ostream & operator<<(std::ostream &os, const Params &p){
  os << p.print();
  return os;
}



struct ParamDerivs{
  double df_by_da;
  double df_by_db;
  double df_by_dc;
  double df_by_dd;

  inline double& operator()(const int i){
    switch(i){
    case 0:
      return df_by_da;
    case 1:
      return df_by_db;
    case 2:
      return df_by_dc;
    case 3:
      return df_by_dd;
    default:
      assert(0);
    }
  }
  inline const double& operator()(const int i) const{ return const_cast<const double&>( const_cast<ParamDerivs*>(this)->operator()(i) ); }

  inline void zero(){ df_by_da=df_by_db=df_by_dc=df_by_dd=0.; }
  
  std::string print() const{
    std::ostringstream os; os << "(d/da="<<df_by_da<<", d/db="<<df_by_db<<", d/dc="<<df_by_dc<<", d/dd="<<df_by_dd <<")";
    return os.str();
  }
  inline int size() const{ return 4; }  
};

//Accessors (views) to pass values to/from the underlying fit functions
struct CoordTview{
  const Coord &c;
  CoordTview(const Coord &_c): c(_c){}
  
  inline double operator[](const int i) const{ //interface that NumericLinearFit wants
    assert(i==0);
    return c.t;
  }
};
struct ParamView{
  Params &pfull;
  int type_idx;
  
  ParamView(Params &_pfull, int _type_idx): pfull(_pfull), type_idx(_type_idx){}
  
  inline double & operator()(const int i){ return pfull(2*type_idx + i); }
  inline const double & operator()(const int i) const{ return const_cast<const double &>( const_cast<ParamView*>(this)->operator()(i)); }
};
struct ParamDerivsView{
  ParamDerivs &pfull;
  int type_idx;
  
  ParamDerivsView(ParamDerivs &_pfull, int _type_idx): pfull(_pfull), type_idx(_type_idx){}
  
  inline double & operator()(const int i){ return pfull(2*type_idx + i); }
  inline const double & operator()(const int i) const{ return const_cast<const double &>( const_cast<ParamDerivsView*>(this)->operator()(i)); }

  void resize(const int sz){ assert(sz == 2); }

  ParamDerivsView & operator=(const NumericVector<double> &sub){
    for(int i=0;i<sub.size();i++)
      (*this)(i) = sub(i);
    return *this;
  }
};




class TwoTypeLinearFit{
public:
  typedef double ValueType;
  typedef Params ParameterType;  //a b c d
  typedef ParamDerivs ValueDerivativeType; //derivative wrt parameters
  typedef Coord GeneralizedCoordinate;

private:
  NumericLinearFit<CoordTview,double,1, ParamView> fitfunc_a; //a + b*t
  NumericLinearFit<CoordTview,double,1, ParamView> fitfunc_b; //c + d*t
public:
  
  ValueType value(const GeneralizedCoordinate &coord, const ParameterType &params) const{
    CoordTview c(coord);
    ParamView p(const_cast<ParameterType &>(params),coord.idx);
    
    switch(coord.idx){
    case 0:
      return fitfunc_a.value(c,p);      
    case 1:
      return fitfunc_b.value(c,p);
    };
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &coord, const ParameterType &params) const{
    ValueDerivativeType yderivs;
    yderivs.zero();
    CoordTview c(coord);
    ParamView p(const_cast<ParameterType &>(params),coord.idx);
    ParamDerivsView pd(yderivs,coord.idx);
    
    switch(coord.idx){
    case 0:
      pd = fitfunc_a.parameterDerivatives(c,p); break;
    case 1:
      pd = fitfunc_b.parameterDerivatives(c,p); break;
    };
    return yderivs;
  }

  inline int Nparams() const{ return fitfunc_a.Nparams() + fitfunc_b.Nparams(); }
};


class filterCoordTrange{
  double t_min;
  double t_max;
  bool _invert;
public:
  filterCoordTrange(const double _t_min, const double _t_max): t_min(_t_min), t_max(_t_max), _invert(false){}

  void invert(){ _invert = !_invert; } 
  
  template<typename T>
  bool accept(const Coord &x, const T &y) const{
    bool cond = x.t >= t_min && x.t <= t_max;
    return _invert ? !cond : cond;
  }
};




int main(void){  
  RNG.initialize(1234);

  int nsample = 100;
  int npoints = 10;
  double t_min = 2;
  double t_max = 8;
  

  //Generate some random data
  double slope_a = 0.3;
  double slope_b = -0.7;
  double off_a = 1.4;
  double off_b = -0.9;
  double sigma_a = 0.7*sqrt(nsample);
  double sigma_b = 1.3*sqrt(nsample);
  
  publicationPrint<> printer;

  typedef dataSeries<Coord, distributionD> rawTimeSeriesType;

  rawTimeSeriesType data(2*npoints, nsample);
  
  for(int i=0;i<npoints;i++){
    data.coord(i).t = i;
    data.coord(i+npoints).t = i;

    data.coord(i).idx = 0;
    data.coord(i+npoints).idx = 1;

    for(int j=0;j<nsample;j++){
      data.value(i).sample(j) = gaussianRandom<double>(off_a + slope_a*i , sigma_a);
      data.value(i+npoints).sample(j) = gaussianRandom<double>(off_b + slope_b*i , sigma_b);
    }
  }
  
  //Setup fit range
  filterCoordTrange tfilter(t_min,t_max);

  typedef filteredDataSeries<rawTimeSeriesType> filteredRawTimeSeriesType;

  filteredRawTimeSeriesType inrange_data(data,tfilter);
  tfilter.invert();
  filteredRawTimeSeriesType outofrange_data(data,tfilter);
  
  //Resample data
  typedef dataSeries<Coord, jackknifeDistributionD> jackknifeTimeSeriesType;

  jackknifeTimeSeriesType inrange_jackknife;   resample(inrange_jackknife, inrange_data);
  jackknifeTimeSeriesType outofrange_jackknife;   resample(outofrange_jackknife, outofrange_data);

  const int ndata_in_range = inrange_jackknife.size();
  
  std::cout << "Resampled data in range:\n";
  
  for(int i=0;i<ndata_in_range;i++){
    printer << "t=" << inrange_jackknife.coord(i).t << " idx=" << inrange_jackknife.coord(i).idx  << " value=" << inrange_jackknife.value(i) << std::endl;    
  }

  //Setup fit function
  TwoTypeLinearFit func;

  //Setup cost function, acts on individual jackknife samples so need accessor
  std::vector<double> sigma(ndata_in_range);
  for(int i=0;i<sigma.size();i++) sigma[i] = inrange_jackknife.value(i).standardError();

  typedef sampleSeries<const jackknifeTimeSeriesType> sampleSeriesConstType; //const access
  typedef UncorrelatedChisqCostFunction<TwoTypeLinearFit, sampleSeriesConstType, double, NumericVector<double> > CostFunctionType;
  typedef CostFunctionType::CostType CostType;
  

  //Setup minimizer
  typedef MarquardtLevenbergMinimizer<CostFunctionType> MinimizerType;

  MarquardtLevenbergParameters<CostType> mlparams;
  mlparams.output = &null_stream;

  jackknifeDistribution<Params> params(nsample);

#pragma omp parallel for
  for(int j=0;j<nsample;j++){    
    sampleSeriesConstType dsample(inrange_jackknife, j);

    CostFunctionType costfunc(func, dsample, sigma);

    MinimizerType fitter(costfunc, mlparams);
    
    Params &pj = params.sample(j);
    //Guesses
    pj.a = 0.;
    pj.b = 1.0;
    pj.c = 0.;
    pj.d = 1.0; 

    CostType cost = fitter.fit(pj);
    assert(fitter.hasConverged());
  }

  std::cout << "Params: " << params.mean() << " " << params.standardError() << std::endl;

  //Plot result
#define USE_MPL_INTERFACE
#if defined(HAVE_PYTHON) && defined(USE_MPL_INTERFACE)
  typedef MatPlotLibInterface Plotter;
#else
  typedef MatPlotLibScriptGenerate Plotter;  
#endif
  
  Plotter mpl;
  typedef SingleTypePlotAccessor<jackknifeTimeSeriesType> PlotDataAccessor;

  typedef Plotter::handleType handleType;
  typedef Plotter::kwargsType kwargsType;

  {
    kwargsType args; args["color"] = 'r';    
    PlotDataAccessor datainterface(inrange_jackknife,0);
    handleType handle = mpl.plotData(datainterface,args);
    mpl.setLegend(handle,"set 0");
  }
  {
    kwargsType args; args["color"] = 'b';    
    PlotDataAccessor datainterface(inrange_jackknife,1);
    mpl.plotData(datainterface,args);    
  }

  
  {
    PlotDataAccessor datainterface(outofrange_jackknife,0);
    kwargsType args; args["hollowsymbol"] = 1; args["color"] = 'r';   
    mpl.plotData(datainterface,args);
  }
  {
    PlotDataAccessor datainterface(outofrange_jackknife,1);
    kwargsType args; args["hollowsymbol"] = 1; args["color"] = 'b';  
    mpl.plotData(datainterface,args);
  }

  //Generate fit curve
  int npt = 100;
  jackknifeTimeSeriesType fit_0(npt,nsample);
  jackknifeTimeSeriesType fit_1(npt,nsample);
  
  double delta = (t_max - t_min)/(npt-1);
#pragma omp parallel for
  for(int i=0;i<npt;i++){
    fit_0.coord(i).t = t_min + i*delta;
    fit_0.coord(i).idx = 0;

    fit_1.coord(i).t = t_min + i*delta;
    fit_1.coord(i).idx = 1;

    for(int j=0;j<nsample;j++){
      fit_0.value(i).sample(j) = func.value(fit_0.coord(i), params.sample(j));
      fit_1.value(i).sample(j) = func.value(fit_1.coord(i), params.sample(j));
    }
  }
  
  {
    PlotDataAccessor curveinterface(fit_0,0);
    kwargsType args; args["alpha"] = 0.3; args["color"] = "r"; args["boundary_lines"] = true;
    handleType handle = mpl.errorBand(curveinterface,args);
    mpl.setLegend(handle,"fit to set 0");
  }
  {
    PlotDataAccessor curveinterface(fit_1,1);
    kwargsType args; args["alpha"] = 0.3; args["color"] = "b"; args["boundary_lines"] = true;
    mpl.errorBand(curveinterface,args);
  }

  mpl.createLegend();
  
  mpl.write("plot.pdf");
  //mpl.write("test.py");


  
  std::cout << "Normal exit\n"; std::cout.flush();
  return 0;
}

