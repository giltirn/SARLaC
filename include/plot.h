#ifndef _CPSFIT_PLOT_H
#define _CPSFIT_PLOT_H

#include <boost/python.hpp>
#include <map>
#include <sstream>

class MatPlotLibInterface{
public:
  typedef boost::python::dict kwargsType;

private:  
  boost::python::object Figure;
  boost::python::object FigureCanvas;

  boost::python::object plt;
  boost::python::object figure;
  boost::python::object axes;
  boost::python::object canvas;
public:
  MatPlotLibInterface(){
    using namespace boost::python;
    try {
    
      if(!Py_IsInitialized())
	Py_Initialize();
    
      PyRun_SimpleString("import matplotlib");

      Figure = object(handle<>(PyImport_ImportModule("matplotlib.figure"))).attr("Figure");
      FigureCanvas = object(handle<>(PyImport_ImportModule("matplotlib.backends.backend_pdf"))).attr("FigureCanvas");

      plt = import("matplotlib.pyplot");
      
      figure = Figure();
      axes = figure.attr("add_subplot")(111);
      canvas = FigureCanvas(figure);
    } catch( error_already_set ) {
      PyErr_Print();
    }
  }

  void write(const std::string &filename){
    using namespace boost::python;
    try {
      canvas.attr("draw")();
      figure.attr("savefig")(filename);
    } catch( error_already_set ) {
      PyErr_Print();
    }
  }

  template<typename Data>
  void plotData(const Data &data){
    boost::python::dict kwargs;
    return plotData(data,kwargs);
  }


  //Data should be an accessor wrapper that has methods
  //double x(const int),
  //double y(const int),
  //double dxm(const int)  [minus-error],
  //double dxp(const int) [plus-error],
  //double dym(const int),
  //double dyp(const int)
  template<typename Data>
  void plotData(const Data &data, boost::python::dict &kwargs ){
    using namespace boost::python;
    try {
      //Default arguments if not set
      if(!kwargs.has_key("linestyle"))
	kwargs["linestyle"] = "";
      if(!kwargs.has_key("marker"))
	kwargs["marker"] = "o";   //cf. http://matplotlib.org/api/markers_api.html
      if(!kwargs.has_key("ms"))
	kwargs["ms"] = 5;  //marker size
      
      //My extra arguments
      object capwidth(2.0);
      object bar_linestyle("solid");
      object color("r");
      bool hollowsymbol = false;
      
      if(kwargs.has_key("capwidth")){
	capwidth = kwargs["capwidth"];
	kwargs["capwidth"].del();
      }
      if(kwargs.has_key("bar_linestyle")){ //ACCEPTS: ['solid' | 'dashed', 'dashdot', 'dotted' | (offset, on-off-dash-seq) ]
	capwidth = kwargs["bar_linestyle"];
	kwargs["bar_linestyle"].del();
      }     
      if(kwargs.has_key("hollowsymbol") && kwargs["hollowsymbol"] == 1){
	hollowsymbol = true;
	kwargs["hollowsymbol"].del();
      }
      if(kwargs.has_key("color")){ //if it has color specified, replace default
	color = kwargs["color"];
	kwargs["color"].del();
      }
      if(!kwargs.has_key("ecolor")) //only replace if not specified
	kwargs["ecolor"] = color;
      if(!kwargs.has_key("mfc"))
	kwargs["mfc"] = color;
      
      
      int sz = data.size();
      list px, py, pdxm, pdxp, pdym, pdyp;
      for(int i=0;i<sz;i++){
      	px.append(data.x(i));
      	py.append(data.y(i));

	pdxm.append(data.dxm(i));
	pdxp.append(data.dxp(i));
	
	pdym.append(data.dym(i));
	pdyp.append(data.dyp(i));
      }

      list pdxpm;
      pdxpm.append(pdxm);
      pdxpm.append(pdxp);

      list pdypm;
      pdypm.append(pdym);
      pdypm.append(pdyp);
      
      list args;
      args.append(px);
      args.append(py);

      kwargs["xerr"]=pdxpm;
      kwargs["yerr"]=pdypm;
      
      object plotset = axes.attr("errorbar")(*tuple(args), **kwargs);

      //Linestyle
      plotset[2][0].attr("set_linestyles")(bar_linestyle);
      plotset[2][1].attr("set_linestyles")(bar_linestyle);

      //Cap size
      {
	int ncap = len(plotset[1]);
	dict setparg; setparg["mew"] = capwidth;
	for(int i=0;i<ncap;i++)
	  plt.attr("setp")(*make_tuple(plotset[1][i]),**setparg);
      }

      //Hollow symbols if required
      if(hollowsymbol){
	object s = plotset[0];
	s.attr("set_markerfacecolor")("None");
	s.attr("set_markeredgecolor")(color);
	if(!kwargs.has_key("markeredgewidth"))
	  s.attr("set_markeredgewidth")(1.25);
      }
      
    } catch( error_already_set ) {
      PyErr_Print();
    }
  }


  template<typename Band>
  void errorBand(const Band &band){
    boost::python::dict kwargs;
    return errorBand(band,kwargs);
  }
    
  //Band is an accessor with methods:
  //double x(const int i)
  //double upper(const int i)
  //double lower(const int i)
  //int size()
  template<typename Band>
  void errorBand(const Band &band, boost::python::dict &kwargs){
    using namespace boost::python;
    try{
      int sz = band.size();

      list x,upper, lower;
      for(int i=0;i<sz;i++){
	x.append(band.x(i));
	upper.append(band.upper(i));
	lower.append(band.lower(i));
      }
    
      object usetransparency(0);
      //object alpha(0.2);
      bool boundary_lines = false;
      object boundary_lines_zorder(10);

      if(kwargs.has_key("usetransparency")){
	usetransparency=kwargs["usetransparency"];
	kwargs["usetransparency"].del();
      }
      // if(kwargs.has_key("alpha")){
      //   alpha = kwargs["alpha"];
      // }
      if(kwargs.has_key("boundary_lines") && kwargs["boundary_lines"] == 1){  //draw dashed lines definining the boundary of the errorband
	boundary_lines = true;
	kwargs["boundary_lines"].del();
      }
      if(kwargs.has_key("boundary_lines_zorder")){ //change the zorder of the boundary lines. Default 10
	boundary_lines_zorder = kwargs["boundary_lines_zorder"];
	kwargs["boundary_lines_zorder"].del();
      }

      object plotband = axes.attr("fill_between")(*make_tuple(x,lower,upper), **kwargs);

      // if(boundary_lines){
      //   object color = plt.attr("getp")(plotband,"facecolors");
      //   chex = ColourPallete.toHex(colour[0][0]*255,colour[0][1]*255,colour[0][2]*255)
      //         self.ax.plot(x,uy,marker='None',linestyle='--',linewidth=1.5,color=chex,zorder=boundary_lines_zorder)
      //         self.ax.plot(x,ly,marker='None',linestyle='--',linewidth=1.5,color=chex,zorder=boundary_lines_zorder)
    } catch( error_already_set ) {
      PyErr_Print();
    }
  }
};


class OstreamHook{
public:
  virtual void write(std::ostream &) const = 0;
};

inline std::ostream & operator<<(std::ostream &os, const OstreamHook &hk){
  hk.write(os);
  return os;
}




class KWargElem: public OstreamHook{
  std::string val;

public:
  KWargElem(){}
    
  template<typename T>
  explicit KWargElem(const T& v){
    *this = v;
  }

  template<typename T>
  KWargElem &operator=(const T &v){
    std::ostringstream os; os << v;
    val = os.str();
    return *this;
  }

  KWargElem &operator=(const char* v){
    val = "\"" + std::string(v) + "\"";
    return *this;
  }
  
  KWargElem &operator=(const std::string &v){
    val = "\"" + v + "\"";
    return *this;
  }
  KWargElem &operator=(const char v){
    std::ostringstream os; os << "\'" << v << "\'";
    val = os.str();
    return *this;
  }
  
  inline void write(std::ostream &os) const{
    os << val;
  }
};


struct ListPrint: public OstreamHook{
  const std::vector<double> &lst;
  ListPrint(const std::vector<double> &_lst): lst(_lst){}
  void write(std::ostream &os) const{
    os << "[";
    for(int i=0;i<lst.size()-1; i++) os << lst[i] << ',';      
    os << lst.back() << "]";
  }
};

class PythonDataContainer{
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> dxm;
  std::vector<double> dxp;
  std::vector<double> dym;
  std::vector<double> dyp;
  std::string set_tag;

public:

  template<typename Data>
  void import(const Data &data, const std::string &tag){
    set_tag = tag;
    
    int sz = data.size();
    for(int i=0;i<sz;i++){
      x.push_back(data.x(i));
      y.push_back(data.y(i));
      
      dxm.push_back(data.dxm(i));
      dxp.push_back(data.dxp(i));
      
      dym.push_back(data.dym(i));
      dyp.push_back(data.dyp(i));
    }
  }

  void write(std::ostream &os) const{
    os << set_tag << "= pyplot.DataSet()\n";
    os << set_tag << ".x = " << ListPrint(x) << '\n';
    os << set_tag << ".y = " << ListPrint(y) << '\n';
    os << set_tag << ".dxm = " << ListPrint(dxm) << '\n';
    os << set_tag << ".dxp = " << ListPrint(dxp) << '\n';
    os << set_tag << ".dym = " << ListPrint(dym) << '\n';
    os << set_tag << ".dyp = " << ListPrint(dyp) << "\n\n";
  }

  const std::string &tag() const{ return set_tag; }
  
};


class PythonErrorBandContainer{
  std::vector<double> x;
  std::vector<double> upper;
  std::vector<double> lower;
  std::string set_tag;

public:
  
  template<typename Band>
  void import(const Band &band, const std::string &tag){
    set_tag = tag;
    
    int sz = band.size();
    for(int i=0;i<sz;i++){
      x.push_back(band.x(i));
      upper.push_back(band.upper(i));
      lower.push_back(band.lower(i));
    }
  }

  void write(std::ostream &os) const{
    os << set_tag << "= pyplot.ErrorBand()\n";
    os << set_tag << ".x = " << ListPrint(x) << '\n';
    os << set_tag << ".upper = " << ListPrint(upper) << '\n';
    os << set_tag << ".lower = " << ListPrint(lower) << '\n';
  }

  const std::string &tag() const{ return set_tag; }  
};




class MatPlotLibScriptGenerate{
public:
  typedef std::map<std::string,KWargElem> kwargsType;
private:
  
  std::vector<PythonDataContainer> plotdata_sets;
  std::vector<kwargsType> plotdata_args;

  std::vector<PythonErrorBandContainer> ploterrorband_sets;
  std::vector<kwargsType> ploterrorband_args;

  struct kwargsPrint: public OstreamHook{
    const kwargsType &args;
    kwargsPrint(const kwargsType &_args): args(_args){}
    void write(std::ostream &os) const{
      for(kwargsType::const_iterator it = args.begin(); it != args.end(); it++){
	os << ',' << it->first << "=" << it->second;
      }
    }
  };
  
public:

  //Data should be an accessor wrapper that has methods
  //double x(const int),
  //double y(const int),
  //double dxm(const int)  [minus-error],
  //double dxp(const int) [plus-error],
  //double dym(const int),
  //double dyp(const int)
  //int size()
  template<typename Data>
  void plotData(const Data &data, kwargsType &kwargs = kwargsType()){
    std::ostringstream os; os << "dset" << plotdata_sets.size();
    plotdata_sets.push_back(PythonDataContainer());
    plotdata_sets.back().import(data, os.str());
    plotdata_args.push_back(kwargs);
  }

  //Band is an accessor with methods:
  //double x(const int i)
  //double upper(const int i)
  //double lower(const int i)
  //int size()
  template<typename Band>
  void errorBand(const Band &band, kwargsType &kwargs = kwargsType()){
    std::ostringstream os; os << "band" << ploterrorband_sets.size();
    ploterrorband_sets.push_back(PythonErrorBandContainer());
    ploterrorband_sets.back().import(band, os.str());
    ploterrorband_args.push_back(kwargs);
  }
  
  void write(std::ostream &os, const std::string &script_gen_filename = "plot.pdf") const{
    os << "import pyplot\n\n";

    for(int i=0;i<plotdata_sets.size();i++)
      plotdata_sets[i].write(os);

    for(int i=0;i<ploterrorband_sets.size();i++)
      ploterrorband_sets[i].write(os);

    os << "\nfig = pyplot.plt.figure()\n";
    os << "ax = fig.add_subplot(1,1,1)\n";

    for(int i=0;i<plotdata_sets.size();i++)
      os << "plot_" << plotdata_sets[i].tag() << " = " << "pyplot.plotDataSet(ax, " << plotdata_sets[i].tag() << kwargsPrint(plotdata_args[i]) << ")\n";

    for(int i=0;i<ploterrorband_sets.size();i++)
      os << "plot_" << ploterrorband_sets[i].tag() << " = " << "pyplot.plotErrorBand(ax, " << ploterrorband_sets[i].tag() << kwargsPrint(ploterrorband_args[i]) << ")\n";    
    
    os << "fig.canvas.draw()\n";
    os << "fig.savefig(\"" << script_gen_filename << "\")";
  }
  void write(const std::string &filename, const std::string &script_gen_filename = "plot.pdf"){
    std::ofstream of(filename.c_str());
    write(of,script_gen_filename);
    of.close();
  }
};




//Example accessor for x,y data in std::vectors with symmetric errors
//Can also be used to define error bands
class DataVectorAccessor{
  const std::vector<double> &_x;
  const std::vector<double> &_y;
  const std::vector<double> &_dx;
  const std::vector<double> &_dy;
  int sz;
public:
  DataVectorAccessor(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &dx, const std::vector<double> &dy): _x(x), _y(y), _dx(dx), _dy(dy), sz(y.size()){
    assert(x.size() == sz && dx.size() == sz && dy.size() == sz);
  }

  inline double x(const int i) const{ return _x[i]; }
  inline double y(const int i) const{ return _y[i]; }
  inline double dxm(const int i) const{ return _dx[i]; }
  inline double dxp(const int i) const{ return _dx[i]; }
  inline double dym(const int i) const{ return _dy[i]; }
  inline double dyp(const int i) const{ return _dy[i]; }

  inline double upper(const int i) const{ return _y[i]+_dy[i]; }
  inline double lower(const int i) const{ return _y[i]-_dy[i]; }
  
  inline int size() const{ return sz; }
};

class BandVectorAccessor{
  const std::vector<double> &_x;
  const std::vector<double> &_upper;
  const std::vector<double> &_lower;
  int sz;
public:
  BandVectorAccessor(const std::vector<double> &x, const std::vector<double> &upper, const std::vector<double> &lower): _x(x), _upper(upper), _lower(lower), sz(upper.size()){
    assert(lower.size() == sz && x.size() == sz);
  }
  inline double x(const int i) const{ return _x[i]; }
  inline double upper(const int i) const{ return _upper[i]; }
  inline double lower(const int i) const{ return _lower[i]; }
  inline int size() const{ return sz; }
};


//CoordinatePolicy is converts the underlying coordinate type into double central values and errors
//ValuePolicy does the same for the y data
template<typename DataSeriesType, typename CoordinatePolicy, typename ValuePolicy>
class DataSeriesAccessor: public ValuePolicy, public CoordinatePolicy{
  const DataSeriesType &series;
public:
  DataSeriesAccessor(const DataSeriesType &_series):series(_series){}
  
  int size() const{ return series.size(); }

  inline double x(const int i) const{ return this->CoordinatePolicy::value(series.coord(i)); }
  inline double dxp(const int i) const{ return this->CoordinatePolicy::errplus(series.coord(i)); }
  inline double dxm(const int i) const{ return this->CoordinatePolicy::errminus(series.coord(i)); }  

  inline double y(const int i) const{ return this->ValuePolicy::value(series.value(i)); }
  inline double dyp(const int i) const{ return this->ValuePolicy::errplus(series.value(i)); }
  inline double dym(const int i) const{ return this->ValuePolicy::errminus(series.value(i)); }  

  inline double upper(const int i) const{ return y(i)+dyp(i); }
  inline double lower(const int i) const{ return y(i)-dym(i); }
};

template<typename DistributionType>
class DistributionPlotAccessor{
public:
  static inline double value(const DistributionType &d){ return d.mean(); }
  static inline double errplus(const DistributionType &d){ return d.standardError(); }
  static inline double errminus(const DistributionType &d){ return d.standardError(); }  
};
//Coordinate accessor that relies on implicit conversion of type to double
template<typename T>
class ScalarCoordinateAccessor{
public:
  static inline double value(const T &d){ return d; }
  static inline double errplus(const T &d){ return 0; } //zero error
  static inline double errminus(const T &d){ return 0; }  
};

#endif
