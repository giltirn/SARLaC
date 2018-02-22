#ifndef _CPSFIT_MPL_INTERFACE_2D_H_
#define _CPSFIT_MPL_INTERFACE_2D_H_

//Interface directly with Matplotlib in real-time using Boost's python interface
#include<config.h>

#ifdef HAVE_PYTHON

#include<utility>
#include<vector>
#include<boost/python.hpp>

#include<utils/macros.h>
#include<plot/plot/datasets_2d.h>

CPSFIT_START_NAMESPACE

class MatPlotLibInterface{
public:
  typedef boost::python::dict kwargsType;
  typedef std::pair<boost::python::object,SetType> handleType;
private:  
  boost::python::object Figure;
  boost::python::object FigureCanvas;

  boost::python::object patches;
  boost::python::object plt;
  boost::python::object figure;
  boost::python::object axes;
  boost::python::object canvas;

  std::vector<handleType> leg_handles;
  std::vector<std::string> legends;

  boost::python::object toHex(const boost::python::object &r, const boost::python::object &g, const boost::python::object &b){
    double rcd = boost::python::extract<double>(r);
    double gcd = boost::python::extract<double>(g);
    double bcd = boost::python::extract<double>(b);
    unsigned int rc(rcd);
    unsigned int gc(gcd);
    unsigned int bc(bcd);
    
    //std::ostringstream os; os << boost::format("#%02x%02x%02x") % rc % gc % bc;
    char buf[100];
    sprintf(buf,"#%02x%02x%02x",rc,gc,bc);    
    return boost::python::object(buf);
  }
  
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
      patches = import("matplotlib.patches");
	
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
  handleType plotData(const Data &data){
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
  handleType plotData(const Data &data, boost::python::dict &kwargs){
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

      return handleType(plotset, DataSetType);
      
    } catch( error_already_set ) {
      PyErr_Print();
    }
  }


  template<typename Band>
  handleType errorBand(const Band &band){
    boost::python::dict kwargs;
    return errorBand(band,kwargs);
  }
    
  //Band is an accessor with methods:
  //double x(const int i)
  //double upper(const int i)
  //double lower(const int i)
  //int size()
  template<typename Band>
  handleType errorBand(const Band &band, boost::python::dict &kwargs){
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

      if(boundary_lines){
        object color = plt.attr("getp")(plotband,"facecolors");
        object chex = toHex(color[0][0]*255,color[0][1]*255,color[0][2]*255);
	boost::python::dict args;
	args["marker"] = "None";
	args["linestyle"] = "--";
	args["linewidth"] = 1.5;
	args["color"] = chex;
	args["zorder"] = boundary_lines_zorder;
	
	axes.attr("plot")(*make_tuple(x,upper), **args);
	axes.attr("plot")(*make_tuple(x,lower), **args);
      }
	  
      return handleType(plotband, ErrorBandType);
    } catch( error_already_set ) {
      PyErr_Print();
    }
  }

  void setLegend(const handleType &handle, const std::string &to){
    leg_handles.push_back(handle);
    legends.push_back(to);    
  }
  void createLegend(){
    boost::python::dict kwargs;
    createLegend(kwargs);
  }
  void createLegend(boost::python::dict &kwargs){
    assert(leg_handles.size() == legends.size());
    using namespace boost::python;
    try{    
      list phandles;
      list plegends;
      for(int i=0;i<leg_handles.size();i++){
	if(leg_handles[i].second == ErrorBandType){
#if 0
	  object color = plt.attr("getp")(*make_tuple(leg_handles[i].first, "facecolors"));
	  color = toHex(255*color[0][0],255*color[0][1],255*color[0][2]);
	  boost::python::dict rargs; rargs["color"] = color;

	  object zz = *make_tuple(0,0);	  
	  object p = patches.attr("Rectangle")(*make_tuple(zz, 5.0, 10.0) , **rargs);
	  phandles.append(p);
#else
	  phandles.append(leg_handles[i].first);
#endif
	}else{
	  phandles.append(leg_handles[i].first[0]);
	}
	plegends.append(object(legends[i].c_str()));      
      }

      if(!kwargs.has_key("numpoints"))
	kwargs["numpoints"] = 1;
      if(!kwargs.has_key("shadow"))
	kwargs["shadow"] = true;
      if(!kwargs.has_key("fancybox"))
	kwargs["fancybox"] = true;
    
      axes.attr("legend")(*make_tuple(phandles,plegends),**kwargs);
    } catch( error_already_set ) {
      PyErr_Print();
    }      
  }	   

		  
  
};

CPSFIT_END_NAMESPACE
#endif
#endif
