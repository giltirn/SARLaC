#ifndef _SARLAC_MPL_SCRIPT_GENERATE_2D_H_
#define _SARLAC_MPL_SCRIPT_GENERATE_2D_H_

//Generate python scripts for generating 2D plots using Matplotlib

#include<iostream>
#include<vector>

#include<config.h>
#include<utils/macros.h>
#include<plot/plot/write_python.h>
#include<plot/plot/datasets_2d.h>
#include<plot/plot/MPLscriptgenerate_base.h>

SARLAC_START_NAMESPACE

//To run the scripts generated by this class, make sure to setup your PYTHONPATH:
//export PYTHONPATH=${SRCDIR}/src:${PYTHONPATH}
//where SRCDIR is the root source directory of this library
class MatPlotLibScriptGenerate: public MatPlotLibScriptGenerateBase{
public:
  typedef MatPlotLibScriptGenerateBase::kwargsType kwargsType;
  typedef std::pair<int,SetType> handleType;
private:
  
  std::vector<PythonDataContainer> plotdata_sets;
  std::vector<kwargsType> plotdata_args;

  std::vector<PythonErrorBandContainer> ploterrorband_sets;
  std::vector<kwargsType> ploterrorband_args;

  std::vector<PythonHistogramContainer> plothistogram_sets;
  std::vector<kwargsType> plothistogram_args;
  
  std::vector<PythonErrorLineContainer> ploterrorlines_sets;
  std::vector<kwargsType> ploterrorlines_args;

  std::vector<PythonErrorCurveContainer> ploterrorcurves_sets;
  std::vector<kwargsType> ploterrorcurves_args;

  std::vector<handleType> leg_handles;
  std::vector<std::string> legends;
  std::string leg_py;
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
  handleType plotData(const Data &data, const kwargsType &kwargs, const std::string &tag = ""){
    plotdata_sets.push_back(PythonDataContainer());   
    if(tag == ""){
      std::ostringstream os; os << "dset" << plotdata_sets.size();
      plotdata_sets.back().import(data, os.str());
    }else{
      plotdata_sets.back().import(data, tag);
    }
    plotdata_args.push_back(kwargs);
    return handleType(plotdata_sets.size()-1, DataSetType);
  }
  template<typename Data>
    inline handleType plotData(const Data &data, const std::string &tag = ""){
    kwargsType kwargs;
    return plotData(data,kwargs,tag);
  }
    
  //Band is an accessor with methods:
  //double x(const int i)
  //double upper(const int i)
  //double lower(const int i)
  //int size()
  template<typename Band>
  handleType errorBand(const Band &band, const kwargsType &kwargs, const std::string &tag = ""){
    ploterrorband_sets.push_back(PythonErrorBandContainer());
    if(tag == ""){
      std::ostringstream os; os << "band" << ploterrorband_sets.size();
      ploterrorband_sets.back().import(band, os.str());
    }else{
      ploterrorband_sets.back().import(band, tag);
    }
    ploterrorband_args.push_back(kwargs);
    return handleType(ploterrorband_sets.size()-1,ErrorBandType);
  }
  template<typename Band>
  inline handleType errorBand(const Band &band, const std::string &tag = ""){
    kwargsType kwargs;
    return errorBand(band,kwargs,tag);
  }

  //Lines is an accessor with methods:
  //int size()
  //PythonTuple<double> start(const int i)
  //PythonTuple<double> end(const int i)
  template<typename Lines>
  handleType errorLines(const Lines &el, const kwargsType &kwargs, std::string tag = ""){
    ploterrorlines_sets.push_back(PythonErrorLineContainer());
    if(tag == "") tag = stringize("eline%d", ploterrorlines_sets.size());
    ploterrorlines_sets.back().import(el, tag);
    ploterrorlines_args.push_back(kwargs);
    return handleType(ploterrorlines_sets.size()-1,ErrorLineType);
  }
  template<typename Lines>
  inline handleType errorLines(const Lines &el, std::string tag = ""){
    kwargsType kwargs;
    return errorLines(el,kwargs,tag);
  }

  //Curves is an accessor with methods:
  //int ncurves()   - number of curves to plot
  //std::vector<pythonTuple<double> >  curves(const int i)  //coordinates of chosen points along curve i
  //int nmarkers()   - number of markers to plot
  //PythonTuple<double>  markers(const int i)  //coordinates of marker i
  template<typename Curves>
  handleType errorCurves(const Curves &el, const kwargsType &kwargs, std::string tag = ""){
    ploterrorcurves_sets.push_back(PythonErrorCurveContainer());
    if(tag == "") tag = stringize("ecurve%d", ploterrorcurves_sets.size());
    ploterrorcurves_sets.back().import(el, tag);
    ploterrorcurves_args.push_back(kwargs);
    return handleType(ploterrorcurves_sets.size()-1,ErrorCurveType);
  }
  template<typename Curves>
  inline handleType errorCurves(const Curves &el, std::string tag = ""){
    kwargsType kwargs;
    return errorCurves(el,kwargs,tag);
  }


  //Plot a histogram of the scalar data provided
  //Data is an accessors with methods:
  //double y(const int)
  //int size()
  template<typename Data>
  handleType histogram(const Data &data, const kwargsType &kwargs, const std::string &tag = ""){
    plothistogram_sets.push_back(PythonHistogramContainer());
    if(tag == ""){
      std::ostringstream os; os << "histogram" << plothistogram_sets.size();
      plothistogram_sets.back().import(data, os.str());
    }else{
      plothistogram_sets.back().import(data, tag);
    }
    plothistogram_args.push_back(kwargs);
    return handleType(plothistogram_sets.size()-1,HistogramType);
  }
  template<typename Data>
  inline handleType histogram(const Data &data, const std::string &tag = ""){
    kwargsType kwargs;
    return histogram(data,kwargs,tag);
  }

  
  void write(std::ostream &os, const std::string &script_gen_filename = "plot.pdf") const{
    os << "import pyplot\n";
    os << "import matplotlib\n\n";

    for(int i=0;i<plotdata_sets.size();i++)
      plotdata_sets[i].write(os);

    for(int i=0;i<ploterrorband_sets.size();i++)
      ploterrorband_sets[i].write(os);

    for(int i=0;i<plothistogram_sets.size();i++)
      plothistogram_sets[i].write(os);

    for(int i=0;i<ploterrorlines_sets.size();i++)
      ploterrorlines_sets[i].write(os);

    for(int i=0;i<ploterrorcurves_sets.size();i++)
      ploterrorcurves_sets[i].write(os);

    os << preamble.str(); //user code
    
    os << "\nif __name__ == '__main__':\n";

    os << "\tfig = pyplot.plt.figure()\n";
    os << "\tax = fig.add_subplot(1,1,1)\n";

    for(int i=0;i<plotdata_sets.size();i++)
      os << "\tplot_" << plotdata_sets[i].tag() << " = " << "pyplot.plotDataSet(ax, " << plotdata_sets[i].tag() << kwargsPrint(plotdata_args[i]) << ")\n";

    for(int i=0;i<ploterrorband_sets.size();i++)
      os << "\tplot_" << ploterrorband_sets[i].tag() << " = " << "pyplot.plotErrorBand(ax, " << ploterrorband_sets[i].tag() << kwargsPrint(ploterrorband_args[i]) << ")\n";    

    for(int i=0;i<plothistogram_sets.size();i++)
      os << "\tplot_" << plothistogram_sets[i].tag() << " = " << "pyplot.plotHistogram(ax, " << plothistogram_sets[i].tag() << kwargsPrint(plothistogram_args[i]) << ")\n";
    
    for(int i=0;i<ploterrorlines_sets.size();i++)
      os << "\tplot_" << ploterrorlines_sets[i].tag() << " = " << "pyplot.plotErrorLines(ax, " << ploterrorlines_sets[i].tag() << kwargsPrint(ploterrorlines_args[i]) << ")\n";

    for(int i=0;i<ploterrorcurves_sets.size();i++)
      os << "\tplot_" << ploterrorcurves_sets[i].tag() << " = " << "pyplot.plotErrorCurves(ax, " << ploterrorcurves_sets[i].tag() << kwargsPrint(ploterrorcurves_args[i]) << ")\n";

    os << user.str(); //user code
    
    os << leg_py;    
    os << "\tfig.canvas.draw()\n";
    os << "\tfig.savefig(\"" << script_gen_filename << "\",bbox_inches=\"tight\")";
  }
  void write(const std::string &filename, const std::string &script_gen_filename = "plot.pdf"){
    std::ofstream of(filename.c_str());
    of.precision(16);
    write(of,script_gen_filename);
    of.close();
  }

  void setLegend(const handleType &handle, const std::string &to){
    leg_handles.push_back(handle);
    legends.push_back(to);    
  }
  inline void createLegend(){
    kwargsType kwargs;
    createLegend(kwargs);
  }
  
  void createLegend(kwargsType kwargs){
    if(legends.size() == 0) return;
    assert(leg_handles.size() == legends.size());
    std::vector<std::string> handles_str;
    
    for(int i=0;i<leg_handles.size();i++){
      int set_idx = leg_handles[i].first;
      if(leg_handles[i].second == ErrorBandType){
	const std::string pyhandle = "plot_" + ploterrorband_sets[set_idx].tag();
	handles_str.push_back(pyhandle);
      }else if(leg_handles[i].second == HistogramType){
	const std::string pyhandle = "plot_" + plothistogram_sets[set_idx].tag() + "[2][0]";
	handles_str.push_back(pyhandle);	
      }else{
	const std::string pyhandle = "plot_" + plotdata_sets[set_idx].tag() + "[0]";
	handles_str.push_back(pyhandle);
      }
    }

    std::ostringstream py;
    py << "\tphandles = " << ListPrint<std::string>(handles_str) << '\n';
    py << "\tplegends = " << ListPrint<std::string>(legends,"\"") << '\n';
    
    if(!kwargs.count("numpoints"))
      kwargs["numpoints"] = 1;
    if(!kwargs.count("shadow"))
      kwargs["shadow"] = true;
    if(!kwargs.count("fancybox"))
      kwargs["fancybox"] = true;

    py << "\tax.legend(phandles,plegends" << kwargsPrint(kwargs) << ")\n";
    leg_py = py.str();
  }
  
  void setXlabel(const std::string &to, const std::string size = "x-large"){
    invoke() << "\tax.set_xlabel(r'" << to << "',size='" << size << "')\n";
  }
  void setYlabel(const std::string &to, const std::string size = "x-large"){
    invoke() << "\tax.set_ylabel(r'" << to << "',size='" << size << "')\n";
  }
  void setXaxisBounds(const double min, const double max){
    invoke() << "\tax.set_xlim(" << min << "," << max << ")\n";
  }
  void setYaxisBounds(const double min, const double max){
    invoke() << "\tax.set_ylim(" << min << "," << max << ")\n";
  }
  void setXaxisMajorTickSpacing(const double val){
    invoke() << "\tax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator("<<val<<"))" << std::endl;
  }
  void setYaxisMajorTickSpacing(const double val){
    invoke() << "\tax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator("<<val<<"))" << std::endl;
  }
  void setXaxisMinorTickSpacing(const double val){
    invoke() << "\tax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator("<<val<<"))" << std::endl;
  }
  void setYaxisMinorTickSpacing(const double val){
    invoke() << "\tax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator("<<val<<"))" << std::endl;
  }


};

SARLAC_END_NAMESPACE
#endif
