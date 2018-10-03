#ifndef _CPSFIT_MPL_SCRIPT_GENERATE_3D_H_
#define _CPSFIT_MPL_SCRIPT_GENERATE_3D_H_

//Generate python scripts for generating 3D plots using Matplotlib

#include<iostream>
#include<vector>

#include<config.h>
#include<utils/macros.h>
#include<plot/plot/write_python.h>
#include<plot/plot/datasets_3d.h>
#include<plot/plot/MPLscriptgenerate_base.h>

CPSFIT_START_NAMESPACE

//A python script generator for 3D plots
class MatPlotLib3DscriptGenerate: public MatPlotLibScriptGenerateBase{
  std::vector<Python3DdataContainer> wireframe_sets;
  std::vector<kwargsType> wireframe_args;

  std::vector<Python3DdataContainer> scatter_sets;
  std::vector<kwargsType> scatter_args;  
public:
  typedef std::pair<int,SetType3D> handleType;
  typedef MatPlotLibScriptGenerateBase::kwargsType kwargsType;

  template<typename Accessor>
  handleType plotWireframe(const Accessor &data, const kwargsType &kwargs){
    std::ostringstream os; os << "wireframe" << wireframe_sets.size();
    wireframe_sets.push_back(Python3DdataContainer());
    wireframe_sets.back().import(data, os.str());
    wireframe_args.push_back(kwargs);
    return handleType(wireframe_sets.size()-1, WireframeType);
  }
  template<typename Accessor>
  inline handleType plotWireframe(const Accessor &data){
    kwargsType kwargs;
    return plotWireframe(data,kwargs);
  }

  template<typename Accessor>
  handleType plotScatter(const Accessor &data, const kwargsType &kwargs){
    std::ostringstream os; os << "scatter" << scatter_sets.size();
    scatter_sets.push_back(Python3DdataContainer());
    scatter_sets.back().import(data, os.str());
    scatter_args.push_back(kwargs);
    return handleType(scatter_sets.size()-1, ScatterType);
  }
  template<typename Accessor>
  inline handleType plotScatter(const Accessor &data){
    kwargsType kwargs;
    return plotScatter(data,kwargs);
  }

  
  void write(std::ostream &os, const std::string &script_gen_filename = "plot.pdf") const{
    os << "import pyplot\n\n";

    for(int i=0;i<wireframe_sets.size();i++)
      wireframe_sets[i].write(os);

    for(int i=0;i<scatter_sets.size();i++)
      scatter_sets[i].write(os);

    os << preamble.str(); //user code    

    os << "\nif __name__ == '__main__':\n";
    os << "\tfig = pyplot.plt.figure()\n";
    os << "\tax = fig.add_subplot(1,1,1, projection='3d')\n";

    for(int i=0;i<wireframe_sets.size();i++)
      os << "\tplot_" << wireframe_sets[i].tag() << " = " << "pyplot.plotWireframe(ax, " << wireframe_sets[i].tag() << kwargsPrint(wireframe_args[i]) << ")\n";
    
    for(int i=0;i<scatter_sets.size();i++)
      os << "\tplot_" << scatter_sets[i].tag() << " = " << "pyplot.plotScatter(ax, " << scatter_sets[i].tag() << kwargsPrint(scatter_args[i]) << ")\n";

    os << user.str(); //user code
    
    os << "\tfig.canvas.draw()\n";
    os << "\tfig.savefig(\"" << script_gen_filename << "\")";
  }
  void write(const std::string &filename, const std::string &script_gen_filename = "plot.pdf"){
    std::ofstream of(filename.c_str());
    write(of,script_gen_filename);
    of.close();
  }
  
};

CPSFIT_END_NAMESPACE
#endif
