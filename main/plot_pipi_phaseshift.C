//Plot Luscher's curves for the pipi-scattering phase shift overlaying phenomenological curves from Schenk et al

#include<parser.h>
#include<plot.h>
#include<physics/luscherzeta.h>
#include<physics/delta_0_pheno.h>

using namespace CPSfit;

struct DiscontinuousCurve{
  std::vector<std::vector<double> > xvalues;
  std::vector<std::vector<double> > yvalues;
  
  int period;
  double prev;
  bool first;

  DiscontinuousCurve(): period(0), xvalues(1), yvalues(1), first(true){}

  void add(const double x, const double y){
    if(!first && fabs(y) > 50 && y/prev < 0){
      xvalues.resize(xvalues.size()+1);
      yvalues.resize(yvalues.size()+1);
      ++period;
    }
    xvalues[period].push_back(x);
    yvalues[period].push_back(y);
    first = false;
    prev = y;
  }
  
  int periods() const{ return xvalues.size(); }
  
  //For plotting
  void setPeriod(const int i){ period = i; }
  
  double x(const int i) const{
    return xvalues[period][i];
  }
  double y(const int i) const{
    return yvalues[period][i];
  }
  double dxm(const int i) const{ return 0; }
  double dxp(const int i) const{ return 0; }
  double dym(const int i) const{ return 0; }
  double dyp(const int i) const{ return 0; }
  int size() const{ return xvalues[period].size(); }
};


struct LuscherCurve: public DiscontinuousCurve{
  LuscherCurve(const double mpi, const int L, const double Emin, const double Emax, const int npoints, const std::vector<int> &twists){
    LuscherZeta zeta({twists[0],twists[1],twists[2]}, {0,0,0});
    double dE = (Emax - Emin)/(npoints-1);
    for(int i=0;i<npoints;i++){
      double E = Emin + i*dE;
      double delta = phaseShiftZ(zeta,E,mpi,L);
      this->add(E, delta);
    }
  }
};

struct SchenkCurve: public DiscontinuousCurve{
  SchenkCurve(const int I, const char curve, const double mpi_latt, const double ainv_GeV, const double Emin_latt, const double Emax_latt, const int npoints){
    double ainv_MeV = ainv_GeV * 1000;
    double mpi_MeV = ainv_MeV * mpi_latt;

    double dE = (Emax_latt - Emin_latt)/(npoints-1);
    
    for(int i=0;i<npoints;i++){
      double E_latt = Emin_latt + i*dE;
      double s_MeV2 = pow(E_latt * ainv_MeV,2);

      this->add(E_latt, PhenoCurveSchenk::compute(s_MeV2, I, mpi_MeV, curve) * 180./M_PI);
    }    
  }
};

//ainv in GeV
//other quantities in lattice units
#define ARGS_MEMBERS \
  ( std::vector<int>, twists )			\
  ( int, L )					\
  ( double, ainv)				\
  ( double, mpi )				\
  ( double, Emin)				\
  ( double, Emax)				\
  ( int, npoints)

struct Args{
  GENERATE_MEMBERS(ARGS_MEMBERS);

  Args(): twists({0,0,0}), L(32), ainv(1.3784), mpi(0.10382), Emin(0.3), Emax(1.), npoints(200){}
};
GENERATE_PARSER(Args, ARGS_MEMBERS);


int main(const int argc, const char *argv[]){
  Args args;
  if(argc < 2){
    std::ofstream of("template.args");
    (std::cout << "No parameter file provided: writing template to 'template.args' and exiting\n").flush();
    of << args;
    return 1;
  }    
  
  parse(args, argv[1]);

  LuscherCurve luscher(args.mpi,args.L,args.Emin,args.Emax,args.npoints,args.twists);
  LuscherCurve luscher_n1(luscher);
  LuscherCurve luscher_nm1(luscher);
  for(int i=0;i<luscher.yvalues.size();i++)
    for(int j=0;j<luscher.yvalues[i].size();j++){
      luscher_n1.yvalues[i][j] += 180;
      luscher_nm1.yvalues[i][j] -= 180;
    }

  SchenkCurve schenk_I0_A(0, 'A', args.mpi, args.ainv, args.Emin, args.Emax,args.npoints);
  SchenkCurve schenk_I0_B(0, 'B', args.mpi, args.ainv, args.Emin, args.Emax,args.npoints);
  SchenkCurve schenk_I0_C(0, 'C', args.mpi, args.ainv, args.Emin, args.Emax,args.npoints);

  SchenkCurve schenk_I2_A(2, 'A', args.mpi, args.ainv, args.Emin, args.Emax,args.npoints);
  SchenkCurve schenk_I2_B(2, 'B', args.mpi, args.ainv, args.Emin, args.Emax,args.npoints);
  SchenkCurve schenk_I2_C(2, 'C', args.mpi, args.ainv, args.Emin, args.Emax,args.npoints);

  SchenkCurve* schenk[2][3] = { {&schenk_I0_A,&schenk_I0_B,&schenk_I0_C},
				{&schenk_I2_A,&schenk_I2_B,&schenk_I2_C} };


  MatPlotLibScriptGenerate plot;
  MatPlotLibScriptGenerate::kwargsType kwargs;
  kwargs["color"] = "r";
  kwargs["linestyle"] = "-";
  kwargs["linewidth"] = 1.5;
  kwargs["marker"] = "None";

  for(int p=0;p<luscher.periods();p++){
    luscher.setPeriod(p);
    plot.plotData(luscher, kwargs);
  }
  for(int p=0;p<luscher.periods();p++){
    luscher_n1.setPeriod(p);
    plot.plotData(luscher_n1, kwargs);
  }
  for(int p=0;p<luscher.periods();p++){
    luscher_nm1.setPeriod(p);
    plot.plotData(luscher_nm1, kwargs);
  }

  
  for(int i=0;i<2;i++){
    kwargs["color"] = i == 0 ? "b" : "g";

    for(int c=0;c<3;c++){
      for(int p=0;p<schenk[i][c]->periods();p++){
	schenk[i][c]->setPeriod(p);
	plot.plotData(*schenk[i][c], kwargs);
      }
    }
  }

  plot.setYaxisBounds(-180,180);
  plot.setXaxisMajorTickSpacing(0.1);
  plot.setXaxisMinorTickSpacing(0.02);
  plot.setYaxisMajorTickSpacing(20);
  plot.setYaxisMinorTickSpacing(5);
  plot.setXlabel(R"($E_{\pi\pi}^{\rm latt}$)");
  plot.setYlabel(R"($\delta$ (deg))");

  plot.write("phase_shift.py","phase_shift.pdf");
  
  return 0;
}
