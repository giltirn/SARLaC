#if 1

//Fixme

int main(void){
  return 0;
}
#else


#include <pipi_common/pipi_common.h>
#include <pipi_common/analyze_chisq.h>

using namespace CPSfit;

#include <fit_pipi_gparity/args.h>
#include <fit_pipi_gparity/cmdline.h>

#include <fit_pipi_sigma_sim_gparity/fitfunc.h>
#include <fit_pipi_sigma_sim_gparity/fit.h>
#include <fit_pipi_sigma_sim_gparity/args.h>
#include <fit_pipi_sigma_sim_gparity/cmdline.h>
#include <fit_pipi_sigma_sim_gparity/read_data.h>


template<typename DistributionType>
void getResampledData(correlationFunction<SimFitCoord, DistributionType> &corr_comb, const rawData &raw_data,
		      const int Lt, const int tsep_pipi, const int bin_size, bool do_vacuum_subtraction,
		      bool include_pipi_2pt, bool include_pipi_to_sigma, bool include_sigma_2pt){
  typedef correlationFunction<double,DistributionType> CorrFuncType;

  //Get double-jack data
  auto pipi_to_sigma = binResample<CorrFuncType>(raw_data.pipi_to_sigma_raw, bin_size);
  auto sigma2pt = binResample<CorrFuncType>(raw_data.sigma2pt_raw, bin_size);
  auto pipi = binResample<CorrFuncType>(raw_data.pipi_raw, bin_size);
  
  //Compute vacuum subtractions
  if(do_vacuum_subtraction){
    { //Pipi->sigma
      PiPiProjectorA1Basis111 proj_pipi;
      bubbleData pipi_self_proj = projectSourcePiPiBubble(raw_data.pipi_self_data, proj_pipi);
      pipi_to_sigma = pipi_to_sigma - computePiPiToSigmaVacSub<CorrFuncType>(raw_data.sigma_self_data, pipi_self_proj, bin_size);
    }
    { //sigma->sigma
      sigma2pt = sigma2pt - computeSigmaVacSub<CorrFuncType>(raw_data.sigma_self_data, bin_size);
    }
    { //Pipi->pipi
      pipi = pipi - computePiPi2ptVacSub<CorrFuncType>(raw_data.pipi_self_data, bin_size, tsep_pipi);
    }
  }
  
  //Fold data
  pipi_to_sigma = fold(pipi_to_sigma, tsep_pipi);

  sigma2pt = fold(sigma2pt, 0);

  pipi = fold(pipi, 2*tsep_pipi);
  
  //Build the combined data set
  CorrFuncType const* dsets[3] = { &pipi, &pipi_to_sigma, &sigma2pt };

  static const SimFitType dtype[3] = { SimFitType::PiPi2pt, SimFitType::PiPiToSigma, SimFitType::Sigma2pt };
  bool dincl[3] = { include_pipi_2pt, include_pipi_to_sigma, include_sigma_2pt };

  for(int d=0;d<3;d++)
    if(dincl[d])
      for(int t=0; t<Lt; t++){ 
	corr_comb.push_back(SimFitCoord(dtype[d],t), dsets[d]->value(t));
      }
}		      
		      
typedef bootstrapDistribution<double> bootstrapDistributionD;
typedef correlationFunction<SimFitCoord,  bootstrapDistributionD> simFitCorrFuncB;

int main(const int argc, const char* argv[]){
  PiPiSigmaSimArgs args;
  if(argc < 2){
    std::ofstream of("template.args");
    of << args;
    std::cout << "Wrote template argument file to template.args\n";
    return 0;
  }
  {
    const std::string arg_file = argv[1];
    parse(args,arg_file);
    std::cout << "Read arguments: \n" << args << std::endl;
  }

  PiPiSigmaSimCMDline cmdline(argc,argv,2);

  bool include_type[3] = {false,false,false};
  include_type[(int)SimFitType::PiPi2pt] = cmdline.include_pipi_2pt;
  include_type[(int)SimFitType::PiPiToSigma] = cmdline.include_pipi_to_sigma;
  include_type[(int)SimFitType::Sigma2pt] = cmdline.include_sigma_2pt;

  simFitCorrFuncJ corr_comb_j_all;
  simFitCorrFuncDJ corr_comb_dj_all;

  if(cmdline.load_resampled_data){
    HDF5reader rd(cmdline.load_resampled_data_file);
    read(rd,corr_comb_j_all, "corr_comb_j_all");
    read(rd,corr_comb_dj_all, "corr_comb_dj_all");
  }else{
    rawData raw_data;
    if(cmdline.load_checkpoint){
      raw_data.readDataFromCheckpoint(cmdline.load_checkpoint_file);
    }else{
      raw_data.readDataFromOrigFiles(args.data_dir, 
				     args.pipi2pt_figure_file_fmt, args.sigma2pt_file_fmt, args.pipitosigma_file_fmt,
				     args.pipi_bubble_file_fmt, args.sigma_bubble_file_fmt,
				     args.tsep_pipi, args.tstep_pipi2pt, args.tstep_pipitosigma,
				     args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan, !cmdline.use_pipitosigma_disconn_complex_prod);
    }
    if(cmdline.save_checkpoint){
      raw_data.saveCheckpoint(cmdline.save_checkpoint_file);
    }
    getResampledData<jackknifeDistributionD>(corr_comb_j_all, raw_data, args.Lt, args.tsep_pipi, args.bin_size, args.do_vacuum_subtraction,
					     cmdline.include_pipi_2pt, cmdline.include_pipi_to_sigma, cmdline.include_sigma_2pt);
    getResampledData<doubleJackknifeDistributionD>(corr_comb_dj_all, raw_data, args.Lt, args.tsep_pipi, args.bin_size, args.do_vacuum_subtraction,
						   cmdline.include_pipi_2pt, cmdline.include_pipi_to_sigma, cmdline.include_sigma_2pt);

    if(cmdline.analyze_population_distribution){
      bootstrapDistributionOptions::defaultBoots() = 2000;
      simFitCorrFuncB corr_comb_b_all;      
      getResampledData<bootstrapDistributionD>(corr_comb_b_all, raw_data, args.Lt, args.tsep_pipi, args.bin_size, args.do_vacuum_subtraction,
					       cmdline.include_pipi_2pt, cmdline.include_pipi_to_sigma, cmdline.include_sigma_2pt);

      simFitCorrFuncB corr_comb_b_inrange;

      for(int i=0;i<corr_comb_b_all.size();i++){
	auto coord = corr_comb_b_all.coord(i);
	if(include_type[(int)coord.type] && 
	   coord.t >= args.t_min && 
	   coord.t <= args.t_max){
	  std::cout << coord.type << " " << coord.t 
		    << " mean " << corr_comb_b_all.value(i).mean() 
		    << " std.dev " << corr_comb_b_all.value(i).standardDeviation() 
		    << " skewness " << corr_comb_b_all.value(i).standardizedMoment(3)
	    	    << " excess kurtosis " << corr_comb_b_all.value(i).standardizedMoment(4) - 3. 
		    << std::endl;
	  corr_comb_b_inrange.push_back(corr_comb_b_all[i]);
	}
      }

      bootstrapDistributionD bzero(corr_comb_b_all.value(0)); zeroit(bzero);
      
      int nfit = corr_comb_b_inrange.size();
      NumericSquareMatrix<double> cov(nfit);
      NumericVector<double> sigma(nfit);
      for(int i=0;i<nfit;i++){
	for(int j=0;j<nfit;j++){
	  cov(i,j) = bootstrapDistributionD::covariance(corr_comb_b_inrange.value(i), corr_comb_b_inrange.value(j));
	}
	sigma(i) = sqrt(cov(i,i));
      }
      NumericSquareMatrix<double> cor(nfit, [&](const int i, const int j){ return cov(i,j)/sigma(i)/sigma(j); });
      
      std::vector<NumericVector<double> > cov_evecs;
      std::vector<double> cov_evals;
      symmetricMatrixEigensolve(cov_evecs, cov_evals, cov);

      //Each data point can be decomposed into the basis of these eigenvectors. We can then examine how important a given mode is in our data
      //d_i = \sum_a  c^a v_i^a
      //\sum_i v_i^b d_i = \sum_a,i  c^a v_i^a v_i^b = c^b
 
      //Also enters in chi^2
      //chi^2 =  \sum_ij  ( d_i - f_i ) C^{-1}_ij ( d_j - f_j )
      //      =  \sum_a \sum_ij  ( d_i - f_i ) v_i^a 1/lambda^a v_j^a ( d_j - f_j )
      //      =  \sum_a   ( \sum_i ( d_i - f_i ) v_i^a )^2/lambda^a
      
      //projected data    c^a = \sum_i d_i v_i^a
      
      std::vector<bootstrapDistributionD> corr_comb_b_inrange_projected(nfit);
      for(int a=0;a<nfit;a++){
	corr_comb_b_inrange_projected[a] = bzero;
	for(int i=0;i<nfit;i++){
	  corr_comb_b_inrange_projected[a] = corr_comb_b_inrange_projected[a] + cov_evecs[a](i) * corr_comb_b_inrange.value(i);
	}
	
	std::cout << "Projected data mode " << a
		  << " bs " << corr_comb_b_inrange_projected[a]
		  << " mean " << corr_comb_b_inrange_projected[a].mean() 
		  << " std.dev " << corr_comb_b_inrange_projected[a].standardDeviation() 
		  << " skewness " << corr_comb_b_inrange_projected[a].standardizedMoment(3)
		  << " excess kurtosis " << corr_comb_b_inrange_projected[a].standardizedMoment(4) - 3. 
		  << std::endl;
      }

      //Compute probability of mode in data:  alpha^a =   |c^a|^2/( \sum_a |c^a|^2 ) cf https://arxiv.org/pdf/1101.2248.pdf eq 15 and below.
      bootstrapDistributionD nrm = bzero;
      for(int a=0;a<nfit;a++)
	nrm = nrm + corr_comb_b_inrange_projected[a] * corr_comb_b_inrange_projected[a];
      
      std::cout << nrm << std::endl;

      std::cout << "Probability of mode in data:\n";
      for(int a=0;a<nfit;a++){
	bootstrapDistributionD alpha = corr_comb_b_inrange_projected[a] * corr_comb_b_inrange_projected[a] / nrm;
	std::cout << a  << " " << alpha << " (" << corr_comb_b_inrange_projected[a] << ")" << std::endl;
      }


    }

    

  }
  if(cmdline.save_resampled_data){
    HDF5writer wr(cmdline.save_resampled_data_file);
    write(wr,corr_comb_j_all, "corr_comb_j_all");
    write(wr,corr_comb_dj_all, "corr_comb_dj_all");
  }

  //Get data in fit range
 
  simFitCorrFuncJ corr_comb_j;
  simFitCorrFuncDJ corr_comb_dj;
  for(int tt=0;tt<corr_comb_j_all.size();tt++){
    auto coord = corr_comb_j_all.coord(tt);
    if(include_type[(int)coord.type] && 
       coord.t >= args.t_min && 
       coord.t <= args.t_max){
      corr_comb_j.push_back(corr_comb_j_all[tt]);
      corr_comb_dj.push_back(corr_comb_dj_all[tt]);
    }
  }

  std::cout << "Data in fit range:\n";
  for(int i=0;i<corr_comb_j.size();i++){
    std::cout << corr_comb_j.coord(i).type << " " << corr_comb_j.coord(i).t << " " << corr_comb_j.value(i) << std::endl;
  }

  //Run the fit
  SimFitArgs fargs;
  args.transfer(fargs);
  cmdline.transfer(fargs);

  fit(corr_comb_j, corr_comb_dj, fargs);

  std::cout << "Done\n";
  return 0;
}

#endif
