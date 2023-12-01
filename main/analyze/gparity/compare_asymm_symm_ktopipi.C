#if 1

//Code needs to be updated and fixed

int main(void){
  return 0;
}


#else

#include <ktopipi_common/ktopipi_common.h>

using namespace CPSfit;

#include <compare_asymm_symm_ktopipi/cmdline.h>
#include <compare_asymm_symm_ktopipi/args.h>
#include <compare_simple_correlators/compare.h>

void errorWeightedAverage(std::vector<correlationFunction<double, jackknifeDistributionD> > &out_asymm,
			  std::vector<correlationFunction<double, jackknifeDistributionD> > &out_symm,
			  const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &in_asymm,
			  const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &in_symm,
			  const int tmin_k_op){
  assert(in_asymm[0].size() > 0);
  const int nsample = in_asymm[0].value(0).size();
  
  for(int q=0;q<10;q++){
    std::map<int, std::vector<int> > dmap;
    for(int d=0;d<in_asymm[q].size();d++){
      const int tsep_k_op = int(in_asymm[q].coord(d).t);
      const int tsep_k_pi = in_asymm[q].coord(d).tsep_k_pi;
      const int tsep_op_pi = tsep_k_pi - tsep_k_op;
      if(tsep_op_pi >= 0 && tsep_k_op >= tmin_k_op) dmap[tsep_op_pi].push_back(d);
    }
    for(auto it = dmap.begin(); it != dmap.end(); it++){
      const int tsep_op_pi = it->first;
      const std::vector<int> &include = it->second;

      double wsum = 0.;
      jackknifeDistributionD wavg_asymm(nsample,0.);
      jackknifeDistributionD wavg_symm(nsample,0.);

      for(int dd=0; dd<include.size();dd++){
	const int d = include[dd];
	const double stderr = in_asymm[q].value(d).standardError(); //use weights from asymm for both for consistent combination
	const double w = 1/stderr/stderr;
	wsum = wsum + w;
	wavg_asymm = wavg_asymm + w*in_asymm[q].value(d);
	wavg_symm = wavg_symm + w*in_symm[q].value(d);
      }
      wavg_asymm = wavg_asymm / wsum;
      wavg_symm = wavg_symm / wsum;

      out_asymm[q].push_back(tsep_op_pi, wavg_asymm);
      out_symm[q].push_back(tsep_op_pi, wavg_symm);
    }
  }
  
}


void getRawData(type1234Data &type1, type1234Data &type2, type1234Data &type3, type1234Data &type4, const int tsep_k_pi_idx,
		const ComparisonArgs &args, const ComparisonCMDline &cmdline, const bool symmetric){
  int tsep_k_pi = args.tsep_k_pi[tsep_k_pi_idx];
  if(cmdline.load_data_checkpoint){
    std::ostringstream file; file << (symmetric ? cmdline.load_data_checkpoint_symm_stub :  cmdline.load_data_checkpoint_asymm_stub ) << "_tsepkpi" << args.tsep_k_pi[tsep_k_pi_idx] << ".hdf5";
    std::cout << "Loading checkpoint data for tsep_k_pi = " << args.tsep_k_pi[tsep_k_pi_idx] << " from " << file.str() << std::endl;
    HDF5reader rd(file.str());
    read(rd,type1,"type1");
    read(rd,type2,"type2");
    read(rd,type3,"type3");
    read(rd,type4,"type4");
  }else{
    const std::string &data_dir = symmetric ? args.data_dir_symm : args.data_dir_asymm;

    std::vector<std::string> data_file_fmt_asymm =
      { "traj_<TRAJ>_type1_deltat_<TSEP_K_PI>_sep_<TSEP_PIPI>_mom<MOM>",
	"traj_<TRAJ>_type2_deltat_<TSEP_K_PI>_sep_<TSEP_PIPI>",
	"traj_<TRAJ>_type3_deltat_<TSEP_K_PI>_sep_<TSEP_PIPI>",
	"traj_<TRAJ>_type4" };

    std::vector<std::string> data_file_fmt_symm(data_file_fmt_asymm);
    for(int i=0;i<data_file_fmt_symm.size();i++) data_file_fmt_symm[i] = data_file_fmt_symm[i] + cmdline.symmetric_quark_momenta_figure_file_extension;

    
    const std::vector<std::string> &data_file_fmt = symmetric ? data_file_fmt_symm : data_file_fmt_asymm;


    std::vector<std::pair<threeMomentum, double> > type1_pimom_proj = {  { {1,1,1}, 1.0/8.0 }, { {-1,-1,-1}, 1.0/8.0 },  { {-1,1,1}, 3.0/8.0 }, { {1,-1,-1}, 3.0/8.0 }  };

    KtoPiPiFilenamePolicyGen fp(data_file_fmt[0], data_file_fmt[1], data_file_fmt[2], data_file_fmt[3]);
    BasicKtoPiPiReadPolicy<KtoPiPiFilenamePolicyGen> rp(data_dir, args.traj_start, args.traj_inc, args.traj_lessthan, fp);
    type1 = readKtoPiPiType(1, tsep_k_pi, args.tsep_pipi, args.Lt, type1_pimom_proj, rp);
    type2 = readKtoPiPiType(2, tsep_k_pi, args.tsep_pipi, args.Lt, type1_pimom_proj, rp);
    type3 = readKtoPiPiType(3, tsep_k_pi, args.tsep_pipi, args.Lt, type1_pimom_proj, rp);
    type4 = readKtoPiPiType(4, tsep_k_pi, args.tsep_pipi, args.Lt, type1_pimom_proj, rp);
  }

  if(cmdline.save_data_checkpoint){
    std::ostringstream file; file << (symmetric ? cmdline.save_data_checkpoint_symm_stub : cmdline.save_data_checkpoint_asymm_stub) << "_tsepkpi" << args.tsep_k_pi[tsep_k_pi_idx] << ".hdf5";
    std::cout << "Saving checkpoint data for tsep_k_pi = " << tsep_k_pi << " to " << file.str() << std::endl;
    HDF5writer wr(file.str());
    write(wr,type1,"type1");
    write(wr,type2,"type2");
    write(wr,type3,"type3");
    write(wr,type4,"type4");
  }
}

void compare(const NumericTensor<jackknifeDistributionD,1> &asymm,
	      const NumericTensor<jackknifeDistributionD,1> &symm,
	      const int tsep_k_pi,
	      const std::string &descr){
  std::cout << descr << ": t  asymm  symm  diff  rel.diff.\n";
  for(int t=0;t<tsep_k_pi;t++){
    jackknifeDistributionD diff = symm(&t) - asymm(&t);
    jackknifeDistributionD reldiff = 2.*(symm(&t)-asymm(&t))/(symm(&t)+asymm(&t));
    std::cout << t << " " << asymm(&t) << " " << symm(&t) << " " << " " << diff << " " << reldiff << std::endl;
  }
}


//Read and prepare the data for a particular tsep_k_pi_idx
void analyze(const ProjectedBubbleData &bubble_data_asymm,
	     const ProjectedBubbleData &bubble_data_symm,
	     const int tsep_k_pi_idx, const ComparisonArgs &args, const ComparisonCMDline &cmdline){
  basic_bin_resampler resampler(args.bin_size);

  std::cout << "Comparing data for tsep_k_pi = " <<  args.tsep_k_pi[tsep_k_pi_idx] << std::endl;
  int tsep_k_pi = args.tsep_k_pi[tsep_k_pi_idx];

  type1234Data type1_asymm, type2_asymm, type3_asymm, type4_asymm;
  type1234Data type1_symm, type2_symm, type3_symm, type4_symm;

  getRawData(type1_asymm, type2_asymm, type3_asymm, type4_asymm, tsep_k_pi_idx, args, cmdline, false);
  getRawData(type1_symm, type2_symm, type3_symm, type4_symm, tsep_k_pi_idx, args, cmdline, true);
  
  std::vector<int> type1_nonzerotK = type1_asymm.getNonZeroKaonTimeslices();
  std::vector<int> type2_nonzerotK = type2_asymm.getNonZeroKaonTimeslices();
  std::vector<int> type3_nonzerotK = type3_asymm.getNonZeroKaonTimeslices();
  std::vector<int> type4_nonzerotK = type4_asymm.getNonZeroKaonTimeslices();

  typedef computeAmplitudeAlltKtensorControls ctrl;
  typedef computeAmplitudeAlltKtensorControls::inputType ctrli;
  
  //Get the type1-4 and mix3/mix4 components of the data.
  //Note that we have not yet multiplied the type4/mix4 data by the pipi bubble, hence these  are just the K->vac component which are needed for computing alpha
  std::cout << "Computing raw data diagrams prior to multiplying by bubble\n";
  NumericTensor<rawDataDistributionD,3> A0_type1_alltK_asymm = computeAmplitudeType1<ctrl>(ctrli(type1_asymm, type1_nonzerotK, tsep_k_pi)  ); //[Qidx][tK][t]
  NumericTensor<rawDataDistributionD,3> A0_type2_alltK_asymm = computeAmplitudeType2<ctrl>(ctrli(type2_asymm, type2_nonzerotK, tsep_k_pi) ); //[Qidx][tK][t]
  NumericTensor<rawDataDistributionD,3> A0_type3_alltK_asymm = computeAmplitudeType3<ctrl>(ctrli(type3_asymm, type3_nonzerotK, tsep_k_pi) ); //[Qidx][tK][t]
  NumericTensor<rawDataDistributionD,3> A0_type4_alltK_nobub_asymm = computeAmplitudeType4<ctrl>(ctrli(type4_asymm, type4_nonzerotK, tsep_k_pi) ); //[Qidx][tK][t]
  
  NumericTensor<rawDataDistributionD,3> A0_type1_alltK_symm = computeAmplitudeType1<ctrl>(ctrli(type1_symm, type1_nonzerotK, tsep_k_pi) );
  NumericTensor<rawDataDistributionD,3> A0_type2_alltK_symm = computeAmplitudeType2<ctrl>(ctrli(type2_symm, type2_nonzerotK, tsep_k_pi) );
  NumericTensor<rawDataDistributionD,3> A0_type3_alltK_symm = computeAmplitudeType3<ctrl>(ctrli(type3_symm, type3_nonzerotK, tsep_k_pi) );
  NumericTensor<rawDataDistributionD,3> A0_type4_alltK_nobub_symm = computeAmplitudeType4<ctrl>(ctrli(type4_symm, type4_nonzerotK, tsep_k_pi) );

  NumericTensor<rawDataDistributionD,2> mix3_alltK_asymm({args.Lt,args.Lt}), mix3_alltK_symm({args.Lt,args.Lt}); //[tK][t]
  NumericTensor<rawDataDistributionD,2> mix4_alltK_nobub_asymm({args.Lt,args.Lt}), mix4_alltK_nobub_symm({args.Lt,args.Lt}); //[tK][t]

  for(int tK=0;tK<args.Lt;tK++)
    for(int t=0;t<args.Lt;t++){
      mix3_alltK_asymm({tK,t}) = type3_asymm(tK,t).mix();
      mix4_alltK_nobub_asymm({tK,t}) = type4_asymm(tK,t).mix();

      mix3_alltK_symm({tK,t}) = type3_symm(tK,t).mix();
      mix4_alltK_nobub_symm({tK,t}) = type4_symm(tK,t).mix();
    }

  type1_asymm.freeData(); type2_asymm.freeData(); type3_asymm.freeData(); type4_asymm.freeData();
  type1_symm.freeData(); type2_symm.freeData(); type3_symm.freeData(); type4_symm.freeData();


  //Compute the type4/mix4 data with the bubble included
  std::cout << "Computing raw type4/mix4 data with bubble included" << std::endl;
  NumericTensor<rawDataDistributionD,3> A0_type4_alltK_asymm({10,args.Lt,args.Lt}), A0_type4_alltK_symm({10,args.Lt,args.Lt});
  NumericTensor<rawDataDistributionD,2> mix4_alltK_asymm({args.Lt,args.Lt}), mix4_alltK_symm({args.Lt,args.Lt});

  for(int tK=0;tK<args.Lt;tK++)
    for(int t=0;t<args.Lt;t++){
      int tB = (tK + tsep_k_pi + args.tsep_pipi) % args.Lt;
      mix4_alltK_asymm({tK,t}) = mix4_alltK_nobub_asymm({tK,t})*bubble_data_asymm.bubble(&tB);
      mix4_alltK_symm({tK,t}) = mix4_alltK_nobub_symm({tK,t})*bubble_data_symm.bubble(&tB);

      for(int q=0;q<10;q++){
	A0_type4_alltK_asymm({q,tK,t}) = A0_type4_alltK_nobub_asymm({q,tK,t})*bubble_data_asymm.bubble(&tB);
	A0_type4_alltK_symm({q,tK,t}) = A0_type4_alltK_nobub_symm({q,tK,t})*bubble_data_symm.bubble(&tB);
      }
    }

  //Bin everything we are going to use henceforth
  A0_type1_alltK_asymm = bin(A0_type1_alltK_asymm, args.bin_size);
  A0_type2_alltK_asymm = bin(A0_type2_alltK_asymm, args.bin_size);
  A0_type3_alltK_asymm = bin(A0_type3_alltK_asymm, args.bin_size);
  A0_type4_alltK_asymm = bin(A0_type4_alltK_asymm, args.bin_size);
  mix3_alltK_asymm = bin(mix3_alltK_asymm, args.bin_size);
  mix4_alltK_asymm = bin(mix4_alltK_asymm, args.bin_size);
  A0_type4_alltK_nobub_asymm = bin(A0_type4_alltK_nobub_asymm, args.bin_size);
  mix4_alltK_nobub_asymm = bin(mix4_alltK_nobub_asymm, args.bin_size);

  A0_type1_alltK_symm = bin(A0_type1_alltK_symm, args.bin_size);
  A0_type2_alltK_symm = bin(A0_type2_alltK_symm, args.bin_size);
  A0_type3_alltK_symm = bin(A0_type3_alltK_symm, args.bin_size);
  A0_type4_alltK_symm = bin(A0_type4_alltK_symm, args.bin_size);
  mix3_alltK_symm = bin(mix3_alltK_symm, args.bin_size);
  mix4_alltK_symm = bin(mix4_alltK_symm, args.bin_size);
  A0_type4_alltK_nobub_symm = bin(A0_type4_alltK_nobub_symm, args.bin_size);
  mix4_alltK_nobub_symm = bin(mix4_alltK_nobub_symm, args.bin_size);

  for(int q=0;q<10;q++){
    std::cout << "Starting Q" << q+1 << std::endl;

    printMem("Starting new Q");
    
    //Compute alpha and type4/mix4 vacuum subtractions
    std::cout << "Computing single-jackknife alpha and vacuum subtractions\n";
    NumericTensor<jackknifeDistributionD,1> alpha_j_asymm({args.Lt}), alpha_j_symm({args.Lt}); //[t]
    NumericTensor<jackknifeDistributionD,1> A0_type4_srcavg_vacsub_j_asymm({args.Lt}), A0_type4_srcavg_vacsub_j_symm({args.Lt}); //[t]
    NumericTensor<jackknifeDistributionD,1> mix4_srcavg_vacsub_j_asymm({args.Lt}), mix4_srcavg_vacsub_j_symm({args.Lt}); //[t]

    computeAlphaAndVacuumSubtractions(alpha_j_asymm,A0_type4_srcavg_vacsub_j_asymm,mix4_srcavg_vacsub_j_asymm,
				      A0_type4_alltK_nobub_asymm,mix4_alltK_nobub_asymm,bubble_data_asymm.bubble_j,q,type4_nonzerotK,tsep_k_pi,args.Lt, resampler);

    computeAlphaAndVacuumSubtractions(alpha_j_symm,A0_type4_srcavg_vacsub_j_symm,mix4_srcavg_vacsub_j_symm,
				      A0_type4_alltK_nobub_symm,mix4_alltK_nobub_symm,bubble_data_symm.bubble_j,q,type4_nonzerotK,tsep_k_pi,args.Lt, resampler);

    compare(alpha_j_asymm, alpha_j_symm, tsep_k_pi, "Alpha");

    std::cout << "Computing single-jackknife tK averages and mix diagrams\n";
    NumericTensor<jackknifeDistributionD,1> A0_type1_srcavg_j_asymm = resampleAverageTypeData<jackknifeDistributionD>(A0_type1_alltK_asymm, q, type1_nonzerotK, args.Lt, resampler); //[t]
    NumericTensor<jackknifeDistributionD,1> A0_type2_srcavg_j_asymm = resampleAverageTypeData<jackknifeDistributionD>(A0_type2_alltK_asymm, q, type2_nonzerotK, args.Lt, resampler); //[t]
    NumericTensor<jackknifeDistributionD,1> A0_type3_srcavg_j_asymm = resampleAverageTypeData<jackknifeDistributionD>(A0_type3_alltK_asymm, q, type3_nonzerotK, args.Lt, resampler); //[t]
    NumericTensor<jackknifeDistributionD,1> A0_type4_srcavg_j_asymm = resampleAverageTypeData<jackknifeDistributionD>(A0_type4_alltK_asymm, q, type4_nonzerotK, args.Lt, resampler); //[t]
    NumericTensor<jackknifeDistributionD,1> mix3_srcavg_j_asymm = resampleAverageMixDiagram<jackknifeDistributionD>(mix3_alltK_asymm, type3_nonzerotK, args.Lt, resampler);
    NumericTensor<jackknifeDistributionD,1> mix4_srcavg_j_asymm = resampleAverageMixDiagram<jackknifeDistributionD>(mix4_alltK_asymm, type4_nonzerotK, args.Lt, resampler);

    NumericTensor<jackknifeDistributionD,1> A0_type1_srcavg_j_symm = resampleAverageTypeData<jackknifeDistributionD>(A0_type1_alltK_symm, q, type1_nonzerotK, args.Lt, resampler); //[t]
    NumericTensor<jackknifeDistributionD,1> A0_type2_srcavg_j_symm = resampleAverageTypeData<jackknifeDistributionD>(A0_type2_alltK_symm, q, type2_nonzerotK, args.Lt, resampler); //[t]
    NumericTensor<jackknifeDistributionD,1> A0_type3_srcavg_j_symm = resampleAverageTypeData<jackknifeDistributionD>(A0_type3_alltK_symm, q, type3_nonzerotK, args.Lt, resampler); //[t]
    NumericTensor<jackknifeDistributionD,1> A0_type4_srcavg_j_symm = resampleAverageTypeData<jackknifeDistributionD>(A0_type4_alltK_symm, q, type4_nonzerotK, args.Lt, resampler); //[t]
    NumericTensor<jackknifeDistributionD,1> mix3_srcavg_j_symm = resampleAverageMixDiagram<jackknifeDistributionD>(mix3_alltK_symm, type3_nonzerotK, args.Lt, resampler);
    NumericTensor<jackknifeDistributionD,1> mix4_srcavg_j_symm = resampleAverageMixDiagram<jackknifeDistributionD>(mix4_alltK_symm, type4_nonzerotK, args.Lt, resampler);

    compare(A0_type1_srcavg_j_asymm, A0_type1_srcavg_j_symm, tsep_k_pi, "Type 1");
    compare(A0_type2_srcavg_j_asymm, A0_type2_srcavg_j_symm, tsep_k_pi, "Type 2");
    compare(A0_type3_srcavg_j_asymm, A0_type3_srcavg_j_symm, tsep_k_pi, "Type 3");
    compare(A0_type4_srcavg_j_asymm, A0_type4_srcavg_j_symm, tsep_k_pi, "Type 4");

    //Subtract the pseudoscalar operators and mix4 vacuum term
    std::cout << "Subtracting pseudoscalar operators and mix4 vacuum term under single-jackknife\n";
    A0_type3_srcavg_j_asymm = A0_type3_srcavg_j_asymm.transform([&](int const* t, const jackknifeDistributionD &from){ return jackknifeDistributionD(from - alpha_j_asymm(t)*mix3_srcavg_j_asymm(t)); }); 
    A0_type4_srcavg_j_asymm = A0_type4_srcavg_j_asymm.transform(
						    [&](int const* t, const jackknifeDistributionD &from){
						      return jackknifeDistributionD(from - alpha_j_asymm(t)*( mix4_srcavg_j_asymm(t) - mix4_srcavg_vacsub_j_asymm(t) ) );
						    }
						    );

    A0_type3_srcavg_j_symm = A0_type3_srcavg_j_symm.transform([&](int const* t, const jackknifeDistributionD &from){ return jackknifeDistributionD(from - alpha_j_symm(t)*mix3_srcavg_j_symm(t)); }); 
    A0_type4_srcavg_j_symm = A0_type4_srcavg_j_symm.transform(
						    [&](int const* t, const jackknifeDistributionD &from){
						      return jackknifeDistributionD(from - alpha_j_symm(t)*( mix4_srcavg_j_symm(t) - mix4_srcavg_vacsub_j_symm(t) ) );
						    }
						    );
    compare(A0_type3_srcavg_j_asymm, A0_type3_srcavg_j_symm, tsep_k_pi, "Type 3 post mix sub");
    compare(A0_type4_srcavg_j_asymm, A0_type4_srcavg_j_symm, tsep_k_pi, "Type 4 post vac-subbed mix sub");

    //Perform the type 4 vacuum subtraction
    std::cout << "Performing type-4 vacuum subtraction\n";
    A0_type4_srcavg_j_asymm = A0_type4_srcavg_j_asymm - A0_type4_srcavg_vacsub_j_asymm;
    A0_type4_srcavg_j_symm = A0_type4_srcavg_j_symm - A0_type4_srcavg_vacsub_j_symm;

    compare(A0_type4_srcavg_j_asymm, A0_type4_srcavg_j_symm, tsep_k_pi, "Type 4 post vac-sub");

    //Get the full double-jackknife amplitude
    std::cout << "Computing full amplitudes\n";
    NumericTensor<jackknifeDistributionD,1> A0_full_srcavg_j_asymm = A0_type1_srcavg_j_asymm + A0_type2_srcavg_j_asymm + A0_type3_srcavg_j_asymm + A0_type4_srcavg_j_asymm;
    NumericTensor<jackknifeDistributionD,1> A0_full_srcavg_j_symm = A0_type1_srcavg_j_symm + A0_type2_srcavg_j_symm + A0_type3_srcavg_j_symm + A0_type4_srcavg_j_symm;

    compare(A0_full_srcavg_j_asymm, A0_full_srcavg_j_symm, tsep_k_pi, "Full amplitude");
  }

}

void analyze(const ComparisonArgs &args, const ComparisonCMDline &cmdline){
  //Read the bubble data
  basic_resampler resampler;

  std::string bubble_file_fmt_asymm = "traj_<TRAJ>_FigureVdis_sep<TSEP_PIPI>_mom<PB>";
  std::string bubble_file_fmt_symm =  bubble_file_fmt_asymm + "_symm";

  std::vector<std::pair<threeMomentum, double> > bubble_pimom_proj =  {  { {1,1,1}, 1.0/8.0 }, { {-1,-1,-1}, 1.0/8.0 },  
								         { {-1,1,1}, 1.0/8.0 }, { {1,-1,-1}, 1.0/8.0 }, 
								         { {1,-1,1}, 1.0/8.0 }, { {-1,1,-1}, 1.0/8.0 }, 
								         { {1,1,-1}, 1.0/8.0 }, { {-1,-1,1}, 1.0/8.0 } };

  readKtoPiPiDataOptions opt_asymm = cmdline.getReadOptions(Asymmetric);
  ProjectedBubbleData bubble_data_asymm(args.data_dir_asymm, bubble_file_fmt_asymm, args.traj_start, args.traj_inc, args.traj_lessthan, args.bin_size, 
					args.Lt, args.tsep_pipi, bubble_pimom_proj, resampler, opt_asymm);

  readKtoPiPiDataOptions opt_symm = cmdline.getReadOptions(Symmetric);
  ProjectedBubbleData bubble_data_symm(args.data_dir_symm, bubble_file_fmt_symm, args.traj_start, args.traj_inc, args.traj_lessthan, args.bin_size,
				       args.Lt, args.tsep_pipi, bubble_pimom_proj, resampler, opt_symm);
 
  for(int tsep_k_pi_idx=0;tsep_k_pi_idx<args.tsep_k_pi.size();tsep_k_pi_idx++){
    analyze(bubble_data_asymm,bubble_data_symm,tsep_k_pi_idx,args,cmdline);  
  }	
}

int main(const int argc, const char* argv[]){
  ComparisonArgs args;
  if(argc < 2){
    std::ofstream of("template.args");
    of << args;
    std::cout << "Wrote template argument file to template.args\n";
    return 0;
  }
  const std::string arg_file = argv[1];
  parse(args, arg_file);
  
  assert( (args.traj_lessthan - args.traj_start) % args.traj_inc == 0 );  

  ComparisonCMDline cmdline(argc,argv,2);
 
  analyze(args,cmdline);

  std::cout << "Done" << std::endl;
  
  return 0;
}

#endif
