#ifndef _CPSFIT_DISTRIBUTION_HDF5IO_CONVENTIONAL_H_
#define _CPSFIT_DISTRIBUTION_HDF5IO_CONVENTIONAL_H_

//Here we define a conventional format similar to UKfit in which all DistributionType<T>, where T is some compound structure, are stored as std::vector<DistributionType<double> >. This enables us to write simple external programs to manipulate results

#include<distribution/distribution_hdf5io_conventional/helper.h>
#include<distribution/distribution_hdf5io_conventional/io_format.h>
#include<distribution/distribution_hdf5io_conventional/POD_io.h>
#include<distribution/distribution_hdf5io_conventional/compound_io.h>
#include<distribution/distribution_hdf5io_conventional/type_info.h>

#endif
