#ifndef _SUPERJACKKNIFE_HDF5_IO_H_
#define _SUPERJACKKNIFE_HDF5_IO_H_

#include<config.h>

#ifdef HAVE_HDF5

#include<distribution/superjackknife/class.h>
#include<distribution/distribution_hdf5io_basic.h>

CPSFIT_START_NAMESPACE

template<typename D>
inline void write(HDF5writer &writer, const superJackknifeDistribution<D> &sj, const std::string &tag){ sj.write(writer,tag); }
template<typename D>
inline void read(HDF5reader &reader, superJackknifeDistribution<D> &sj, const std::string &tag){ sj.read(reader,tag); }

template<typename T> //override for distributions with extra information
struct extraFlattenIO<superJackknifeDistribution<T> >{
  static inline void write(HDF5writer &writer, const std::vector<superJackknifeDistribution<T> > &value){
    assert(value.size() > 0);
    const superJackknifeLayout &layout = value[0].getLayout();
    if(value.size() > 1)
      for(int i=1;i<value.size();i++) assert(value[i].getLayout() == layout);
    CPSfit::write(writer,layout,"layout");
    
    std::vector<T> cen(value.size());
    for(int i=0;i<cen.size();i++) cen[i] = value[i].best();
    CPSfit::write(writer,cen,"central");
  }    
  static inline void read(HDF5reader &reader, std::vector<superJackknifeDistribution<T> > &value){
    static std::vector<std::unique_ptr<superJackknifeLayout> > layouts; //all will be deleted at the end
    layouts.push_back(std::unique_ptr<superJackknifeLayout>(new superJackknifeLayout));
    superJackknifeLayout* layout = layouts.back().get();
    CPSfit::read(reader,*layout,"layout");
    for(int i=0;i<value.size();i++) value[i].setLayout(*layout);
        
    std::vector<T> cen;
    CPSfit::read(reader,cen,"central");
    for(int i=0;i<cen.size();i++) value[i].best() = cen[i];
  }
  static inline void write(HDF5writer &writer, const std::vector<std::vector<superJackknifeDistribution<T> > > &value){
    assert(value.size() > 0 && value[0].size()>0);
    const superJackknifeLayout &layout = value[0][0].getLayout();
    for(int i=0;i<value.size();i++)
      for(int j= (i==0 ? 1:0); j<value[i].size(); j++)
	assert(value[i][j].getLayout() == layout);
    CPSfit::write(writer,layout,"layout");    
    
    int sz = 0; for(int i=0;i<value.size();i++) sz += value[i].size();
    
    std::vector<T> cen(sz);
    int off=0;
    for(int i=0;i<value.size();i++)
      for(int j=0;j<value[i].size();j++)
	cen[off++] = value[i][j].best();

    CPSfit::write(writer,cen,"central");
  }
  static inline void read(HDF5reader &reader, std::vector<std::vector<superJackknifeDistribution<T> > > &value){
    static std::vector<std::unique_ptr<superJackknifeLayout> > layouts; //all will be deleted at the end
    layouts.push_back(std::unique_ptr<superJackknifeLayout>(new superJackknifeLayout));
    superJackknifeLayout* layout = layouts.back().get();
    CPSfit::read(reader,*layout,"layout");
    for(int i=0;i<value.size();i++)
      for(int j=0;j<value[i].size();j++)
	value[i][j].setLayout(*layout);
    
    std::vector<T> cen;
    CPSfit::read(reader,cen,"central");

    int off=0;
    for(int i=0;i<value.size();i++) //has already been resized
      for(int j=0;j<value[i].size();j++)
	value[i][j].best() = cen[off++];
  } 
};


CPSFIT_END_NAMESPACE
#endif

#endif
