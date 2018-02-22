#ifndef _HDF5_READWRITE_OTHER_H___
#define _HDF5_READWRITE_OTHER_H___

//Functions to read and write some other common container types
#include<config.h>

#ifdef HAVE_HDF5
#include<map>
#include<serialize/hdf5_serialize/read_write_basic.h>

CPSFIT_START_NAMESPACE

//Pair
template<typename T, typename U>
inline static void write(HDF5writer &writer, const std::pair<T,U> &value, const std::string &tag){
  writer.enter(tag); //enter a group
  write(writer, value.first, "first");
  write(writer, value.second, "second");
  writer.leave();
}
template<typename T, typename U>
inline static void read(HDF5reader &reader, std::pair<T,U> &value, const std::string &tag){
  reader.enter(tag); //enter a group
  read(reader, value.first, "first");
  read(reader, value.second, "second");
  reader.leave();
}

//Map
template<typename KeyT, typename DataT>
static void write(HDF5writer &writer, const std::map<KeyT,DataT> &value, const std::string &tag){
  writer.enter(tag);
  int sz = value.size();
  write(writer,sz,"sz");
  writer.enter("entries");
  int i=0;
  for(typename std::map<KeyT,DataT>::const_iterator it = value.begin(); it != value.end(); it++, i++){
    const std::pair<KeyT,DataT> &e = *it;
    std::ostringstream os; os << "entry_" << i;
    write(writer,e,os.str());
  }
  writer.leave();
  writer.leave();
}
template<typename KeyT, typename DataT>
static void read(HDF5reader &reader, std::map<KeyT,DataT> &value, const std::string &tag){
  value.clear();
  reader.enter(tag);
  int sz;
  read(reader,sz,"sz");
  reader.enter("entries");
  for(int i=0;i<sz;i++){
    std::pair<KeyT,DataT> e;
    std::ostringstream os; os << "entry_" << i;
    read(reader,e,os.str());
    value.emplace(std::move(e));
  }
  assert(value.size() == sz);
  reader.leave();
  reader.leave();
}

CPSFIT_END_NAMESPACE

#endif

#endif
