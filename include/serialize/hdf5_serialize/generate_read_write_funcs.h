#ifndef _GENERATE_READ_WRITE_FUNCS_H_
#define _GENERATE_READ_WRITE_FUNCS_H_

//Some useful macros to automatically generate the read/write functions for user classes and structs

#include<type_traits>

#include<config.h>
#include<utils/macros.h>

#ifdef HAVE_HDF5

#include<boost/preprocessor.hpp>
#include "hdf5_writer.h"

SARLAC_START_NAMESPACE

template<typename T>
struct _force_external_lookup{
  inline static void fwrite(HDF5writer &writer, const T& elem, const std::string &tag){
    write(writer,elem,tag);
  }
  inline static void fread(HDF5reader &reader, T &elem, const std::string &tag){ 
    read(reader,elem,tag);
  }
};

#define _GENERATE_HDF5_SERIALIZE_METHOD_WRITEIT(r,data,elem) ::SARLaC::_force_external_lookup<typename std::decay<decltype(this->elem)>::type>::fwrite(writer, this-> elem, BOOST_PP_STRINGIZE(elem)); 
#define _GENERATE_HDF5_SERIALIZE_METHOD_READIT(r,data,elem) ::SARLaC::_force_external_lookup<typename std::decay<decltype(this->elem)>::type>::fread(reader, this-> elem, BOOST_PP_STRINGIZE(elem)); 

//Generate HDF5 serialize methods for a class. MEMBERS should be a series of member names in parentheses, eg (member1)(member2)(member3)....
#define GENERATE_HDF5_SERIALIZE_METHOD(MEMBERS)\
  void write(::SARLaC::HDF5writer &writer, const std::string &tag) const{ \
    writer.enter(tag);							\
    BOOST_PP_SEQ_FOR_EACH(_GENERATE_HDF5_SERIALIZE_METHOD_WRITEIT, , MEMBERS); \
    writer.leave();							\
  }									\
  void read(::SARLaC::HDF5reader &reader, const std::string &tag){	\
    reader.enter(tag);							\
    BOOST_PP_SEQ_FOR_EACH(_GENERATE_HDF5_SERIALIZE_METHOD_READIT, , MEMBERS); \
    reader.leave();							\
  }

#define GENERATE_HDF5_SERIALIZE_FUNC(CLASSNAME)\
  inline void write(::SARLaC::HDF5writer &writer, const CLASSNAME &d, const std::string &tag){ d.write(writer,tag); } \
  inline void read(::SARLaC::HDF5reader &reader, CLASSNAME &d, const std::string &tag){ d.read(reader,tag); }


#define GENERATE_HDF5_ENUM_SERIALIZE(ENUMNAME) \
  inline void write(SARLaC::HDF5writer &writer, const ENUMNAME d, const std::string &tag){ write(writer, (int)d, tag); } \
  inline void read(SARLaC::HDF5reader &reader, ENUMNAME &d, const std::string &tag){ int i; read(reader, i, tag); d = (ENUMNAME)i; }


SARLAC_END_NAMESPACE

#else  //if not HAVE_HDF5

SARLAC_START_NAMESPACE

//Empty macros so we don't need to keep adding #ifdef HAVE_HDF5 around everything
#define GENERATE_HDF5_SERIALIZE_METHOD(MEMBERS)
#define GENERATE_HDF5_SERIALIZE_FUNC(CLASSNAME)
#define GENERATE_HDF5_ENUM_SERIALIZE(ENUMNAME)

SARLAC_END_NAMESPACE

#endif

#endif
