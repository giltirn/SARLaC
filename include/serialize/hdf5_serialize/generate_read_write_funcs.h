#ifndef _GENERATE_READ_WRITE_FUNCS_H_
#define _GENERATE_READ_WRITE_FUNCS_H_

//Some useful macros to automatically generate the read/write functions for user classes and structs

#include<type_traits>

#include<config.h>
#include<utils/macros.h>

#ifdef HAVE_HDF5

#include<boost/preprocessor.hpp>

CPSFIT_START_NAMESPACE

template<typename T>
struct _force_external_lookup{
  inline static void fwrite(HDF5writer &writer, const T& elem, const std::string &tag){
    write(writer,elem,tag);
  }
  inline static void fread(HDF5reader &reader, T &elem, const std::string &tag){ 
    read(reader,elem,tag);
  }
};

#define _GENERATE_HDF5_SERIALIZE_METHOD_WRITEIT(r,data,elem) _force_external_lookup<typename std::decay<decltype(this->elem)>::type>::fwrite(writer, this-> elem, BOOST_PP_STRINGIZE(elem)); 
#define _GENERATE_HDF5_SERIALIZE_METHOD_READIT(r,data,elem) _force_external_lookup<typename std::decay<decltype(this->elem)>::type>::fread(reader, this-> elem, BOOST_PP_STRINGIZE(elem)); 

//Generate HDF5 serialize methods for a class. MEMBERS should be a series of member names in parentheses, eg (member1)(member2)(member3)....
#define GENERATE_HDF5_SERIALIZE_METHOD(MEMBERS)\
    void write(HDF5writer &writer, const std::string &tag) const{ \
      writer.enter(tag);\
      BOOST_PP_SEQ_FOR_EACH(_GENERATE_HDF5_SERIALIZE_METHOD_WRITEIT, , MEMBERS); \
      writer.leave(); \
    }\
    void read(HDF5reader &reader, const std::string &tag){ \
      reader.enter(tag); \
      BOOST_PP_SEQ_FOR_EACH(_GENERATE_HDF5_SERIALIZE_METHOD_READIT, , MEMBERS); \
      reader.leave(); \
    }

#define GENERATE_HDF5_SERIALIZE_FUNC(CLASSNAME)\
  inline void write(HDF5writer &writer, const CLASSNAME &d, const std::string &tag){ d.write(writer,tag); }\
  inline void read(HDF5reader &reader, CLASSNAME &d, const std::string &tag){ d.read(reader,tag); }

CPSFIT_END_NAMESPACE

#else  //if not HAVE_HDF5

CPSFIT_START_NAMESPACE

//Empty macros so we don't need to keep adding #ifdef HAVE_HDF5 around everything
#define GENERATE_HDF5_SERIALIZE_METHOD(MEMBERS)
#define GENERATE_HDF5_SERIALIZE_FUNC(CLASSNAME)

CPSFIT_END_NAMESPACE

#endif

#endif
