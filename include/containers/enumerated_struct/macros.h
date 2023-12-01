#ifndef _ENUMERATED_STRUCT_MACROS_H_
#define _ENUMERATED_STRUCT_MACROS_H_

//Because who doesn't love preprocessor macros?

#include<iostream>
#include<cassert>
#include<boost/preprocessor/punctuation/comma_if.hpp>
#include<boost/preprocessor/seq/elem.hpp>
#include<boost/preprocessor/seq/size.hpp>
#include<boost/preprocessor/seq/for_each.hpp>
#include<boost/preprocessor/seq/for_each_i.hpp>
#include<boost/preprocessor/tuple/elem.hpp>
#include<boost/preprocessor/cat.hpp>
#include<boost/preprocessor/stringize.hpp>
#include<boost/preprocessor/control/if.hpp>

#include<config.h>
#include<utils/utils.h>
#include<ET/generic_ET.h>
#include<parser/parser.h>
#include<serialize/hdf5_serialize.h>

SARLAC_START_NAMESPACE

//Extract info from DEF
#define _ENUMERATED_STRUCT_DEF_STRUCT(DEF)  BOOST_PP_TUPLE_ELEM(0,DEF)
#define _ENUMERATED_STRUCT_DEF_TYPE(DEF) BOOST_PP_TUPLE_ELEM(1,DEF)
#define _ENUMERATED_STRUCT_DEF_SEQ(DEF) BOOST_PP_TUPLE_ELEM(2,DEF)
#define _ENUMERATED_STRUCT_DEF_DEFAULT(DEF) BOOST_PP_TUPLE_ELEM(3,DEF)

//Define the members
#define _ENUMERATED_STRUCT_DEF_MEMBER(R, TYPE, ELEM) TYPE ELEM;

#define _ENUMERATED_STRUCT_DEF_MEMBERS(DEF)			\
  BOOST_PP_SEQ_FOR_EACH(_ENUMERATED_STRUCT_DEF_MEMBER, _ENUMERATED_STRUCT_DEF_TYPE(DEF), _ENUMERATED_STRUCT_DEF_SEQ(DEF))

//Constructor (default, value)
#define _ENUMERATED_STRUCT_DEF_CONSTRUCT_ELEM_IN(R,DEF,IDX,ELEM) BOOST_PP_COMMA_IF(IDX) _ENUMERATED_STRUCT_DEF_TYPE(DEF) BOOST_PP_CAT(_,ELEM) = BOOST_PP_SEQ_ELEM(IDX,_ENUMERATED_STRUCT_DEF_DEFAULT(DEF) )
#define _ENUMERATED_STRUCT_DEF_CONSTRUCT_ELEM_SET(R,DUMMY,IDX,ELEM) BOOST_PP_COMMA_IF(IDX) ELEM(BOOST_PP_CAT(_,ELEM))

#define _ENUMERATED_STRUCT_DEF_CONSTRUCTOR(DEF) _ENUMERATED_STRUCT_DEF_STRUCT(DEF)(BOOST_PP_SEQ_FOR_EACH_I(_ENUMERATED_STRUCT_DEF_CONSTRUCT_ELEM_IN, DEF, _ENUMERATED_STRUCT_DEF_SEQ(DEF)) ): BOOST_PP_SEQ_FOR_EACH_I(_ENUMERATED_STRUCT_DEF_CONSTRUCT_ELEM_SET, , _ENUMERATED_STRUCT_DEF_SEQ(DEF)) {}

//Accessors
#define _ENUMERATED_STRUCT_DEF_MEM_ACCESSOR(R,DUMMY,IDX,ELEM) case IDX: return ELEM;

#define _ENUMERATED_STRUCT_DEF_ACCESSOR(DEF)			\
  _ENUMERATED_STRUCT_DEF_TYPE(DEF) & operator()(const int i){	\
  switch(i){					\
  BOOST_PP_SEQ_FOR_EACH_I(_ENUMERATED_STRUCT_DEF_MEM_ACCESSOR, , _ENUMERATED_STRUCT_DEF_SEQ(DEF))	\
  default:								\
    SARLaC::error_exit(std::cout << BOOST_PP_STRINGIZE(_ENUMERATED_STRUCT_DEF_STRUCT(DEF)) << "::operator() invalid index " << i <<std::endl); \
  }									\
  }									\
  inline const _ENUMERATED_STRUCT_DEF_TYPE(DEF) & operator()(const int i) const{ return const_cast<_ENUMERATED_STRUCT_DEF_STRUCT(DEF)*>(this)->operator()(i); } \
  inline _ENUMERATED_STRUCT_DEF_TYPE(DEF) & operator[](const int i){ return this->operator()(i); }	\
  inline const _ENUMERATED_STRUCT_DEF_TYPE(DEF) & operator[](const int i) const{ return this->operator()(i); }

//Size method
#define _ENUMERATED_STRUCT_DEF_SIZEFUNC(DEF)				\
  inline int size() const{ return BOOST_PP_SEQ_SIZE(_ENUMERATED_STRUCT_DEF_SEQ(DEF)); } \
  inline void resize(const int sz){ assert(sz == this->size()); }

//Zero method
#define _ENUMERATED_STRUCT_DEF_ZERO_MEM(R,DUMMY,ELEM) SARLaC::zeroit(ELEM);

#define _ENUMERATED_STRUCT_DEF_ZEROFUNC(DEF)			\
  inline void zero(){				\
  BOOST_PP_SEQ_FOR_EACH(_ENUMERATED_STRUCT_DEF_ZERO_MEM, , _ENUMERATED_STRUCT_DEF_SEQ(DEF))\
    }

//Print method
#define _ENUMERATED_STRUCT_DEF_PRINT_MEM(R,DUMMY,ELEM) << BOOST_PP_STRINGIZE(ELEM) << "=" << ELEM << " "

#define _ENUMERATED_STRUCT_DEF_PRINTFUNC(DEF) \
  inline std::string print() const{ std::ostringstream os; os BOOST_PP_SEQ_FOR_EACH(_ENUMERATED_STRUCT_DEF_PRINT_MEM, , _ENUMERATED_STRUCT_DEF_SEQ(DEF)); return os.str(); }

//Equivalence method
#define _ENUMERATED_STRUCT_DEF_EQUIV_ELEM(R,DUMMY,IDX,ELEM) BOOST_PP_IF(IDX, &&, ) r.ELEM == ELEM

#define _ENUMERATED_STRUCT_DEF_EQUIV(DEF) \
  inline bool operator==(const _ENUMERATED_STRUCT_DEF_STRUCT(DEF) &r) const{ return BOOST_PP_SEQ_FOR_EACH_I(_ENUMERATED_STRUCT_DEF_EQUIV_ELEM, , _ENUMERATED_STRUCT_DEF_SEQ(DEF)); }

//Get member name by index as string
#define _ENUMERATED_STRUCT_DEF_MEM_STRING_VAL(R,DUMMY,IDX,ELEM) case IDX: return BOOST_PP_STRINGIZE(ELEM);

#define _ENUMERATED_STRUCT_DEF_MEM_STRING(DEF)	\
  std::string memberName(const int i) const{	\
    switch(i){								\
      BOOST_PP_SEQ_FOR_EACH_I(_ENUMERATED_STRUCT_DEF_MEM_STRING_VAL, , _ENUMERATED_STRUCT_DEF_SEQ(DEF)) \
    default:							\
      SARLaC::error_exit(std::cout << BOOST_PP_STRINGIZE(_ENUMERATED_STRUCT_DEF_STRUCT(DEF)) << "::memberName(int) invalid index " << i <<std::endl); \
    };									\
    }


//Put the body together....
#define _ENUMERATED_STRUCT_DEF_STRUCT_BODY(DEF)			\
    _ENUMERATED_STRUCT_DEF_MEMBERS(DEF)				\
      _ENUMERATED_STRUCT_DEF_CONSTRUCTOR(DEF)			\
      _ENUMERATED_STRUCT_DEF_ACCESSOR(DEF)				\
      _ENUMERATED_STRUCT_DEF_SIZEFUNC(DEF)				\
      _ENUMERATED_STRUCT_DEF_ZEROFUNC(DEF)				\
      _ENUMERATED_STRUCT_DEF_PRINTFUNC(DEF)			\
    _ENUMERATED_STRUCT_DEF_EQUIV(DEF)				\
    _ENUMERATED_STRUCT_DEF_MEM_STRING(DEF)

//Interface with parser generator
#define _ENUMERATED_STRUCT_DEF_PARSER_ARG_TRANSFORM_ELEM(R,TYPE,ELEM)(TYPE,ELEM)

#define _ENUMERATED_STRUCT_DEF_PARSER_ARG_TRANSFORM(DEF) BOOST_PP_SEQ_FOR_EACH(_ENUMERATED_STRUCT_DEF_PARSER_ARG_TRANSFORM_ELEM, _ENUMERATED_STRUCT_DEF_TYPE(DEF), _ENUMERATED_STRUCT_DEF_SEQ(DEF))

//Final calls
#define _ENUMERATED_STRUCT_DEF_ENUMERATED_STRUCT(DEF)		\
  struct _ENUMERATED_STRUCT_DEF_STRUCT(DEF){			\
    _ENUMERATED_STRUCT_DEF_STRUCT_BODY(DEF);				\
    ENABLE_GENERIC_ET(_ENUMERATED_STRUCT_DEF_STRUCT(DEF),_ENUMERATED_STRUCT_DEF_STRUCT(DEF),_ENUMERATED_STRUCT_DEF_STRUCT(DEF)); \
    GENERATE_HDF5_SERIALIZE_METHOD(_ENUMERATED_STRUCT_DEF_SEQ(DEF));	\
  };									\
  GENERATE_PARSER(_ENUMERATED_STRUCT_DEF_STRUCT(DEF), _ENUMERATED_STRUCT_DEF_PARSER_ARG_TRANSFORM(DEF)); \
  GENERATE_HDF5_SERIALIZE_FUNC(_ENUMERATED_STRUCT_DEF_STRUCT(DEF));

#define _ENUMERATED_STRUCT_DEF_ENUMERATED_STRUCT_MEMBER(SCOPE,DEF)	\
  struct _ENUMERATED_STRUCT_DEF_STRUCT(DEF){			\
    _ENUMERATED_STRUCT_DEF_STRUCT_BODY(DEF);				\
    ENABLE_GENERIC_ET(_ENUMERATED_STRUCT_DEF_STRUCT(DEF),_ENUMERATED_STRUCT_DEF_STRUCT(DEF),_ENUMERATED_STRUCT_DEF_STRUCT(DEF)); \
    GENERATE_HDF5_SERIALIZE_METHOD(_ENUMERATED_STRUCT_DEF_SEQ(DEF));	\
  };	

#define _ENUMERATED_STRUCT_DEF_ENUMERATED_STRUCT_MEMBER_EXTERNAL(SCOPE,DEF) \
  GENERATE_PARSER_GM(SCOPE::_ENUMERATED_STRUCT_DEF_STRUCT(DEF), BOOST_PP_SEQ_CAT( (SCOPE)(_)(_ENUMERATED_STRUCT_DEF_STRUCT(DEF))(_grammar) ),_ENUMERATED_STRUCT_DEF_PARSER_ARG_TRANSFORM(DEF)); \
  GENERATE_HDF5_SERIALIZE_FUNC(SCOPE::_ENUMERATED_STRUCT_DEF_STRUCT(DEF));


SARLAC_END_NAMESPACE

#endif
