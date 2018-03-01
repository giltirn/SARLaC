#ifndef _ENUMERATED_STRUCT_H_
#define _ENUMERATED_STRUCT_H_

//Macros for automagically generating structs with operator()(const int) and operator[](const int) accessors, etc as well as expression template, hdf5 serialization, and parser hooks.
//Useful for containing fit function parameters

#include<containers/enumerated_struct/macros.h> //<---avert thine eyes!

//Generate a struct with members of a fixed type with operator()(const int) and operator[](const int) accessors, etc as well as expression template and parser hooks.
//Useful for containing fit function parameters
//Call with an argument DEF of the form  ( structname, element type, (elem name 0)(elem name 1)....,  (default val 0)(default val 1)..... )
//where the ellipsis indicates an arbitrary number of further () sequence entries
#define DEF_ENUMERATED_STRUCT(DEF) _ENUMERATED_STRUCT_DEF_ENUMERATED_STRUCT(DEF)

//Same as above but for a struct defined within the body of another enclosing class, eg
//class Enclosing{
//struct MyStruct;
//...
//}
#define DEF_ENUMERATED_STRUCT_MEMBER(ENCLOSING_CLASS,DEF) _ENUMERATED_STRUCT_DEF_ENUMERATED_STRUCT_MEMBER(ENCLOSING_CLASS,DEF)

//If a parser and/or HDF5 serialization is required for this internal struct, the following should be called with the same arguments outside the class scope
#define DEF_ENUMERATED_STRUCT_MEMBER_EXTERNAL(ENCLOSING_CLASS,DEF) _ENUMERATED_STRUCT_DEF_ENUMERATED_STRUCT_MEMBER_EXTERNAL(ENCLOSING_CLASS,DEF)

//cf examples/enumerated_struct_example.C  for an example!

#endif
