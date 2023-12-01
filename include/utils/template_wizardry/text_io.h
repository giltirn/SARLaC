#ifndef _SARLAC_TEMPLATE_WIZARDRY_TEXT_IO_
#define _SARLAC_TEMPLATE_WIZARDRY_TEXT_IO_

//Metaprogramming constructs for obtaining information about types that can be read or written in ascii
#include<cstdlib>
#include<iostream>

#include<config.h>
#include<utils/macros.h>
#include<utils/template_wizardry/types.h>

SARLAC_START_NAMESPACE

//Check if a class T looks like an ostream
template<class T, class Fallback = void>
struct isStreamType{ enum{ value = 0 }; };

template<class T>
struct isStreamType<T, typename Void<decltype( ((T*)(NULL))->operator<<(std::endl) )  >::type >{ enum{ value = 1 }; };


//Check if a class T has a method 'parse'
template<typename T>
struct hasParseMember { 
    struct Fallback { int x; };
    struct Derived : T, Fallback { };
    template<typename C, C> struct ChT; 

    template<typename C> static char (&f(ChT<int Fallback::*, &C::parse>*))[1]; 
    template<typename C> static char (&f(...))[2]; 

    static bool const value = sizeof(f<Derived>(0)) == 2;
}; 

SARLAC_END_NAMESPACE
#endif
