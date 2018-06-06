#ifndef _CPSFIT_MPL_SCRIPT_GENERATE_BASE_H_
#define _CPSFIT_MPL_SCRIPT_GENERATE_BASE_H_

//Base class for classes that generate python scripts for generating plots using Matplotlib

#include<iostream>
#include<vector>

#include<config.h>
#include<utils/macros.h>
#include<plot/plot/write_python.h>

CPSFIT_START_NAMESPACE

class MatPlotLibScriptGenerateBase{
public:
  typedef std::map<std::string,KWargElem> kwargsType;
protected:
  struct kwargsPrint: public OstreamHook{
    const kwargsType &args;
    kwargsPrint(const kwargsType &_args): args(_args){}
    void write(std::ostream &os) const{
      for(kwargsType::const_iterator it = args.begin(); it != args.end(); it++){
	os << ',' << it->first << "=" << it->second;
      }
    }
  };

  std::ostringstream user; //extra python code added by user
  std::ostringstream preamble; //extra preamble python code added by user
public:

  //User code placed outside of __name__ = '__main__' scope, eg functions or global vars
  inline std::ostringstream & preinvoke(){ return preamble; }

  //Access an internal stream that can be used to provide arbitrary python commands that are included in the output file
  //IMPORTANT: User must apply a tab \t in front of each line of python code as it should be executed in the if __name__ == '__main__' scope
  inline std::ostringstream & invoke(){ return user; }  
};

CPSFIT_END_NAMESPACE
#endif
