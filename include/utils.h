#ifndef _UTILS_H__
#define _UTILS_H__

class OstreamHook{
public:
  virtual void write(std::ostream &) const = 0;
};

inline std::ostream & operator<<(std::ostream &os, const OstreamHook &hk){
  hk.write(os);
  return os;
}


#endif
