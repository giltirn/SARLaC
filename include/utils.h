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

template<typename T>
inline T & real_ref(std::complex<T> &z){
  return reinterpret_cast<T(&)[2]>(z)[0];
}
template<typename T>
inline T & imag_ref(std::complex<T> &z){
  return reinterpret_cast<T(&)[2]>(z)[1];
}

//Substitute substring '%d' with configuration idx
inline std::string subsIdx(const std::string fmt, const int idx){
  std::string::size_type off = fmt.find("%d");
  if(off == std::string::npos){
    std::cout << "Could not find substring \"%d\" in format string " << fmt << std::endl;
    std::cout.flush();
    exit(-1);
  }
  std::ostringstream os; os << idx;
  std::string out(fmt);
  out.replace(off,2,os.str());
  return out;
}


#endif
