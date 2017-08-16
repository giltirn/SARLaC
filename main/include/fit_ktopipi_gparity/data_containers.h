#ifndef _FIT_KTOPIPI_DATA_CONTAINERS_H_
#define _FIT_KTOPIPI_DATA_CONTAINERS_H_

struct contractionData{
  std::vector<rawDataDistributionD> d;

  contractionData(const int nsample): d(8,rawDataDistributionD(nsample,0.)){
    mapping(0,'V',0,'A'); //touch mapping
  }

  void parse(std::istream &in, const int sample){
    double im; //discard imaginary part
    for(int i=0;i<8;i++)
      if(!(in >> d[i].sample(sample) >> im)) error_exit(std::cout << "type1234Data::contractions::contractionData::parse failed to parse contraction data for config " << sample << "\n");
  }
    
  inline int map_off(const int i, const char A, const int j, const char B) const{
    int iA = A == 'V' ? 0 : 1;
    int iB = B == 'V' ? 0 : 1;
    return i + 2*(iA + 2*(j + 2*iB));      
  }
      
  int mapping(const int i, const char A, const int j, const char B) const{
    static int map[16];
    static bool initted = false;
    if(!initted){
      map[ map_off(0,'V',0,'A') ] = 0;
      map[ map_off(0,'A',0,'V') ] = 1;
      map[ map_off(0,'V',1,'A') ] = 2;
      map[ map_off(0,'A',1,'V') ] = 3;
      map[ map_off(1,'V',0,'A') ] = 4;
      map[ map_off(1,'A',0,'V') ] = 5;
      map[ map_off(1,'V',1,'A') ] = 6;
      map[ map_off(1,'A',1,'V') ] = 7;
      initted = true;
    }
    return map[map_off(i,A,j,B)];
  }
    
  inline rawDataDistributionD & operator()(const int i, const char A, const int j, const char B){
    return d[mapping(i,A,j,B)];
  }
  inline const rawDataDistributionD & operator()(const int i, const char A, const int j, const char B) const{
    return d[mapping(i,A,j,B)];
  }
  inline contractionData operator+(const contractionData &r) const{
    contractionData out(*this);
    for(int i=0;i<8;i++) out.d[i] = out.d[i] + r.d[i];
    return out;
  }
  inline contractionData operator*(const double d) const{
    contractionData out(*this);
    for(int i=0;i<8;i++) out.d[i] = out.d[i] * d;
    return out;
  }

  bool isZero() const{
    for(int i=0;i<8;i++)
      for(int s=0;s<d[i].size();s++)
	if(d[i].sample(s) != 0.) return false;
    return true;
  }
    
};
    
struct contractions{
  std::vector<contractionData> con;
  std::vector<rawDataDistributionD> extra; //F0 g5, -F1 g5

  contractions(){}
  contractions(const int nsample, const int ncontract, const int nextra = 0): con(ncontract, contractionData(nsample)), extra(nextra,rawDataDistributionD(nsample,0.)){}

  //Parse the ncontract * 8 * 2 double entries for a given tK, t, sample
  void parse(std::istream &in, const int sample){
    for(int c=0;c<con.size();c++) con[c].parse(in,sample);

    double im; //discard imaginary part
    for(int e=0;e<extra.size();e++)
      if(!(in >> extra[e].sample(sample) >> im) ) error_exit(std::cout << "type1234Data::contractions::parse failed to parse extra idx " << e << " data for config " << sample << "\n");    
  }

  const contractionData & C(const int i) const{ assert(i>=0 && i<con.size()); return con[i]; }
  const rawDataDistributionD & mix() const{ return extra[0]; } //only use the F0 g5 spin-flavor structure in practise 
  
  inline bool isZero() const{ 
    for(int i=0;i<con.size();i++) if(!con[i].isZero()) return false;

    for(int e=0;e<extra.size();e++)
      for(int s=0;s<extra[e].size();s++)
	if(extra[e].sample(s) != 0.) return false;
    return true;    
  }

  contractions operator+(const contractions &r){
    if(r.con.size() != con.size() || r.extra.size() != extra.size()) error_exit(std::cout << "contractions::operator+ cannot sum two contractions of different sizes!\n");
    contractions out(*this);
    for(int i=0;i<con.size();i++) out.con[i] = out.con[i] + r.con[i];
    for(int i=0;i<extra.size();i++) out.extra[i] = out.extra[i] + r.extra[i];
    return out;
  }
  contractions operator*(const double d){
    contractions out(*this);
    for(int i=0;i<con.size();i++) out.con[i] = out.con[i] * d;
    for(int i=0;i<extra.size();i++) out.extra[i] = out.extra[i] * d;
    return out;
  }      
};



class type1234Data{
  typedef contractions tdata; //data associated with given t, i.e. the contractions
  typedef std::vector<tdata> tKdata;  //data associated with a given tK, i.e. Lt * tdata
  
  std::vector<tKdata> data; //[tK][t]

  int type;
  int Lt;
  int nsample;
  int ncontract;
  int nextra;
public:
  type1234Data(const int _type, const int _Lt, const int _nsample): type(_type),Lt(_Lt), nsample(_nsample), nextra(0){
    switch(type){
    case 1:
    case 2:
      ncontract = 6; break;
    case 3:
    case 4:
      ncontract = 10;
      nextra = 2; break;  
    default:
      error_exit(std::cout << "type1234Data::type1234Data invalid type " << type << std::endl);
    };

    data.resize(Lt, tKdata(Lt, tdata(nsample,ncontract,nextra)));
  }
  inline int getType() const{ return type; }
  inline int getLt() const{ return Lt; }
  inline int getNsample() const{ return nsample; }
  inline int getNcontract() const{ return ncontract; }
  inline int getNextra() const{ return nextra; }

  const contractions & operator()(const int tK, const int t) const{
    return data[tK][t];
  }

  std::vector<int> getNonZeroKaonTimeslices() const{
    std::vector<int> out;
    for(int tK=0;tK<Lt;tK++){
      bool allzero = true;
      for(int t=0;t<Lt;t++) if(!data[tK][t].isZero()){ allzero=false; break; }
      if(!allzero) out.push_back(tK);
    }
    return out;
  }  
  void parse(std::istream &in, const int sample){
    int tK, t;
    for(int tK_expect=0;tK_expect<getLt();tK_expect++){      
      for(int t_expect=0;t_expect<getLt();t_expect++){
	if(!(in >> tK)) error_exit(std::cout << "type1234Data::parse failed to read tK for config " << sample << "\n");
	if(tK != tK_expect) error_exit(std::cout << "type1234Data::parse tK doesn't match expectations: " << tK << ":" << tK_expect << " for config " << sample << "\n");

	if(!(in >> t)) error_exit(std::cout << "type1234Data::parse failed to read t for config " << sample << "\n");
	if(t != t_expect) error_exit(std::cout << "type1234Data::parse t doesn't match expectations: " << t << ":" << t_expect << " for config " << sample << "\n");

	data[tK][t].parse(in,sample);
      }
    }
  }
  void parse(const std::string &filename, const int sample){
    std::cout << "type1234Data parsing " << filename << " to sample " << sample << std::endl;
    std::ifstream is(filename.c_str());
    if(!is.good()){ std::cout << "Could not open file \"" << filename << "\"\n"; std::cout.flush(); exit(-1); }
    parse(is,sample);
    if(is.fail() || is.bad()){ std::cout << "Error reading file \"" << filename << "\"\n"; std::cout.flush(); exit(-1); }
    is.close();
  }  
  inline type1234Data operator+(const type1234Data &r) const{
    type1234Data out(*this);
#pragma omp parallel for
    for(int tt=0;tt<Lt*Lt;tt++){
      int tK = tt / Lt,  t = tt % Lt;
      out.data[tK][t] = out.data[tK][t] + r.data[tK][t];
    }
    return out;
  }
  inline type1234Data operator*(const double r) const{
    type1234Data out(*this);
#pragma omp parallel for
    for(int tt=0;tt<Lt*Lt;tt++){
      int tK = tt / Lt,  t = tt % Lt;
      out.data[tK][t] = out.data[tK][t] * r;
    }
    return out;
  }  
    
};


#endif
