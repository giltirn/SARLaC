#ifndef _GLOBAL_DATA_MAP
#define _GLOBAL_DATA_MAP

#include<regex>
#include<set>
#include<fstream>

#include<config.h>
#include<utils/macros.h>
#include<utils/utils/string.h>

CPSFIT_START_NAMESPACE

//Expect a file with line of the following format:
//<OUTER CONFIG> <DIRECTORY> <INNER_CONFIG> <TAG>
//The directory and inner config along with an arbitrary filename format specify a data file. Outer config is an index that is used to specify data files that belong to the same
//configuration (they should be on the same gauge configuration!). There can be multiple tags for different data types or groups of data that can potentially lie in different directories.
struct DataLocationInfo{
  int outer_config;
  std::string directory;
  int inner_config;
  std::string tag;

  void parse(const std::string &line){
    std::regex r(R"((\d+)\s([^\s]+)\s(\d+)\s([^\s]+))");
    std::smatch m;
    if(std::regex_match(line,m,r)){
      outer_config = strToAny<int>(m[1]);
      directory = m[2];
      inner_config = strToAny<int>(m[3]);
      tag = m[4];
    }else error_exit(std::cout << "DataLocationInfo::parse could not match line \"" << line << "\"\n");
  }
  DataLocationInfo(){}
  DataLocationInfo(const std::string &line){
    parse(line);
  }
};


std::ostream & operator<<(std::ostream &os, const DataLocationInfo &info){
  os << info.outer_config << " " << info.directory << " " << info.inner_config << " " << info.tag;
  return os;
}

class DataMap{
  std::vector< std::vector<DataLocationInfo> > loc_info; //[outer_idx][meas]
  std::map<std::string, std::vector< std::pair<int,int> > > tag_map; //map  tag -> i,j 

  void generateTagMap(){
    for(int i=0;i<loc_info.size();i++){
      std::set<std::string> tags_found_i;
      for(int j=0;j<loc_info[i].size();j++){
	if(!tags_found_i.count(loc_info[i][j].tag)){
	  tag_map[loc_info[i][j].tag].push_back( std::pair<int,int>(i,j) );
	  tags_found_i.insert(loc_info[i][j].tag);
	}
      }
    }
  }
public:
  void parse(const std::string &filename){
    std::ifstream in(filename);
    assert(in.good());
    std::string line;
    
    int prev_outer = -1;
    while(std::getline(in,line)){
      DataLocationInfo info(line);
      if(info.outer_config == prev_outer){
	loc_info.back().push_back(info);
      }else{
	loc_info.push_back( std::vector<DataLocationInfo>(1,info) );
      }
      prev_outer = info.outer_config;
      assert(!in.fail() && !in.bad());
    }

    for(int i=0;i<loc_info.size();i++){
      int outer = loc_info[i].front().outer_config;
      for(int j=1;j<loc_info[i].size();j++) assert(loc_info[i][j].outer_config == outer);
    }

    generateTagMap();
  }

  std::vector<std::string> getTags() const{    
    std::vector<std::string> tags;
    for(auto it=tag_map.begin();it!=tag_map.end();++it) tags.push_back(it->first);
    return tags;
  }

  inline bool hasTag(const std::string &tag) const{ return tag_map.find(tag) != tag_map.end(); }

  const std::vector< std::pair<int,int> > & getLocationsWithTag(const std::string &tag) const{
    auto it = tag_map.find(tag);
    if(it == tag_map.end()) error_exit(std::cout << "DataMap::getLocationsWithTag could not find tag " << tag << std::endl);
    return it->second;
  }

  void getLocationsWithTag(std::set< std::pair<int,int> > &out, const std::string &tag) const{
    if(!hasTag(tag)) return;
    const std::vector< std::pair<int,int> > &data = getLocationsWithTag(tag);
    for(auto it=data.begin(); it!= data.end(); ++it) out.insert(*it);
  }

  std::set<int> getOuterConfigsWithTag(const std::string &tag) const{
    std::set<int> out;
    if(!hasTag(tag)) return out; //empty set
    const std::vector< std::pair<int,int> > &data = getLocationsWithTag(tag);
    for(auto it=data.begin(); it!= data.end(); ++it) out.insert(it->first);
    return out;
  }

  DataMap(){}
  DataMap(const std::string &filename){ parse(filename); }
  
  size_t size() const{ return loc_info.size(); }
  const std::vector<DataLocationInfo> & operator[](const size_t i) const{ return loc_info[i]; } 
  
  const DataLocationInfo & operator()(const size_t i, const std::string &tag) const{
    for(int j=0;j<loc_info[i].size();j++) if(loc_info[i][j].tag == tag) return loc_info[i][j];
    error_exit(std::cout << "DataMap::operator() could not find tag \"" << tag << "\" for outer config " << i << std::endl);
  }

  DataLocationInfo const* getInfoPtr(const size_t i, const std::string &tag) const{
    for(int j=0;j<loc_info[i].size();j++) if(loc_info[i][j].tag == tag) return &loc_info[i][j];
    return NULL;
  }

  const DataLocationInfo & operator()(const std::pair<int,int> &cp) const{ return loc_info[cp.first][cp.second]; }

  //For each outer configuration oconf in 'outer_confs' produce a mapping of oconf to the first DataLocationInfo found for oconf matching a tag in the 'tags' list
  std::map<int, DataLocationInfo const*> getInfoPtrMapping(const std::set<int> &outer_confs, const std::vector<std::string> &tags) const{
    std::map<int, DataLocationInfo const*> out;
    
    DataLocationInfo const* data_loc;
    for(auto oconf=outer_confs.begin(); oconf != outer_confs.end(); ++oconf){
      bool found = false;
      for(auto tag = tags.begin(); tag != tags.end(); tag++){ //search over tags for this data type to find the data for this outer config
	if( (data_loc = this->getInfoPtr(*oconf, *tag) ) != NULL ){
	  out[*oconf] = data_loc;
	  found=true; break;
	}
      }
      if(!found){
	std::cout << "DataMap::getInfoPtrMapping Could not find data with outer config " << *oconf << " with any of the tags: " << tags << "\n";
	std::cout << "Available tags on this config are:"; 
	for(auto t=(*this)[*oconf].begin(); t!=(*this)[*oconf].end();t++) std::cout << " " << t->tag;
	error_exit(std::cout << std::endl);
      }    
    }
    return out;
  }

};

//For set of outer configs, return the (directory, inner_config) pairs for a given tag
std::vector<std::pair<std::string, int> > getDirectoryConfigList(const std::set<int> &outer_configs, const std::string &tag, const DataMap &dmap){
  std::vector<std::pair<std::string, int> > toread;
  for(auto it=outer_configs.begin(); it != outer_configs.end(); ++it){
    const DataLocationInfo &info = dmap(*it, tag);
    toread.push_back( std::pair<std::string, int>(info.directory, info.inner_config) );
  }
  return toread;
}

//Prune set of outer samples such that they align and fill bins of size 'bin_size'
std::set<int> binningPrune(const std::set<int> &in, const int bin_size){
  int osample_min = INT_MAX;
  int osample_max = INT_MIN;
  for(auto it= in.begin(); it != in.end(); it++){
    osample_min = std::min(osample_min, *it);
    osample_max = std::max(osample_max, *it);
  }  

  std::cout << "binningPrune binning by " << bin_size << " set of size " << in.size() << " osample range " << osample_min << ":" << osample_max << std::endl;
  
  std::set<int> out(in);

  //Align start and end of binning region
  int bin_start = int( floor(float(osample_min)/bin_size) ) * bin_size;
  int bin_end = int( ceil(float(osample_max)/bin_size) ) * bin_size;

  std::cout << "bin range " << bin_start << ":" << bin_end << std::endl;

  for(int b = bin_start; b <= bin_end; b += bin_size){
    bool full_bin = true;
    for(int o=b; o < b+bin_size; o++)
      if(!out.count(o)){
	std::cout << "bin " << b << ":" << b+bin_size-1 << " missing " << o << "\n";
	full_bin = false; //break;
      }
    if(!full_bin){
      std::cout << "bin " << b << ":" << b+bin_size-1 << " incomplete\n";
      for(int o=b; o < b+bin_size; o++)
	out.erase(o);
    }
  }
  std::cout << "binningPrune orig size " << in.size() << " output size " << out.size() << std::endl;

  return out;
}

CPSFIT_END_NAMESPACE

#endif
