#ifndef _CPSFIT_UTILS_MEMORY_H_
#define _CPSFIT_UTILS_MEMORY_H_

#include<iostream>
#include<fstream>
#include<limits>
#include<sys/sysinfo.h>

#include<config.h>
#include<utils/macros.h>


CPSFIT_START_NAMESPACE

//Print memory usage
inline void printMem(const std::string &reason = "", FILE* stream = stdout){
  if(reason != "") fprintf(stream, "printMem called with reason '%s': ", reason.c_str());
  else fprintf(stream, "printMem: ");
  
  struct sysinfo myinfo;
  sysinfo(&myinfo);  //cf http://man7.org/linux/man-pages/man2/sysinfo.2.html
  double total_mem = myinfo.mem_unit * myinfo.totalram;
  total_mem /= (1024.*1024.);
  double free_mem = myinfo.mem_unit * myinfo.freeram;
  free_mem /= (1024.*1024.);

  fprintf(stream,"Memory: total: %.2f MB, avail: %.2f MB, used %.2f MB\n",total_mem, free_mem, total_mem-free_mem);

#define REPORT_MEMINFO
#ifdef REPORT_MEMINFO
  {
    double total_mem2;
    double free_mem2;
  
    std::string token;
    std::ifstream file("/proc/meminfo");
    while(file >> token) {
        if(token == "MemTotal:") {
	  file >> total_mem2;
	}else if(token == "MemAvailable:"){
	  file >> free_mem2;
	}
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');        // ignore rest of the line
    }
    file.close();

    fprintf(stream,"Memory (/proc/meminfo): total: %.2f MB, avail: %.2f MB, used %.2f MB\n",total_mem2/1024., free_mem2/1024., (total_mem2-free_mem2)/1024.);
  }
#endif
  
  fflush(stream);
}

CPSFIT_END_NAMESPACE
#endif
