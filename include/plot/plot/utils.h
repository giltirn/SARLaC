#ifndef MPL_PLOT_UTILS_H_
#define MPL_PLOT_UTILS_H_

#include<config.h>
#include<string>
#include<vector>
#include<iterator>

CPSFIT_START_NAMESPACE

//Get 'nhues' colors that are visually distinct
std::vector<std::string> plotColorPallete(const int nhues){
  static std::vector<std::string> colors = { "r", "b", "g", "c", "m", "k", "tab:orange", "tab:gray", "saddlebrown", "gold", "yellowgreen", "pink", "aquamarine" };
  assert(nhues <= colors.size());
  return std::vector<std::string>(colors.begin(), std::next(colors.begin(),nhues) );
}

//Get 'nsymb' symbols that are visually distinct
std::vector<std::string> plotMarkers(const int nsymb){
  static std::vector<std::string> symbols = { "o", "s", "^", "D", "v", "*", "P", "X", "<", ">", "1", "2", "3", "4" };
  assert(nsymb <= symbols.size());
  return std::vector<std::string>(symbols.begin(), std::next(symbols.begin(),nsymb) );
}


CPSFIT_END_NAMESPACE


#endif
