#include <distribution/jackknife.h>
#include <distribution/boot_jackknife.h>
#include <utils/template_wizardry.h>
using namespace CPSfit;

int main(void){

  {
    //Test get base type
    static_assert(std::is_same< typename getBaseType<jackknifeDistribution<double> >::type, double >::value == true );
    static_assert(std::is_same< typename getBaseType<bootJackknifeDistribution<double> >::type, double >::value == true );
  }

  {
    //Test Complexify
    static_assert(std::is_same< Complexify<double>, std::complex<double> >::value == true);
    static_assert(std::is_same< Complexify<jackknifeDistribution<double> >, jackknifeDistribution<std::complex<double> > >::value == true);
    static_assert(std::is_same< Complexify<bootJackknifeDistribution<double> >, bootJackknifeDistribution<std::complex<double> > >::value == true);
  }

  {
    //Test Realify
    static_assert(std::is_same< Realify<std::complex<double> >, double >::value == true);
    static_assert(std::is_same< Realify<jackknifeDistribution<std::complex<double> > >, jackknifeDistribution<double> >::value == true);
    static_assert(std::is_same< Realify<bootJackknifeDistribution<std::complex<double> > >, bootJackknifeDistribution<double> >::value == true);
  }

  return 0;
}
