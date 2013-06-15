#include <boost/python.hpp>
#include <boost/numpy.hpp>

#include "../../fitfun.h"

namespace p = boost::python;
namespace np = boost::numpy;

double doubleit(double x) {
    return 2.*x;
}

void fill_array(np::ndarray &y) {
    int i;
    double *d = (double*) y.get_data();

    for (i=0; i<y.shape(0); i++) {
        d[i] = i;
    }
}

BOOST_PYTHON_MODULE(fitfun_boost_numpy) {
    p::def("doubleit", doubleit);
    
    np::initialize();
    p::def("fill_array", fill_array);
}

