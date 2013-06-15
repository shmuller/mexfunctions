#include <boost/python.hpp>
#include <numpy/arrayobject.h>

#include "../../fitfun.h"

using namespace boost::python;

double doubleit(double x) {
    return 2.*x;
}

void fill_array(numeric::array &y) {
    int i;
    double *d = (double*) PyArray_DATA(y.ptr());

    for (i=0; i<3; i++) {
        d[i] = i;
    }
}

BOOST_PYTHON_MODULE(fitfun_boost) {
    def("doubleit", doubleit);
    
    import_array();
    numeric::array::set_module_and_type("numpy", "ndarray");
    def("fill_array", fill_array);
}

