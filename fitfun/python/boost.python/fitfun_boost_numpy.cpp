#include <boost/python.hpp>
#include <boost/numpy.hpp>

#include "../../fitfun.h"

namespace p = boost::python;
namespace np = boost::numpy;

np::ndarray None;

void parse_args(data *D, np::ndarray &P, np::ndarray &x, np::ndarray &y, 
                         np::ndarray &ydata, np::ndarray &a)
{
    D->P = (double*) P.get_data();
    D->n = P.shape(0);

    D->x = (double*) x.get_data();
    D->m = x.shape(0);

    D->y = (double*) y.get_data();

    if (ydata != None) D->ydata = (double*) ydata.get_data();
    if (a != None) D->a = (double*) a.get_data();
}


np::ndarray meth_e2(np::ndarray &P, np::ndarray &x, np::ndarray &y, 
                    np::ndarray &ydata=None, np::ndarray &a=None)
{
    data D = {0};
    parse_args(&D, P, x, y, ydata, a);
    e2(&D);
    return y;
}

np::ndarray meth_e2_diff(np::ndarray &P, np::ndarray &x, np::ndarray &y, 
                         np::ndarray &ydata=None, np::ndarray &a=None)
{
    data D = {0};
    parse_args(&D, P, x, y, ydata, a);
    e2_diff(&D);
    return y;
}

double meth_e2_rms(np::ndarray &P, np::ndarray &x, np::ndarray &y, 
                   np::ndarray &ydata=None, np::ndarray &a=None)
{
    data D = {0};
    parse_args(&D, P, x, y, ydata, a);
    return e2_rms(&D);
}

np::ndarray meth_e2_fit(np::ndarray &P, np::ndarray &x, np::ndarray &y, 
                        np::ndarray &ydata=None, np::ndarray &a=None)
{
    data D = {0};
    parse_args(&D, P, x, y, ydata, a);
    e2_fit(&D);
    return P;
}

BOOST_PYTHON_FUNCTION_OVERLOADS(meth_e2_overloads, meth_e2, 3, 5)
BOOST_PYTHON_FUNCTION_OVERLOADS(meth_e2_diff_overloads, meth_e2_diff, 3, 5)
BOOST_PYTHON_FUNCTION_OVERLOADS(meth_e2_rms_overloads, meth_e2_rms, 3, 5)
BOOST_PYTHON_FUNCTION_OVERLOADS(meth_e2_fit_overloads, meth_e2_fit, 3, 5)


BOOST_PYTHON_MODULE(fitfun_boost_numpy) {
    np::initialize();

    p::def("e2", meth_e2, meth_e2_overloads());
    p::def("e2_diff", meth_e2_diff, meth_e2_diff_overloads());
    p::def("e2_rms", meth_e2_rms, meth_e2_rms_overloads());
    p::def("e2_fit", meth_e2_fit, meth_e2_fit_overloads());
}

