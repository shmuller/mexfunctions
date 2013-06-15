#include <boost/python.hpp>
#include <boost/numpy.hpp>

#include "../../fitfun.h"

namespace p = boost::python;
namespace np = boost::numpy;

void parse_args(data *D, np::ndarray &P, np::ndarray &x, np::ndarray &y)
{
    D->P = (double*) P.get_data();
    D->n = P.shape(0);

    D->x = (double*) x.get_data();
    D->m = x.shape(0);

    D->y = (double*) y.get_data();
}

void parse_args_a(data *D, np::ndarray &P, np::ndarray &x, np::ndarray &y, 
                           np::ndarray &a)
{
    parse_args(D, P, x, y);
    D->a = (double*) a.get_data();
}

void parse_args_ydata(data *D, np::ndarray &P, np::ndarray &x, np::ndarray &y, 
                               np::ndarray &ydata)
{
    parse_args(D, P, x, y);
    D->ydata = (double*) ydata.get_data();
}


np::ndarray meth_e2(np::ndarray &P, np::ndarray &x, np::ndarray &y)
{
    data D = {0};
    parse_args(&D, P, x, y);
    e2(&D);
    return y;
}

np::ndarray meth_e2_diff(np::ndarray &P, np::ndarray &x, np::ndarray &y, np::ndarray &ydata)
{
    data D = {0};
    parse_args_ydata(&D, P, x, y, ydata);
    e2_diff(&D);
    return y;
}

double meth_e2_rms(np::ndarray &P, np::ndarray &x, np::ndarray &y, np::ndarray &ydata)
{
    data D = {0};
    parse_args_ydata(&D, P, x, y, ydata);
    return e2_rms(&D);
}

np::ndarray meth_e2_fit(np::ndarray &P, np::ndarray &x, np::ndarray &y, np::ndarray &ydata)
{
    data D = {0};
    parse_args_ydata(&D, P, x, y, ydata);
    e2_fit(&D);
    return P;
}


BOOST_PYTHON_MODULE(fitfun_boost_numpy) {
    np::initialize();
    p::def("e2", meth_e2);
    p::def("e2_diff", meth_e2_diff);
    p::def("e2_rms", meth_e2_rms);
    p::def("e2_fit", meth_e2_fit);
}

