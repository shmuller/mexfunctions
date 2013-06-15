#include <boost/python.hpp>
#include <boost/python/raw_function.hpp>
#include <boost/numpy.hpp>

#include "../../fitfun.h"

namespace p = boost::python;
namespace np = boost::numpy;

void parse_args(p::tuple &args, data *D)
{
    np::ndarray P = np::array(args[0]);
    D->P = (double*) P.get_data();
    D->n = P.shape(0);

    np::ndarray x = np::array(args[1]);
    D->x = (double*) x.get_data();
    D->m = x.shape(0);

    np::ndarray y = np::array(args[2]);
    D->y = (double*) y.get_data();
}

void parse_args_a(p::tuple &args, data *D)
{
    parse_args(args, D);
    np::ndarray a = np::array(args[3]);
    D->a = (double*) a.get_data();
}

void parse_args_ydata(p::tuple &args, data *D)
{
    parse_args(args, D);
    np::ndarray ydata = np::array(args[3]);
    D->ydata = (double*) ydata.get_data();
}

void parse_args_a_ydata(p::tuple &args, data *D)
{
    parse_args_ydata(args, D);
    np::ndarray a = np::array(args[4]);
    D->a = (double*) a.get_data();
}


#define meth_template_passthru(fun, parser, i)             \
static np::ndarray meth_##fun(p::tuple args, p::dict kw) { \
    data D = {0};                                          \
    parser(args, &D);                                      \
    fun(&D);                                               \
    return np::array(args[i]);                             \
}

#define meth_template_double(fun, parser)             \
static double meth_##fun(p::tuple args, p::dict kw) { \
    data D = {0};                                     \
    parser(args, &D);                                 \
    return fun(&D);                                   \
}

#define meth_template_passthru_fit(fun, parser, i)         \
static np::ndarray meth_##fun(p::tuple args, p::dict kw) { \
    data D = {0};                                          \
    parser(args, &D);                                      \
    if (kw.has_key("do_var")) {                            \
        np::ndarray do_var = np::array(kw.get("do_var"));  \
        D.do_var = (int*) do_var.get_data();               \
    }                                                      \
    fun(&D);                                               \
    return np::array(args[i]);                             \
}


#define meth_template_all(fun, parser)                   \
meth_template_passthru(fun, parser, 2)                   \
meth_template_passthru(fun##_diff, parser##_ydata, 2)    \
meth_template_passthru_fit(fun##_fit, parser##_ydata, 0) \
meth_template_double(fun##_rms, parser)

meth_template_all(e2, parse_args)
meth_template_all(IV3, parse_args)
meth_template_all(IV4, parse_args_a)
meth_template_all(IV5, parse_args_a)
meth_template_all(IV6, parse_args_a)
meth_template_all(IVdbl, parse_args)
meth_template_all(IVdbl2, parse_args_a)


#define def(fun)                                      \
p::def(#fun, raw_function(meth_##fun));               \
p::def(#fun"_diff", raw_function(meth_##fun##_diff)); \
p::def(#fun"_rms", raw_function(meth_##fun##_rms));   \
p::def(#fun"_fit", raw_function(meth_##fun##_fit));


BOOST_PYTHON_MODULE(fitfun_boost_numpy) {
    np::initialize();
    
    def(e2)
    def(IV3)
    def(IV4)
    def(IV5)
    def(IV6)
    def(IVdbl)
    def(IVdbl2)
}



