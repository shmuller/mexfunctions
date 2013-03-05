#include <minpack.h>

#define template_proto(fun) \
void fun(data *D);          \
void fun##_diff(data *D);   \
double fun##_rms(data *D);  \
void fun##_fit(data *D);

template_proto(e2)
template_proto(IV3)
template_proto(IV4)
template_proto(IV5)
template_proto(IV6)
template_proto(IVdbl)
template_proto(IVdbl2)

