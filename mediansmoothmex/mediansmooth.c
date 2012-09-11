#include <float.h>

#define TYPE double
#define TYPE_MIN (-DBL_MAX)
#define TYPE_MAX DBL_MAX
#define elem_type elem_type_double
#define elem_ptr elem_ptr_double
#define print print_double
#define compar compar_double
#define median_init median_init_double
#define find_spot find_spot_double
#define median_add median_add_double
#define median_remove median_remove_double
#define median_replace median_replace_double
#define median_filt median_filt_double

#include "mediansmooth_template.c"
#include "undef.h"

#define TYPE float
#define TYPE_MIN (-FLT_MAX)
#define TYPE_MAX FLT_MAX
#define elem_type elem_type_float
#define elem_ptr elem_ptr_float
#define print print_float
#define compar compar_float
#define median_init median_init_float
#define find_spot find_spot_float
#define median_add median_add_float
#define median_remove median_remove_float
#define median_replace median_replace_float
#define median_filt median_filt_float

#include "mediansmooth_template.c"
#include "undef.h"


