/* specfunmex(name, x)
 * Octave interface for specfun.f90. Compile with:
 *
 * mkoctfile specfunoct.cc -lcommon -lspecfun
 *
 * ln -s specfunoct.oct specfunmex.oct
 *
 * S. H. Muller, 2012/02/14
 */

#include <octave/oct.h>

#include "../specfun.h"

DEFUN_DLD(specfunmex, args, nargout, "Special functions")
{
	octave_value retval;

	const charMatrix C = args(0).char_matrix_value();
	const char *name = C.data();

	// trick octave to not make a copy
	const NDArray X = args(1).array_value();
	double *x = (double*) X.data();

	func *fun = (func*) kv_select(KV_LEN(kv_specfun), kv_specfun, name);
	if (fun == NULL) {
		error("specfunmex: Unknown function name");
		return retval;
	}

	for(int i=X.nelem(); i--; ) {
		*x++ = fun(x);
	}

	return retval;
}


