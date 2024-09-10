/* The batman package: fast computation of exoplanet transit light curves
 * Copyright (C) 2015 Laura Kreidberg	 
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Minor edits for use with limb-brighted, X-ray transits distributed under the same terms
 * Copyright (C) 2024 George King
 */

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include "numpy/arrayobject.h"
#include "common.h"

static PyObject *_custom_ld(PyObject *self, PyObject *args);

/*
	- The intensity function returns the stellar intensity at a radius x; where 0 <= x <= 1
	- The function arguments are the normalized radius (x), and limb darkening coefficients c1, ..., un
	- see below for an example intensity profile for I(x) \propto 1 - c1*(1-sqrt(1-x^2)) - c2*ln((sqrt(1-x^2)+c)/(1+c))
	- The normalization constant is calculated by constraining the integrated intensity to equal 1:
		\int_x \int_theta {I(x)*x*dx*dtheta}/norm = 1
*/
inline double intensity(double x, double intVal[10001])
{
	if(x > 0.99995) x = 0.99995;
	if(x < 0.00005) x = 0.00005;
	
	int i=0;
        double j=0;
	double x1;

        x1 = x * 10000.;
	i = ceil( x1 );
	j = i / 10000.;

	double interpol = (intVal[i] - intVal[i-1]) * ( x - (j-(1./10000.)) ) / (1. / 10000.);
	return intVal[i-1] + interpol;
}


static PyObject *_custom_ld(PyObject *self, PyObject *args)
{
	double rprs, fac, c2, c3, c4, c5, c6;
	int nthreads;
	npy_intp dims[1];
	PyArrayObject *ds, *flux, *c1;

  	if(!PyArg_ParseTuple(args,"OdOddddddi", &ds, &rprs, &c1, &c2, &c3, &c4, &c5, &c6, &fac, &nthreads)) return NULL; //parses input arguments
	
	dims[0] = PyArray_DIMS(ds)[0]; 
	flux = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_TYPE(ds));	//creates numpy array to store return flux values

	double *f_array = PyArray_DATA(flux);
	double *d_array = PyArray_DATA(ds);
	double *intVals = PyArray_DATA(c1);

	/*
		NOTE:  the safest way to access numpy arrays is to use the PyArray_GETITEM and PyArray_SETITEM functions.
		Here we use a trick for faster access and more convenient access, where we set a pointer to the 
		beginning of the array with the PyArray_DATA (e.g., f_array) and access elements with e.g., f_array[i].
		Success of this operation depends on the numpy array storing data in blocks equal in size to a C double.
		If you run into trouble along these lines, I recommend changing the array access to something like:
			d = PyFloat_AsDouble(PyArray_GETITEM(ds, PyArray_GetPtr(ds, &i))); 
		where ds is a numpy array object.


		Laura Kreidberg 07/2015
	*/

	#pragma acc data copyin(intensity_args)
	calc_limb_darkening(f_array, d_array, dims[0], rprs, fac, nthreads, intVals);
	return PyArray_Return((PyArrayObject *)flux);
} 

static char _custom_ld_doc[] = "This extension module returns a limb darkened light curve for a custom stellar intensity profile.";


static PyMethodDef _custom_ld_methods[] = {
  {"_custom_ld", _custom_ld, METH_VARARGS, _custom_ld_doc},{NULL}};

#if PY_MAJOR_VERSION >= 3
	static struct PyModuleDef _custom_ld_module = {
		PyModuleDef_HEAD_INIT,
		"_custom_ld",
		_custom_ld_doc,
		-1, 
		_custom_ld_methods
	};

	PyMODINIT_FUNC
	PyInit__custom_ld(void)
	{
		PyObject* module = PyModule_Create(&_custom_ld_module);
		if(!module)
		{
			return NULL;
		}
		import_array(); 
		return module;
	}
#else

	void init_custom_ld(void)
	{
		Py_InitModule("_custom_ld", _custom_ld_methods);
		import_array(); 
	}
#endif

