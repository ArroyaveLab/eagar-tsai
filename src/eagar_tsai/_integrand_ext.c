/*
 * _integrand_ext.c — Python C extension exposing the Eagar–Tsai integrand
 * as a scipy.integrate.LowLevelCallable PyCapsule.
 *
 * The capsule name "double (int, double *)" is one of the four signatures
 * recognised by SciPy's QUADPACK wrapper, enabling zero-Python-overhead
 * evaluation during numerical integration.
 *
 * Build: handled automatically by setup.py / uv sync.
 *
 * Reference: T. W. Eagar and N.-S. Tsai, Welding Journal, December 1983.
 * Reformulation: Sasha Rubenchik, LLNL, 2015.
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>

/* ------------------------------------------------------------------ */
/* Core integrand                                                       */
/* ------------------------------------------------------------------ */

/*
 * _et_integrand(n, args)
 *
 * args layout (indices 0–4):
 *   args[0] = t   integration variable (dimensionless time), t > 0
 *   args[1] = x   non-dimensional x (scan direction)
 *   args[2] = y   non-dimensional y (cross-scan direction)
 *   args[3] = z   non-dimensional z (depth)
 *   args[4] = p   alpha / (v * sigma)
 *
 * Returns:
 *   f = 1 / ((4p*t + 1) * sqrt(t))
 *       * exp( -z^2/(4t) - (y^2 + (x - t)^2) / (4p*t + 1) )
 */
static double _et_integrand(int n, double *args)
{
    double t       = args[0];
    double x       = args[1];
    double y       = args[2];
    double z       = args[3];
    double p       = args[4];
    double denom_xy = 4.0 * p * t + 1.0;
    double dx      = x - t;

    return (1.0 / (denom_xy * sqrt(t)))
         * exp(-(z * z) / (4.0 * t)
               - ((y * y) + dx * dx) / denom_xy);
}

/* ------------------------------------------------------------------ */
/* Python-callable: return a PyCapsule wrapping _et_integrand          */
/* ------------------------------------------------------------------ */

static PyObject *
get_integrand_capsule(PyObject *self, PyObject *args)
{
    /* The capsule name must match one of SciPy's accepted LLC signatures. */
    return PyCapsule_New(
        (void *)_et_integrand,
        "double (int, double *)",
        NULL
    );
}

/* ------------------------------------------------------------------ */
/* Module boilerplate                                                   */
/* ------------------------------------------------------------------ */

static PyMethodDef _methods[] = {
    {
        "get_integrand_capsule",
        get_integrand_capsule,
        METH_NOARGS,
        "Return a PyCapsule wrapping the Eagar-Tsai integrand.\n\n"
        "Pass the capsule to scipy.integrate.LowLevelCallable for\n"
        "zero-overhead numerical integration.\n\n"
        "Returns:\n"
        "    capsule: PyCapsule with signature 'double (int, double *)'\n"
    },
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef _module = {
    PyModuleDef_HEAD_INIT,
    "_integrand_ext",
    "C extension exposing the Eagar-Tsai integrand as a LowLevelCallable.",
    -1,
    _methods
};

PyMODINIT_FUNC
PyInit__integrand_ext(void)
{
    return PyModule_Create(&_module);
}
