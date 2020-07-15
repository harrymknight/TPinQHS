#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include <stdbool.h>
#include <numpy/ndarraytypes.h>
#include <numpy/ndarrayobject.h>
#include "brent.h"
#include "ode.h"
#include "capsule.h"

//For ode subroutine
const int neqn;
double *y;
double *yp;
double *x;
double xout;
const double relerr;
const double abserr; //Absolute error of eigenvalue to be found
int *iflag;
double *work;
int *iwork;

//Properties of psi
double E;
double V;
double x_0; //Cyclotron orbit centre
int level;

double temp;
int dpoints;


double potential(double x_0, const double grad)
{
    return x_0 * grad;
}

void dpsi_dx(double x, double *y, double *yp, double E, double V)
{
    yp[0] = y[1];
    yp[1] = ((x * x) - 2 * (E - V)) * y[0];
}

void psi(void (*f)(double x, double *y, double *yp, double E, double V), double *x, double xout, double E, double V, double x_0, double *y, double *yp, const int neqn, double relerr, double abserr, int *iflag, double *work, int *iwork)
{
    //If at/beyond edge, crop wavefunction's domain. Extend right edge for large x_0 to allow decay of wavefunctions.
    if (x_0 < 0)
    {
        *x = -5 - x_0;
        xout = 5 - (x_0 * 0.35);
    }
    else
    {
        *x = -(x_0 * 0.35) - 5;
        xout = 5 - x_0;
    }
    f(*x, y, yp, E, V);
    ode(f, E, V, neqn, y, x, xout, relerr, abserr, iflag, work, iwork);
}

double try_eigenvalue(void (*f)(double x, double *y, double *yp, double E, double V), double *x, double xout, double E, double V, double x_0, const double igrad, double *y, double *yp, const int neqn, double relerr, double abserr, int *iflag, double *work, int *iwork)
{
    y[0] = 0;
    y[1] = igrad;
    *iflag = 1;
    psi(f, x, xout, E, V, x_0, y, yp, neqn, relerr, abserr, iflag, work, iwork);
    return y[0];
}

static PyObject *
get_eigenvalue(PyObject *self, PyObject *args)
{
    if (!PyArg_ParseTuple(args, "ddi", &x_0, &temp, &level))
        return NULL;
    const double grad = temp;
    if(level)
    {
        temp = -0.001;
    }
    else
    {
        temp = 0.001;
    }
    const double wfgrad = temp;

    //For ode subroutine
    const int neqn = 2;
    double y[neqn];
    double yp[neqn];
    double x[1];
    double xout = 0;
    const double relerr = 1e-5;
    const double abserr = 1e-6;
    int iflag = 1;
    double work[100 + 21 * neqn];
    int iwork[5];

    //Loop index
    int j;

    //Solving for eigenvalue (E)
    double a = 0;
    double b = 0;
    
    for (j = 0; j < 200; j++)
    {
        y[0] = 0;
        y[1] = wfgrad;
        iflag = 1;
        if (a < b && a > 0)
        {
            break;
        }
        psi(dpsi_dx, x, xout, j*0.1, potential(x_0, grad), x_0, y, yp, neqn, relerr, abserr, &iflag, work, iwork);
        if (y[0] > 0)
        {
            if (!a)
            {
                a = j*0.1;
            }
        }
        else
        {
            b = j*0.1;
        }
    }
    E = zero(a, b, 2.22e-16, abserr, try_eigenvalue, dpsi_dx, x, xout, potential(x_0, grad), x_0, wfgrad, y, yp, neqn, relerr, abserr, &iflag, work, iwork);
    return PyFloat_FromDouble(E);
}

static PyObject*
get_eigenvalues(PyObject *self, PyObject *args)
{
    PyObject *objx_k;
    if (!PyArg_ParseTuple(args, "Odi", &objx_k, &temp, &level))
        return NULL;
    const double grad = temp;
    if(level)
    {
        temp = -0.001;
    }
    else
    {
        temp = 0.001;
    }
    const double wfgrad = temp;

    //For ode subroutine
    const int neqn = 2;
    double y[neqn];
    double yp[neqn];
    double x[1];
    double xout = 0;
    const double relerr = 1e-5;
    const double abserr = 1e-6;
    int iflag = 1;
    double work[100 + 21 * neqn];
    int iwork[5];

    //Loop index
    int i;
    int j;

    int size = PyArray_SIZE(objx_k);
    double *x_k = PyArray_DATA(objx_k);
    double *eigenvalues = malloc(sizeof(double)*size);

    //solving for eigenvalues
    double a;
    double b;
    double past_E = 0;
    for(i=0; i<size; i++)
    {
        //Solving for eigenvalue corresponding to cyclotron orbit buffer[i]
        a = 0;
        b = 0;
        
        for (j = past_E; j < 200; j++)
        {
            y[0] = 0;
            y[1] = wfgrad;
            iflag = 1;
            if (a < b && a > 0)
            {
                break;
            }
            psi(dpsi_dx, x, xout, j*0.1, potential(x_k[i], grad), x_k[i], y, yp, neqn, relerr, abserr, &iflag, work, iwork);
            if (y[0] > 0)
            {
                if (!a)
                {
                    a = j*0.1;
                }
            }
            else
            {
                b = j*0.1;
            }
        }
        eigenvalues[i] = zero(a, b, 2.22e-16, abserr, try_eigenvalue, dpsi_dx, x, xout, potential(x_k[i], grad), x_k[i], wfgrad, y, yp, neqn, relerr, abserr, &iflag, work, iwork);
        past_E = (eigenvalues[i] - 1) > 0 ? : 0;
    }
    npy_intp const dims[1] = {size};
    PyObject *objeigenvalues = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, eigenvalues);
    PyObject *capeigenvalues = PyCapsule_New(eigenvalues, NULL, capsule_cleanup);
    PyArray_SetBaseObject((PyArrayObject *) objeigenvalues, capeigenvalues);    
    return objeigenvalues;
}

static PyObject *
get_wavefunc(PyObject *self, PyObject *args)
{
    if (!PyArg_ParseTuple(args, "dddin", &x_0, &E, &temp, &level, &dpoints))
        return NULL;
    const double grad = temp;

    typedef struct {
        double x;
        double y;
    } Wavefunc;

    Wavefunc wf[dpoints];

    if(level)
    {
        temp = -0.001;
    }
    else
    {
        temp = 0.001;
    }
    double wfgrad;

    //For ode subroutine
    const int neqn = 2;
    double y[neqn];
    double yp[neqn];
    double x[1];
    double xout = 0;
    const double relerr = 1e-5;
    const double abserr = 1e-6;
    int iflag = 1;
    double work[100 + 21 * neqn];
    int iwork[5];

    //loop indices
    int i = 0;

    double xmax;
    double interval;
    if (x_0 < 0)
    {
        *x = 5 - (x_0 * 0.35);
        xmax = -5 - x_0;
        interval = *x - xmax;
        wfgrad = -temp;
    }
    else
    {
        *x = -(x_0 * 0.35) - 5;
        xmax = 5 - x_0;
        interval = xmax - *x;
        wfgrad = temp;
    }
    double spacing = interval/((dpoints-1));
    xout = *x;
    wf[0].y = 0;
    y[0] = 0;
    y[1] = wfgrad;
    double sum;
    if(x_0 < 0)
    {
        for(i=dpoints-1;i>-1;i--)
        {
            wf[i].x = xout;
            dpsi_dx(*x, y, yp, E, V);
            xout -= spacing;
            iflag = 1;
            ode(dpsi_dx, E, potential(x_0, grad), neqn, y, x, xout, relerr, abserr, &iflag, work, iwork);
            sum += y[0]*y[0]*spacing;
            wf[i].y = y[0];
        }
    }
    else
    {
        for(i=0;i<dpoints;i++)
        {
            wf[i].x = xout;
            dpsi_dx(*x, y, yp, E, V);
            xout += spacing;
            iflag = 1;
            ode(dpsi_dx, E, potential(x_0, grad), neqn, y, x, xout, relerr, abserr, &iflag, work, iwork);
            sum += y[0]*y[0]*spacing;
            wf[i].y = y[0];
        }
    }
    double _sum = 1/sum;
    _sum = sqrt(_sum);
    for(i=0;i<dpoints;i++)
    {
        wf[i].y *= _sum;
    }
    double *flatwf = malloc(sizeof(double)*(2*dpoints));
    for(i=0;i<2*dpoints;i+=2)
    {
        flatwf[i] = wf[i/2].x;
        flatwf[i+1] = wf[i/2].y;
    }
    npy_intp const dims[] = {dpoints};
    PyObject *temp;
    PyArray_Descr *descr;
    temp = Py_BuildValue("[(s, s), (s, s)]", "x", "f8", "y", "f8");
    PyArray_DescrConverter(temp, &descr);
    Py_DECREF(temp);
    PyObject *objwf = PyArray_NewFromDescr(&PyArray_Type, descr, 1, dims, NULL, flatwf, NPY_ARRAY_CARRAY, NULL);
    PyObject *capwf = PyCapsule_New(flatwf, NULL, capsule_cleanup);
    PyArray_SetBaseObject((PyArrayObject *) objwf, capwf);
    return objwf;
}

static PyObject *
get_x_expectation(PyObject *self, PyObject *args)
{
    PyObject *objwf;
    if (!PyArg_ParseTuple(args, "O", &objwf))
        return NULL;

    PyObject *objwfx, *objwfy;
    objwfx = PyArray_GetField((PyArrayObject *) objwf, PyArray_DescrFromType(NPY_FLOAT64), 0);
    objwfy = PyArray_GetField((PyArrayObject *) objwf, PyArray_DescrFromType(NPY_FLOAT64), 8);

    int i;
    int size = PyArray_SIZE(objwf);
    double *iwfy, *iwfx0, *iwfx1;
    double sum = 0;
    for(i=1;i<size;i++)
    {
        iwfy = PyArray_GETPTR1((PyArrayObject *) objwfy, i);
        iwfx0 = PyArray_GETPTR1((PyArrayObject *) objwfx, i-1);
        iwfx1 = PyArray_GETPTR1((PyArrayObject *) objwfx, i);
        sum += (*iwfy) * (*iwfy) * (*iwfx1) * (*iwfx1-*iwfx0);
    }
    return PyFloat_FromDouble(sum);
}

static PyMethodDef WavefuncMethods[] = {
    {
     "get_eigenvalues",
      get_eigenvalues, 
      METH_VARARGS,
     "Get energies (eigenvalue) of a set of wavefunctions given an iterable containing cyclotron orbit centres and the Landau level they occupy."
    },
    {
     "get_eigenvalue",
      get_eigenvalue, 
      METH_VARARGS,
     "Get energy (eigenvalue) of wavefunction given its cyclotron orbit centre and the Landau level it occupies."
    },
    {
     "get_wavefunc",
      get_wavefunc, 
      METH_VARARGS,
     "Get wavefunction given its cyclotron orbit centre, energy, and number of points to plot."
    },
    {
     "get_x_expectation",
      get_x_expectation, 
      METH_VARARGS,
     "Get average position of wavefunction given its cyclotron orbit centre, energy, and number of points to plot."
    },
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef wavefuncmodule = {
    PyModuleDef_HEAD_INIT,
    "wavefunc",
    "Python interface for routines which solve for properties of wavefunctions",
    -1,
    WavefuncMethods
};

PyMODINIT_FUNC PyInit_wavefunc(void) {
    import_array();
    return PyModule_Create(&wavefuncmodule);
}
