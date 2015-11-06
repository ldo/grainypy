/*
    Lower-level code for grainy module.
*/

#include <iso646.h>
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdbool.h>

/*
    Miscellaneous useful stuff
*/

typedef unsigned int
    uint;

typedef unsigned char
    cpnt_t; /* pixel component type */
enum
  {
    CPNT_MAX = 255, /* max value of cpnt_t */
  };

static void make_bayer_slice
  (
    uint order,
    cpnt_t * coeffs,
    uint offset,
    uint stride
  )
  /* recursively constructs a submatrix of a Bayer matrix. */
  {
    static const cpnt_t d2[4] = {3, 1, 0, 2};
    if (order > 1)
      {
        if (order > 2)
          {
            uint row, col;
            make_bayer_slice(order / 2, coeffs, offset, stride);
            make_bayer_slice(order / 2, coeffs, offset + order / 2, stride);
            make_bayer_slice(order / 2, coeffs, offset + order / 2 * stride, stride);
            make_bayer_slice(order / 2, coeffs, offset + order / 2 * (stride + 1), stride);
            for (row = 0; row != order; ++row)
              {
                for (col = 0; col != order; ++col)
                  {
                    coeffs[offset + row * stride + col] =
                            4 * coeffs[offset + row * stride + col]
                        +
                            d2[(row * 2 / order) * 2 + col * 2 / order];
                  } /*for*/
              } /*for*/
          }
        else
          {
            coeffs[offset] = d2[0];
            coeffs[offset + 1] = d2[1];
            coeffs[offset + stride] = d2[2];
            coeffs[offset + stride + 1] = d2[3];
          } /*if*/
      }
    else
      {
        coeffs[offset] = 1;
      } /*if*/
  } /*make_bayer_slice*/

static void make_bayer
  (
    uint order, /* must be power of 2 */
    cpnt_t * coeffs /* array[order * order] */
  )
  /* fills coeffs with a Bayer matrix of the specified order. */
  {
    make_bayer_slice(order, coeffs, 0, order);
  } /*make_bayer*/

/*
    User-visible stuff
*/

static PyObject * grainyx_bayer
  (
    PyObject * self,
    PyObject * args
  )
  /* returns a tuple of coefficients making up a Bayer ordered-dither matrix. */
  {
    uint order;
    PyObject * result = 0;
    cpnt_t * coeffs = 0;
    PyObject * result_temp = 0;
    uint row, col;
    do /*once*/
      {
        if (not PyArg_ParseTuple(args, "I", &order))
            break;
        if ((order - 1 & order) != 0)
          {
            PyErr_SetString(PyExc_ValueError, "order must be power of 2");
            break;
          } /*if*/
        coeffs = calloc(order * order, sizeof(cpnt_t));
        if (coeffs == 0)
          {
            PyErr_NoMemory();
            break;
          } /*if*/
        make_bayer(order, coeffs);
        result_temp = PyTuple_New(order * order);
        if (result_temp == 0)
            break;
        for (row = 0;;)
          {
            if (row == order)
                break;
            for (col = 0;;)
              {
                if (col == order)
                    break;
                PyTuple_SET_ITEM
                  (
                    result_temp,
                    row * order + col,
                    PyLong_FromLong(coeffs[row * order + col])
                  );
                if (PyErr_Occurred())
                    break;
                ++col;
              } /*for*/
            if (PyErr_Occurred())
                break;
            ++row;
          } /*for*/
        if (PyErr_Occurred())
            break;
      /* all done */
        result = result_temp;
        result_temp = 0; /* so I don't free it yet */
      }
    while (false);
    Py_XDECREF(result_temp);
    free(coeffs);
    return result;
  } /*grainyx_bayer*/

static PyMethodDef grainyx_methods[] =
  {
    {"bayer", grainyx_bayer, METH_VARARGS,
        "bayer(order)\n"
        "returns a tuple of (order * order) coefficients making up a Bayer ordered-dither"
        " matrix. order must be a power of 2."
    },
    {0, 0, 0, 0} /* marks end of list */
  };
static PyModuleDef grainyx_module =
  {
    PyModuleDef_HEAD_INIT,
    "grainyx", /* module name */
    "lower-level functions for grainy", /* docstring */
    -1, /* size of per-interpreter state, -1 if entirely global */
    grainyx_methods,
  };

PyMODINIT_FUNC PyInit_grainyx(void)
  {
    return
        PyModule_Create(&grainyx_module);
  } /*PyInit_grainyx*/
