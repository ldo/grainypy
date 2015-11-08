/*
    Lower-level code for grainy module.
*/

#include <iso646.h>
#include <cairo/cairo.h>
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdbool.h>
#include <stdint.h>

/*
    Miscellaneous useful stuff
*/

typedef unsigned int
    uint;

typedef unsigned char
    cpnt_t; /* pixel component type */
enum
  {
    CPNT_BITS = 8,
    CPNT_MAX = (1 << CPNT_BITS) - 1, /* max value of cpnt_t */
  };
typedef uint32_t
    pix_type;

typedef struct
  { /* identifies a pixel component within a pixemap */
    uint8_t * baseaddr; /* the base address of the pixels */
    uint width, height; /* the dimensions of the pixmap in pixels */
    uint depth; /* the number of bytes per pixel */
    uint stride; /* the number of bytes per row of pixels */
    uint shiftoffset; /* the offset in bits of the component from the bottom of the pixel */
    uint bitwidth; /* the width in bits of the pixel component */
  } channel;

static uint make_bayer_slice
  (
    uint order,
    uint * coeffs,
    uint offset,
    uint stride
  )
  /* recursively constructs a submatrix of a Bayer matrix. */
  {
    static const uint d2[4] = {3, 1, 0, 2};
    uint coeff_bits;
    if (order > 1)
      {
        if (order > 2)
          {
            uint row, col;
            make_bayer_slice(order / 2, coeffs, offset, stride);
            make_bayer_slice(order / 2, coeffs, offset + order / 2, stride);
            make_bayer_slice(order / 2, coeffs, offset + order / 2 * stride, stride);
            coeff_bits = make_bayer_slice(order / 2, coeffs, offset + order / 2 * (stride + 1), stride) + 2;
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
            coeff_bits = 2;
          } /*if*/
      }
    else
      {
        coeffs[offset] = 1;
        coeff_bits = 1;
      } /*if*/
    return
        coeff_bits;
  } /*make_bayer_slice*/

static uint make_bayer
  (
    uint order, /* must be power of 2 */
    uint * coeffs /* array[order * order] */
  )
  /* fills coeffs with a Bayer matrix of the specified order. */
  {
    return
        make_bayer_slice(order, coeffs, 0, order);
  } /*make_bayer*/

static void get_channel
  (
    PyObject * desc,
    channel * result
  )
  /* extracts the fields of desc, which is expected to be an instance of grainy.Channel.
    Does some basic consistency checks, and returns the info in the fields of result. */
  {
    PyObject * val = 0;
    do /*once*/
      {
        val = PyObject_GetAttrString(desc, "baseaddr");
        if (PyErr_Occurred())
            break;
        result->baseaddr = (uint8_t *)PyLong_AsSize_t(val);
        if (PyErr_Occurred())
            break;
#define get_field(name) \
        Py_XDECREF(val); \
        val = PyObject_GetAttrString(desc, #name); \
        if (PyErr_Occurred()) \
            break; \
        result->name = PyLong_AsLong(val); \
        if (PyErr_Occurred()) \
            break;
        get_field(width)
        get_field(height)
        get_field(depth)
        get_field(stride)
        get_field(shiftoffset)
        get_field(bitwidth)
#undef get_field
        if (result->shiftoffset + result->bitwidth > result->depth * 8)
          {
            PyErr_SetString(PyExc_IndexError, "pixmap channel position is outside pixel");
            break;
          } /*if*/
        if (result->width * result->depth > result->stride)
          {
            PyErr_SetString(PyExc_IndexError, "pixmap row width exceeds stride");
            break;
          } /*if*/
        if (result->shiftoffset % 8 != 0 or result->bitwidth != 8)
          {
            PyErr_SetString(PyExc_ValueError, "only aligned single-byte channels currently supported");
            break;
          } /*if*/
      }
    while (false);
    Py_XDECREF(val);
  } /*get_channel*/

/*
    User-visible stuff
*/

static PyObject * grainyx_bayer
  (
    PyObject * self,
    PyObject * args
  )
  {
    uint order;
    PyObject * result = 0;
    uint * coeffs = 0;
    uint coeff_bits;
    PyObject * coeffs_tuple = 0;
    PyObject * result_tuple = 0;
    uint row, col;
    do /*once*/
      {
        if (not PyArg_ParseTuple(args, "I", &order))
            break;
        if (order == 0 or (order - 1 & order) != 0)
          {
            PyErr_SetString(PyExc_ValueError, "order must be power of 2");
            break;
          } /*if*/
        coeffs = calloc(order * order, sizeof(uint));
        if (coeffs == 0)
          {
            PyErr_NoMemory();
            break;
          } /*if*/
        coeff_bits = make_bayer(order, coeffs);
        coeffs_tuple = PyTuple_New(order * order);
        if (coeffs_tuple == 0)
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
                    coeffs_tuple,
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
        result_tuple = PyTuple_New(2);
        if (result_tuple == 0)
            break;
        PyTuple_SET_ITEM(result_tuple, 0, PyLong_FromLong(coeff_bits));
        if (PyErr_Occurred())
            break;
        PyTuple_SET_ITEM(result_tuple, 1, coeffs_tuple);
        coeffs_tuple = 0; /* so I don't dispose of it yet */
      /* all done */
        result = result_tuple;
        result_tuple = 0; /* so I don't free it yet */
      }
    while (false);
    Py_XDECREF(coeffs_tuple);
    free(coeffs);
    return
        result;
  } /*grainyx_bayer*/

static PyObject * grainyx_ordered_dither
  (
    PyObject * self,
    PyObject * args
  )
  {
    PyObject * result = 0;
    cairo_surface_t * pix = 0;
    cairo_format_t pix_fmt;
    uint new_depth, dither_order, coeff_bits;
    bool do_a, do_r, do_g, do_b;
    uint * coeffs = 0;
    do /*once*/
      {
          {
            unsigned long pixaddr; /* why is there no Py_size_t? */
            uint i_do_a, i_do_r, i_do_g, i_do_b;
            do /*once*/
              {
                if
                  (
                    not PyArg_ParseTuple
                      (
                        args,
                        "kIIpppp",
                        &pixaddr, &new_depth, &dither_order, &i_do_a,  &i_do_r, &i_do_g, &i_do_b
                      )
                  )
                    break;
                if (new_depth > CPNT_BITS)
                  {
                    PyErr_Format(PyExc_ValueError, "depth cannot exceed %d", CPNT_BITS);
                    break;
                  } /*if*/
                if (dither_order > 1 and (dither_order - 1 & dither_order) != 0)
                  {
                    PyErr_SetString(PyExc_ValueError, "order must be power of 2");
                    break;
                  } /*if*/
                pix = (cairo_surface_t *)pixaddr;
                cairo_surface_reference(pix);
                if (cairo_surface_get_type(pix) != CAIRO_SURFACE_TYPE_IMAGE)
                  {
                    PyErr_SetString(PyExc_TypeError, "not an image surface");
                    break;
                  } /*if*/
                pix_fmt = cairo_image_surface_get_format(pix);
                if (pix_fmt != CAIRO_FORMAT_ARGB32 and pix_fmt != CAIRO_FORMAT_RGB24)
                  {
                    PyErr_SetString(PyExc_TypeError, "unsupported image surface format");
                    break;
                  } /*if*/
                do_a = i_do_a;
                if (do_a and pix_fmt != CAIRO_FORMAT_ARGB32)
                  {
                    PyErr_SetString(PyExc_TypeError, "cannot do_a with no alpha");
                    break;
                  } /*if*/
                do_r = i_do_r;
                do_g = i_do_g;
                do_b = i_do_b;
              }
            while (false);
          }
        if (PyErr_Occurred())
            break;
        if (not (do_a or do_r or do_g or do_b) or new_depth == CPNT_BITS)
          {
          /* nothing to do */
            Py_INCREF(Py_None);
            result = Py_None; /* return success */
            break;
          } /*if*/
          {
            uint drop_bits;
            cpnt_t * pix_base;
            uint width, height, row, col;
            cpnt_t drop_mask, keep_mask;
            size_t stride;
            if (dither_order > 1)
              {
                coeffs = calloc(dither_order * dither_order, sizeof(uint));
                if (coeffs == 0)
                  {
                    PyErr_NoMemory();
                    break;
                  } /*if*/
                coeff_bits = make_bayer(dither_order, coeffs);
              } /*if*/
            Py_BEGIN_ALLOW_THREADS
            drop_bits = CPNT_BITS - new_depth;
            drop_mask = (1 << drop_bits) - 1;
            keep_mask = ~drop_mask;
            cairo_surface_flush(pix);
            pix_base = cairo_image_surface_get_data(pix);
            stride = cairo_image_surface_get_stride(pix);
            width = cairo_image_surface_get_width(pix);
            height = cairo_image_surface_get_height(pix);
            for (row = 0; row != height; ++row)
              {
                pix_type * const row_base = (pix_type *)(pix_base + row * stride);
                for (col = 0; col != width; ++col)
                  {
                    pix_type pixel = row_base[col];
                  /* note: the following algorithm really only works correctly
                    when new_depth = 1. Otherwise it will produce colours that
                    are almost, but not quite, the same. This is because the
                    threshold comparison can only produce one of 2 values for
                    the bottommost rounded component bit. It needs to produce
                    (1 << new_depth) different values to map to all the rounded
                    component bits. This will require that number of threshold
                    coefficients at each pixel position. */
#define cond_do(doit, shift) \
                    if (doit) \
                      { \
                        pix_type val = pixel >> shift & CPNT_MAX; \
                        if (dither_order > 1) \
                          { \
                            uint const threshold = coeffs[row % dither_order * dither_order + col % dither_order]; \
                            uint frac = val & (drop_mask << 1) + 1; \
                            frac = \
                                drop_bits > coeff_bits ? \
                                    frac >> drop_bits - coeff_bits \
                                : \
                                    frac << coeff_bits - drop_bits; \
                            frac = frac > threshold ? CPNT_MAX : 0; \
                            val = val & keep_mask | frac & drop_mask; \
                          } \
                        else \
                          { \
                            val = (val & keep_mask) * CPNT_MAX / (CPNT_MAX - drop_mask); \
                              /* normalize truncated value to full brightness */ \
                          } /*if*/ \
                        pixel = pixel & ~(CPNT_MAX << shift) | val << shift; \
                      } /*if*/
                    cond_do(do_a, 3 * CPNT_BITS)
                    cond_do(do_r, 2 * CPNT_BITS)
                    cond_do(do_g, CPNT_BITS)
                    cond_do(do_b, 0)
#undef cond_do
                    row_base[col] = pixel;
                  } /*for*/
              } /*for*/
            cairo_surface_mark_dirty(pix);
            Py_END_ALLOW_THREADS
          }
      /* all successfully done */
        Py_INCREF(Py_None);
        result = Py_None;
      }
    while (false);
    free(coeffs);
    cairo_surface_destroy(pix);
    return
        result;
  } /*grainyx_ordered_dither*/

static PyObject * grainyx_copy_channel
  (
    PyObject * self,
    PyObject * args
  )
  {
    PyObject * result = 0;
    channel src, dst;
    do /*once*/
      {
          {
            PyObject * srcobj = 0;
            PyObject * dstobj = 0;
            do /*once*/
              {
                if (not PyArg_ParseTuple(args, "OO", &srcobj, &dstobj))
                    break;
                Py_INCREF(srcobj);
                Py_INCREF(dstobj);
                get_channel(srcobj, &src);
                if (PyErr_Occurred())
                    break;
                get_channel(dstobj, &dst);
                if (PyErr_Occurred())
                    break;
              }
            while (false);
            Py_XDECREF(srcobj);
            Py_XDECREF(dstobj);
          }
        if (PyErr_Occurred())
            break;
        if (src.width != dst.width or src.height != dst.height)
          {
            PyErr_SetString(PyExc_ValueError, "src/dst pixel dimensions mismatch");
            break;
          } /*if*/
        if (src.bitwidth != dst.bitwidth)
          {
            PyErr_SetString(PyExc_ValueError, "change in channel depth not currently supported");
            break;
          } /*if*/
        if (src.depth != 1 and src.depth != 4 or dst.depth != 1 and dst.depth != 4)
          {
            PyErr_SetString(PyExc_ValueError, "only pixel depths of 1 and 4 bytes currently supported");
            break;
          } /*if*/
        Py_BEGIN_ALLOW_THREADS
          {
            uint row, col;
            for (row = 0; row != src.height; ++row)
              {
                const uint8_t * const srcrow = src.baseaddr + row * src.stride;
                uint8_t * const dstrow = dst.baseaddr + row * dst.stride;
                for (col = 0; col != src.width; ++col)
                  {
                    uint val;
                    if (src.depth == 4)
                      {
                        val =
                                ((uint32_t *)srcrow)[col] >> src.shiftoffset
                            &
                                (1 << src.bitwidth) - 1;
                      }
                    else /* src.depth = 1 */
                      {
                      /* assume src.shiftoffset = 0 and src.bitwidth = 8 */
                        val = srcrow[col];
                      } /*if*/
                    if (dst.depth == 4)
                      {
                        uint dstval = ((uint32_t *)dstrow)[col];
                        dstval =
                                dstval & ~((1 << dst.bitwidth) - 1 << dst.shiftoffset)
                            |
                                val << dst.shiftoffset;
                        ((uint32_t *)dstrow)[col] = dstval;
                      }
                    else /* dst.depth = 1 */
                      {
                      /* assume dst.shiftoffset = 0 and dst.bitwidth = 8 */
                        dstrow[col] = val;
                      } /*if*/
                  } /*for*/
              } /*for*/
          }
        Py_END_ALLOW_THREADS
      /* all successfully done */
        Py_INCREF(Py_None);
        result = Py_None;
      }
    while (false);
    return
        result;
  } /*grainyx_copy_channel*/

static PyMethodDef grainyx_methods[] =
  {
    {"bayer", grainyx_bayer, METH_VARARGS,
        "bayer(order)\n"
        "returns a tuple(coeff_bits, coeffs) where coeffs is a tuple of (order * order)"
        " coefficients making up a Bayer ordered-dither matrix, and coeff_bits is the"
        " number of bits needed to hold each coefficient value. order must be a power of 2."
    },
    {"ordered_dither", grainyx_ordered_dither, METH_VARARGS,
        "ordered_dither(pix, depth, order, do_a, do_r, do_g, do_b)\n"
        "reduces the pixels in pix, which must be a pointer to a cairo_surface_t"
        " of FORMAT_RGB24 or FORMAT_ARGB32, down to depth bits per component."
        " order is size of the dither matrix to use; it must be a power of 2."
    },
    {"copy_channel", grainyx_copy_channel, METH_VARARGS,
        "copy_channel(src_channel, dst_channel)\n"
        "copies the contents of one grainy.Channel instance to another. You may copy"
        " between two channels of the same pixmap, but the pixmaps should not otherwise"
        " overlap."
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
