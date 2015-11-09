/*
    Lower-level code for grainy module. This C code implements the most
    CPU-intensive parts of the algorithms. It is written to make as few
    assumptions about data structures as possible, leaving it to the
    Python layer to make things easier for the caller.
*/

#include <iso646.h>
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
  } channel_t;

typedef struct
  {
    uint order, bits;
    uint coeffs[/*order * order */];
  } dither_matrix_t;

static bool get_channel
  (
    PyObject * desc,
    bool required, /* desc cannot be None for no channel */
    channel_t * same_pix_as, /* if non-null, then result must have same baseaddr, width, height, stride and depth as this */
    channel_t * result
  )
  /* extracts the fields of desc, which is expected to be an instance of grainy.Channel.
    Does some basic consistency checks, and returns the info in the fields of result. */
  {
    bool got_channel = false; /* to begin with */
    if (required or desc != 0 and desc != Py_None)
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
            if
              (
                    same_pix_as != 0
                and
                    (
                        result->baseaddr != same_pix_as->baseaddr
                    or
                        result->width != same_pix_as->width
                    or
                        result->height != same_pix_as->height
                    or
                        result->stride != same_pix_as->stride
                    or
                        result->depth != same_pix_as->depth
                    )
              )
              {
                PyErr_SetString(PyExc_ValueError, "inconsistent pixmap for channel");
                break;
              } /*if*/
            got_channel = true;
          }
        while (false);
        Py_XDECREF(val);
      }
    return
        got_channel;
  } /*get_channel*/

static void get_dither_matrix
  (
    PyObject * obj,
    bool required, /* obj cannot be None for no matrix */
    dither_matrix_t ** mat
  )
  /* extracts the fields from a grainy.DitherMatrix instance and allocates and fills in
    a new dither_matrix_t structure with the results. */
  {
    *mat = 0; /* to begin with */
    if (required or obj != 0 and obj != Py_None)
      {
        dither_matrix_t * result = 0;
        PyObject * field = 0;
        uint order, row, col;
        do /*once*/
          {
            field = PyObject_GetAttrString(obj, "order");
            if (PyErr_Occurred())
                break;
            order = PyLong_AsUnsignedLong(field);
            if (PyErr_Occurred())
                break;
            result = malloc(sizeof(dither_matrix_t) + sizeof(uint) * order * order);
            if (result == 0)
              {
                PyErr_NoMemory();
                break;
              } /*if*/
            result->order = order;
            Py_XDECREF(field);
            field = PyObject_GetAttrString(obj, "bits");
            if (PyErr_Occurred())
                break;
            result->bits = PyLong_AsUnsignedLong(field);
            if (PyErr_Occurred())
                break;
            Py_XDECREF(field);
            field = PyObject_GetAttrString(obj, "coeffs");
            if (PyErr_Occurred())
                break;
            for (row = 0;;)
              {
                if  (row == order)
                    break;
                for (col = 0;;)
                  {
                    if (col == order)
                        break;
                    const uint index = row * order + col;
                    PyObject * const elt = PyTuple_GetItem(field, index);
                      /* borrowed reference, valid as long as field is valid */
                    if (PyErr_Occurred())
                        break;
                    result->coeffs[index] = PyLong_AsUnsignedLong(elt);
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
            *mat = result;
            result = 0; /* so I don't dispose of it yet */
          }
        while (false);
        free(result);
        Py_XDECREF(field);
      } /*if*/
  } /*get_dither_matrix*/

/*
    User-visible stuff
*/

static PyObject * grainyx_ordered_dither
  (
    PyObject * self,
    PyObject * args
  )
  {
    PyObject * result = 0;
    dither_matrix_t * dither = 0;
    uint to_depth;
    channel_t srcchan1, dstchan1, srcchan2, dstchan2, srcchan3, dstchan3, srcchan4, dstchan4;
    bool gotchan1, gotchan2, gotchan3, gotchan4;
    channel_t * some_src;
    channel_t * some_dst;
    do /*once*/
      {
          { /* get args */
            PyObject * matrixobj = 0;
            PyObject * srcchan1obj = 0;
            PyObject * dstchan1obj = 0;
            PyObject * srcchan2obj = 0;
            PyObject * dstchan2obj = 0;
            PyObject * srcchan3obj = 0;
            PyObject * dstchan3obj = 0;
            PyObject * srcchan4obj = 0;
            PyObject * dstchan4obj = 0;
            do /*once*/
              {
                if
                  (
                    not PyArg_ParseTuple
                      (
                        args,
                        "OIOO|OOOOOO",
                        &matrixobj,
                        &to_depth,
                        &srcchan1obj, &dstchan1obj,
                        &srcchan2obj, &dstchan2obj,
                        &srcchan3obj, &dstchan3obj,
                        &srcchan4obj, &dstchan4obj
                      )
                  )
                    break;
              /* replace omitted arguments with None to simplify checks */
                if (srcchan2obj == 0)
                  {
                    srcchan2obj = Py_None;
                  } /*if*/
                if (dstchan2obj == 0)
                  {
                    dstchan2obj = Py_None;
                  } /*if*/
                if (srcchan3obj == 0)
                  {
                    srcchan3obj = Py_None;
                  } /*if*/
                if (dstchan3obj == 0)
                  {
                    dstchan3obj = Py_None;
                  } /*if*/
                if (srcchan4obj == 0)
                  {
                    srcchan4obj = Py_None;
                  } /*if*/
                if (dstchan4obj == 0)
                  {
                    dstchan4obj = Py_None;
                  } /*if*/
                Py_INCREF(matrixobj);
                Py_INCREF(srcchan1obj);
                Py_INCREF(dstchan1obj);
                Py_INCREF(srcchan2obj);
                Py_INCREF(dstchan2obj);
                Py_INCREF(srcchan3obj);
                Py_INCREF(dstchan3obj);
                Py_INCREF(srcchan4obj);
                Py_INCREF(dstchan4obj);
                if
                  (
                        (srcchan1obj != Py_None) != (dstchan1obj != Py_None)
                    or
                        (srcchan2obj != Py_None) != (dstchan2obj != Py_None)
                    or
                        (srcchan3obj != Py_None) != (dstchan3obj != Py_None)
                    or
                        (srcchan4obj != Py_None) != (dstchan4obj != Py_None)
                  )
                  {
                    PyErr_SetString
                      (
                        PyExc_ValueError,
                        "src and dst channels must be specified/omitted in pairs"
                      );
                    break;
                  } /*if*/
                get_dither_matrix(matrixobj, true, &dither);
                if (PyErr_Occurred())
                    break;
                some_src = 0;
                some_dst = 0;
                gotchan1 = get_channel(srcchan1obj, false, 0, &srcchan1);
                if (PyErr_Occurred())
                    break;
                if (gotchan1)
                  {
                    get_channel(dstchan1obj, true, 0, &dstchan1);
                    if (PyErr_Occurred())
                        break;
                    some_src = &srcchan1;
                    some_dst = &dstchan1;
                  } /*if*/
                gotchan2 = get_channel(srcchan2obj, false, some_src, &srcchan2);
                if (PyErr_Occurred())
                    break;
                if (gotchan2)
                  {
                    get_channel(dstchan2obj, true, some_dst, &dstchan2);
                    if (PyErr_Occurred())
                        break;
                    some_src = &srcchan2; /* in case not already set */
                    some_dst = &dstchan2;
                  } /*if*/
                gotchan3 = get_channel(srcchan3obj, false, some_src, &srcchan3);
                if (PyErr_Occurred())
                    break;
                if (gotchan3)
                  {
                    get_channel(dstchan3obj, true, some_dst, &dstchan3);
                    if (PyErr_Occurred())
                        break;
                    some_src = &srcchan3; /* in case not already set */
                    some_dst = &dstchan3;
                  } /*if*/
                gotchan4 = get_channel(srcchan4obj, false, some_src, &srcchan4);
                if (PyErr_Occurred())
                    break;
                if (gotchan4)
                  {
                    get_channel(dstchan4obj, true, some_dst, &dstchan4);
                    if (PyErr_Occurred())
                        break;
                    some_src = &srcchan4; /* in case not already set */
                    some_dst = &dstchan4;
                  } /*if*/
                if
                  (
                       (gotchan1 or gotchan2 or gotchan3 or gotchan4)
                   and
                        (some_src->width != some_dst->width or some_src->height != some_dst->height)
                  )
                  {
                    PyErr_SetString
                      (
                        PyExc_ValueError,
                        "src and dst channel pixmaps must have same dimensions"
                      );
                    break;
                  } /*if*/
              }
            while (false);
            Py_XDECREF(matrixobj);
            Py_XDECREF(srcchan1obj);
            Py_XDECREF(dstchan1obj);
            Py_XDECREF(srcchan2obj);
            Py_XDECREF(dstchan2obj);
            Py_XDECREF(srcchan3obj);
            Py_XDECREF(dstchan3obj);
            Py_XDECREF(srcchan4obj);
            Py_XDECREF(dstchan4obj);
          }
        if (PyErr_Occurred())
            break;
        if (not (gotchan1 or gotchan2 or gotchan3 or gotchan4) or to_depth == CPNT_BITS)
          {
          /* nothing to do */
            Py_INCREF(Py_None);
            result = Py_None; /* return success */
            break;
          } /*if*/
        /* assert some_src and some_dst both not null */
          { /* do the work */
            Py_BEGIN_ALLOW_THREADS
            uint drop_bits;
            const uint width = some_src->width;
            const uint height = some_src->height;
            const uint8_t * const src_pix_base = some_src->baseaddr;
            uint8_t * const dst_pix_base = some_dst->baseaddr;
            uint row, col;
            cpnt_t drop_mask, keep_mask;
            drop_bits = CPNT_BITS - to_depth;
            drop_mask = (1 << drop_bits) - 1;
            keep_mask = ~drop_mask;
            for (row = 0; row != height; ++row)
              {
              /* FIXME: currently always assuming pixel depth is 4 and bitwidth is 8! */
                const pix_type * const src_row_base = (const pix_type *)(src_pix_base + row * some_src->stride);
                pix_type * const dst_row_base = (pix_type *)(dst_pix_base + row * some_dst->stride);
                for (col = 0; col != width; ++col)
                  {
                    pix_type const srcpixel = src_row_base[col];
                    pix_type dstpixel = dst_row_base[col];
                  /* note: the following algorithm really only works correctly
                    when to_depth = 1. Otherwise it will produce colours that
                    are almost, but not quite, the same. This is because the
                    threshold comparison can only produce one of 2 values for
                    the bottommost rounded component bit. It needs to produce
                    (1 << to_depth) different values to map to all the rounded
                    component bits. This will require that number of threshold
                    coefficients at each pixel position. */
#define cond_do(gotchan, srcchan, dstchan) \
                    if (gotchan) \
                      { \
                        pix_type val = srcpixel >> srcchan.shiftoffset & CPNT_MAX; \
                        if (dither->order > 1) \
                          { \
                            uint const threshold = dither->coeffs[row % dither->order * dither->order + col % dither->order]; \
                            uint frac = val & (drop_mask << 1) + 1; \
                            frac = \
                                drop_bits > dither->bits ? \
                                    frac >> drop_bits - dither->bits \
                                : \
                                    frac << dither->bits - drop_bits; \
                            frac = frac > threshold ? CPNT_MAX : 0; \
                            val = val & keep_mask | frac & drop_mask; \
                          } \
                        else \
                          { \
                            val = (val & keep_mask) * CPNT_MAX / (CPNT_MAX - drop_mask); \
                              /* normalize truncated value to full brightness */ \
                          } /*if*/ \
                        dstpixel = dstpixel & ~(CPNT_MAX << dstchan.shiftoffset) | val << dstchan.shiftoffset; \
                      } /*if*/
                    cond_do(gotchan1, srcchan1, dstchan1)
                    cond_do(gotchan2, srcchan2, dstchan2)
                    cond_do(gotchan3, srcchan3, dstchan3)
                    cond_do(gotchan4, srcchan4, dstchan4)
#undef cond_do
                    dst_row_base[col] = dstpixel;
                  } /*for*/
              } /*for*/
            Py_END_ALLOW_THREADS
          }
      /* all successfully done */
        Py_INCREF(Py_None);
        result = Py_None;
      }
    while (false);
    free(dither);
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
    channel_t src, dst;
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
                get_channel(srcobj, true, 0, &src);
                if (PyErr_Occurred())
                    break;
                get_channel(dstobj, true, 0, &dst);
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
    {"ordered_dither", grainyx_ordered_dither, METH_VARARGS,
        "ordered_dither(matrix, depth, srcchan1 = None, dstchan1 = None, srcchan2 = None, dstchan2 = None, srcchan3 = None, dstchan3 = None, srcchan4 = None, dstchan4 = None)\n"
        "applies the specified matrix, which must be a DitherMatrix instance, to up to"
        " 4 pairs of source and destination Channel instances. All the source Channels must"
        " have the same baseaddr, width, height, depth and stride, and similarly all the"
        " destination Channels."
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
