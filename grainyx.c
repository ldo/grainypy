/*
    Lower-level code for grainy module. This C code implements the most
    CPU-intensive parts of the algorithms. It is written to make as few
    assumptions about data structures as possible, leaving it to the
    Python layer to make things easier for the caller.

    Copyright 2015 Lawrence D'Oliveiro <ldo@geek-central.gen.nz>.
    Licensed under the GNU Lesser General Public License v2.1 or later.
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
  { /* identifies a pixel component within a pixmap */
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
    uint coeffs[/*order * order*/];
  } ordered_dither_matrix_t;

static bool get_channel
  (
    PyObject * desc,
    bool required, /* indicates desc cannot be None for no channel */
    const channel_t * same_pix_as, /* if non-null, then result must have same baseaddr, width, height, stride and depth as this */
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

static void get_ordered_dither_matrix
  (
    PyObject * obj,
    bool required, /* indicates obj cannot be None for no matrix */
    ordered_dither_matrix_t ** mat
  )
  /* extracts the fields from a grainy.DitherMatrix instance and allocates and fills in
    a new ordered_dither_matrix_t structure with the results. */
  {
    *mat = 0; /* to begin with */
    if (required or obj != 0 and obj != Py_None)
      {
        ordered_dither_matrix_t * result = 0;
        PyObject * field = 0;
        uint order;
        do /*once*/
          {
            field = PyObject_GetAttrString(obj, "order");
            if (PyErr_Occurred())
                break;
            order = PyLong_AsUnsignedLong(field);
            if (PyErr_Occurred())
                break;
            result = malloc(sizeof(ordered_dither_matrix_t) + sizeof(uint) * order * order);
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
            for (uint row = 0;;)
              {
                if (row == order)
                    break;
                for (uint col = 0;;)
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
  } /*get_ordered_dither_matrix*/

/*
    User-visible stuff

    Note on PyArg_ParseTuple calls: this may have partially succeeded even if
    it returned false, with some arguments returned but not others. To deal
    with this, I need to ensure I correctly allocate references to all returned
    objects, so they get correctly released at the end. You’ll note that this is
    done before checking whether the call succeeded or not.
*/

static PyObject * grainyx_ordered_dither
  (
    PyObject * self,
    PyObject * args
  )
  {
    PyObject * result = 0;
    ordered_dither_matrix_t * dither = 0;
    uint to_depth;
    channel_t srcchan[4], dstchan[4];
    bool gotchan[4];
    channel_t * some_src;
    channel_t * some_dst;
    do /*once*/
      {
          { /* get args */
            PyObject * matrixobj = 0;
            PyObject * srcchanobj[4] = {0, 0, 0, 0};
            PyObject * dstchanobj[4] = {0, 0, 0, 0};
            do /*once*/
              {
                const bool success = PyArg_ParseTuple
                  (
                    args,
                    "OIOO|OOOOOO",
                    &matrixobj,
                    &to_depth,
                    srcchanobj + 0, dstchanobj + 0,
                    srcchanobj + 1, dstchanobj + 1,
                    srcchanobj + 2, dstchanobj + 2,
                    srcchanobj + 3, dstchanobj + 3
                  );
                Py_XINCREF(matrixobj);
                Py_XINCREF(srcchanobj[0]);
                Py_XINCREF(dstchanobj[0]);
                for (uint chan = 1; chan != 4; ++chan)
                  {
                    if (srcchanobj[chan] == 0)
                      {
                        srcchanobj[chan] = Py_None;
                      } /*if*/
                    if (dstchanobj[chan] == 0)
                      {
                        dstchanobj[chan] = Py_None;
                      } /*if*/
                    Py_INCREF(srcchanobj[chan]);
                    Py_INCREF(dstchanobj[chan]);
                  } /*if*/
                if (not success)
                    break;
                if
                  (
                        (srcchanobj[0] != Py_None) != (dstchanobj[0] != Py_None)
                    or
                        (srcchanobj[1] != Py_None) != (dstchanobj[1] != Py_None)
                    or
                        (srcchanobj[2] != Py_None) != (dstchanobj[2] != Py_None)
                    or
                        (srcchanobj[3] != Py_None) != (dstchanobj[3] != Py_None)
                  )
                  {
                    PyErr_SetString
                      (
                        PyExc_ValueError,
                        "src and dst channels must be specified/omitted in pairs"
                      );
                    break;
                  } /*if*/
                get_ordered_dither_matrix(matrixobj, true, &dither);
                if (PyErr_Occurred())
                    break;
                some_src = 0;
                some_dst = 0;
                for (uint chan = 0;;)
                  {
                    if (chan == 4)
                        break;
                    gotchan[chan] = get_channel(srcchanobj[chan], false, some_src, srcchan + chan);
                    if (PyErr_Occurred())
                        break;
                    if (gotchan[chan])
                      {
                        get_channel(dstchanobj[chan], true, some_dst, dstchan + chan);
                        if (PyErr_Occurred())
                            break;
                        some_src = srcchan + chan; /* in case not already set */
                        some_dst = dstchan + chan;
                      } /*if*/
                    ++chan;
                  } /*for*/
                if (PyErr_Occurred())
                    break;
                if
                  (
                       (gotchan[0] or gotchan[1] or gotchan[2] or gotchan[3])
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
            for (uint chan = 0; chan != 4; ++chan)
              {
                Py_XDECREF(srcchanobj[chan]);
                Py_XDECREF(dstchanobj[chan]);
              } /*for*/
          }
        if (PyErr_Occurred())
            break;
        if (not (gotchan[0] or gotchan[1] or gotchan[2] or gotchan[3]) or to_depth == CPNT_BITS)
          {
          /* nothing to do */
            Py_INCREF(Py_None);
            result = Py_None; /* return success */
            break;
          } /*if*/
        /* assert some_src and some_dst both not null */
          { /* do the work */
            Py_BEGIN_ALLOW_THREADS
            const uint drop_bits = CPNT_BITS - to_depth;
            const cpnt_t drop_mask = (1 << drop_bits) - 1;
            const cpnt_t keep_mask = ~drop_mask;
            const uint width = some_src->width;
            const uint height = some_src->height;
            const uint8_t * const src_pix_base = some_src->baseaddr;
            uint8_t * const dst_pix_base = some_dst->baseaddr;
            uint row, col;
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
                    for (uint chan = 0; chan != 4; ++chan)
                      {
                        if (gotchan[chan])
                          {
                            pix_type val = srcpixel >> srcchan[chan].shiftoffset & CPNT_MAX;
                            if (dither->order > 1)
                              {
                                uint threshold = dither->coeffs
                                  [
                                    row % dither->order * dither->order + col % dither->order
                                  ];
                                uint frac = val & drop_mask;
                                if (drop_bits > dither->bits)
                                  {
                                    threshold <<= drop_bits - dither->bits;
                                  }
                                else if (drop_bits < dither->bits)
                                  {
                                    frac <<= dither->bits - drop_bits;
                                  } /*if*/
                                frac = frac > threshold ? CPNT_MAX : 0;
                                val = val & keep_mask | frac & drop_mask;
                              }
                            else
                              {
                                val = (val & keep_mask) * CPNT_MAX / (CPNT_MAX - drop_mask);
                                  /* normalize truncated value to full brightness */
                              } /*if*/
                            dstpixel =
                                    dstpixel & ~(CPNT_MAX << dstchan[chan].shiftoffset)
                                |
                                    val << dstchan[chan].shiftoffset;
                          } /*if*/
                      } /*for*/
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
                const bool success = PyArg_ParseTuple(args, "OO", &srcobj, &dstobj);
                Py_XINCREF(srcobj);
                Py_XINCREF(dstobj);
                if (not success)
                    break;
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

static PyObject * grainyx_channel_op
  (
    PyObject * self,
    PyObject * args
  )
  {
    PyObject * result = 0;
    cpnt_t * table = 0;
    channel_t srcl[4], srcr[4], dst[4];
    bool gotchan[4];
    const channel_t * some_srcl;
    const channel_t * some_srcr;
    const channel_t * some_dst;
    do /*once*/
      {
          { /* get args */
            PyObject * tableobj = 0;
            PyObject * srclobj[4] = {0, 0, 0, 0};
            PyObject * srcrobj[4] = {0, 0, 0, 0};
            PyObject * dstobj[4] = {0, 0, 0, 0};
            uint table_size;
            do /*once*/
              {
                const bool success = PyArg_ParseTuple
                  (
                    args, "O|OOOOOOOOOOOO",
                    &tableobj,
                    srclobj + 0, srcrobj + 0, dstobj + 0,
                    srclobj + 1, srcrobj + 1, dstobj + 1,
                    srclobj + 2, srcrobj + 2, dstobj + 2,
                    srclobj + 3, srcrobj + 3, dstobj + 3
                  );
                Py_XINCREF(tableobj);
              /* replace omitted arguments with None to simplify checks */
                for (uint chan = 0; chan != 4; ++chan)
                  {
                    if (srclobj[chan] == 0)
                      {
                        srclobj[chan] = Py_None;
                      } /*if*/
                    if (srcrobj[chan] == 0)
                      {
                        srcrobj[chan] = Py_None;
                      } /*if*/
                    if (dstobj[chan] == 0)
                      {
                        dstobj[chan] = Py_None;
                      } /*if*/
                    Py_INCREF(srclobj[chan]);
                    Py_INCREF(srcrobj[chan]);
                    Py_INCREF(dstobj[chan]);
                  } /*for*/
                if (not success)
                    break;
                for (uint chan = 0;;)
                  {
                    if (chan == 4)
                        break;
                    if
                      (
                            (srclobj[chan] != Py_None) != (dstobj[chan] != Py_None)
                        or
                            (srcrobj[chan] != Py_None) != (dstobj[chan] != Py_None)
                      )
                      {
                        PyErr_SetString
                          (
                            PyExc_ValueError,
                            "srcl, srcr and dst channels must be specified/omitted together"
                          );
                        break;
                      } /*if*/
                    ++chan;
                  } /*for*/
                if (PyErr_Occurred())
                    break;
                some_srcl = 0;
                some_srcr = 0;
                some_dst = 0;
                for (uint chan = 0;;)
                  {
                    if (chan == 4)
                        break;
                    gotchan[chan] = get_channel(srclobj[chan], false, some_srcl, srcl + chan);
                    if (PyErr_Occurred())
                        break;
                    if (gotchan[chan])
                      {
                        get_channel(srcrobj[chan], true, some_srcr, srcr + chan);
                        if (PyErr_Occurred())
                            break;
                        get_channel(dstobj[chan], true, some_dst, dst + chan);
                        if (PyErr_Occurred())
                            break;
                        some_srcl = srcl + chan; /* in case not already set */
                        some_srcr = srcr + chan;
                        some_dst = dst + chan;
                      } /*if*/
                    ++chan;
                  } /*for*/
                if (PyErr_Occurred())
                    break;
                if (gotchan[0] or gotchan[1] or gotchan[2] or gotchan[3])
                  {
                    if
                      (
                            some_srcl->bitwidth != some_srcr->bitwidth
                        or
                            some_srcr->bitwidth != some_dst->bitwidth
                      )
                      {
                        PyErr_SetString
                          (
                            PyExc_ValueError,
                            "srcl, srcr and dst must all have same bitwidth"
                          );
                        break;
                      } /*if*/
                    if (some_srcl->bitwidth > CPNT_BITS)
                      {
                        PyErr_Format(PyExc_ValueError, "bitwidth cannot exceed %d for now", CPNT_BITS);
                        break;
                      } /*if*/
                    if
                      (
                            some_srcl->depth != 1 and some_srcl->depth != 4
                        or
                            some_srcr->depth != 1 and some_srcr->depth != 4
                        or
                            some_dst->depth != 1 and some_dst->depth != 4
                      )
                      {
                        PyErr_SetString
                          (
                            PyExc_ValueError,
                            "only pixel depths of 1 and 4 bytes currently supported"
                          );
                        break;
                      } /*if*/
                    table_size = 1 << some_srcl->bitwidth * 2;
                    table = calloc(table_size, sizeof(cpnt_t));
                    if (table == 0)
                      {
                        PyErr_NoMemory();
                        break;
                      } /*if*/
                    for (uint index = 0;;)
                      {
                        if (index == table_size)
                            break;
                        PyObject * const elt = PyTuple_GetItem(tableobj, index);
                          /* borrowed reference, valid as long as tableobj is valid */
                        if (PyErr_Occurred())
                            break;
                        table[index] = PyLong_AsUnsignedLong(elt);
                        if (PyErr_Occurred())
                            break;
                        ++index;
                      } /*for*/
                    if (PyErr_Occurred())
                        break;
                  } /*if*/
              }
            while (false);
            Py_XDECREF(tableobj);
            for (uint chan = 0; chan != 4; ++chan)
              {
                Py_XDECREF(srclobj[chan]);
                Py_XDECREF(srcrobj[chan]);
                Py_XDECREF(dstobj[chan]);
              } /*for*/
          }
        if (PyErr_Occurred())
            break;
        if (not (gotchan[0] or gotchan[1] or gotchan[2] or gotchan[3]))
          {
          /* nothing to do */
            Py_INCREF(Py_None);
            result = Py_None; /* return success */
            break;
          } /*if*/
        Py_BEGIN_ALLOW_THREADS
          { /* do the work */
            for (uint row = 0; row != some_srcl->height; ++row)
              {
                const cpnt_t * const srclrow = some_srcl->baseaddr + row * some_srcl->stride;
                const cpnt_t * const srcrrow = some_srcr->baseaddr + row * some_srcr->stride;
                cpnt_t * const dstrow = some_dst->baseaddr + row * some_dst->stride;
                for (uint col = 0; col != some_srcl->width; ++col)
                  {
                    for (uint chan = 0; chan != 4; ++chan)
                      {
                        if (gotchan[chan])
                          {
                            const uint srclmask = (1 << srcl[chan].bitwidth) - 1; 
                            const uint srcrmask = (1 << srcr[chan].bitwidth) - 1;
                            const uint dstmask =
                                (1 << dst[chan].bitwidth) - 1 << dst[chan].shiftoffset;
                            const uint srclpix =
                                srcl[chan].depth == 4 ?
                                    ((uint32_t *)srclrow)[col]
                                : /* srcl.depth = 1 */
                                    srclrow[col];
                            const uint srcrpix =
                                srcr[chan].depth == 4 ?
                                    ((uint32_t *)srcrrow)[col]
                                : /* srcr[chan].depth = 1 */
                                    srcrrow[col];
                            uint dstpix =
                                dst[chan].depth == 4 ?
                                    ((uint32_t *)dstrow)[col]
                                : /* dst[chan].depth = 1 */
                                    dstrow[col];
                            dstpix =
                                    dstpix & ~dstmask
                                |
                                            table
                                                [
                                                        (
                                                            srclpix >> srcl[chan].shiftoffset
                                                        &
                                                            srclmask
                                                        )
                                                    <<
                                                        srcl[chan].bitwidth
                                                |
                                                    srcrpix >> srcr[chan].shiftoffset & srcrmask
                                                ]
                                        <<
                                            dst[chan].shiftoffset
                                    &
                                        dstmask;
                            if (dst[chan].depth == 4)
                              {
                                ((uint32_t *)dstrow)[col] = dstpix;
                              }
                           else /* dst[chan].depth = 1 */
                             {
                                dstrow[col] = dstpix;
                             } /*if*/
                          } /*if*/
                      } /*for*/
                  } /*for*/
              } /*for*/
          }
        Py_END_ALLOW_THREADS
      /* all successfully done */
        Py_INCREF(Py_None);
        result = Py_None;
      }
    while (false);
    free(table);
    return
        result;
  } /*grainyx_channel_op*/

static PyMethodDef grainyx_methods[] =
  {
    {"ordered_dither", grainyx_ordered_dither, METH_VARARGS,
        "ordered_dither(matrix, depth, srcchan1 = None, dstchan1 = None, srcchan2 = None, dstchan2 = None, srcchan3 = None, dstchan3 = None, srcchan4 = None, dstchan4 = None)\n"
        "applies the specified matrix, which must be a DitherMatrix instance, to up to"
        " 4 pairs of source and destination Channel instances. All the source Channels must"
        " come from the same pixmap, and similarly all the destination Channels."
    },
    {"copy_channel", grainyx_copy_channel, METH_VARARGS,
        "copy_channel(src_channel, dst_channel)\n"
        "copies the contents of one Channel instance to another. You may copy"
        " between two channels of the same pixmap, but the pixmaps should not otherwise"
        " overlap."
    },
    {"channel_op", grainyx_channel_op, METH_VARARGS,
        "channel_op(table, srcl1 = None, srcr1 = None, dst1 = None, srcl2 = None,"
            " srcr2 = None, dst2 = None, srcl3 = None, srcr3 = None, dst3 = None,"
            " srcl4 = None, srcr4 = None, dst4 = None)\n"
        "performs a general functional operation on component values from srcln and srcrn,"
        " putting the result into dstn, for up to 4 sets of sources and destinations. All the"
        " srcln channels must come the same pixmap; similarly all the srcr from the same pixmap,"
        " and all the dst channels from the same pixmap. All pixmaps must have the same width"
        " and height, and all channels must have the same bitwidth. table is a tuple of"
        " integers, being the new destination component values corresponding to each possible"
        " pair of source component values. It is indexed by"
        " srcl_component << bitwidth | srcr_component, so its size must be 1 << 2 * bitwidth."
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
