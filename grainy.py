#+
# grainy -- a Python module for making high-quality images look low-quality.
# This code offers functions for doing dithering to make images look like
# they are being shown on old, 1980s-vintage, colour-limited display hardware.
#-

import enum
import qahirah
from qahirah import \
    CAIRO
import grainyx

class Channel :
    "describes a pixel component. Attributes are:\n" \
    "    baseaddr    -- the base address of the pixmap\n" \
    "    width       -- the width in pixels\n" \
    "    height      -- the height in pixels\n" \
    "    depth       -- the number of bytes per component\n" \
    "    stride      -- the numberof bytes per row of pixels\n" \
    "    shiftoffset -- the offset in bits of the component from the bottom of the pixel\n" \
    "    bitwidth    -- the width in bits of the component"

    def __init__(self, baseaddr, width, height, depth, stride, shiftoffset, bitwidth) :
        self.baseaddr = baseaddr
        self.width = width
        self.height = height
        self.depth = depth
        self.stride = stride
        self.shiftoffset = shiftoffset
        self.bitwidth = bitwidth
    #end __init__

    def subrect(self, rect) :
        "returns a Channel descriptor for the specified rectangular portion of this Channel." \
        " rect should be a qahirah.Rect with integer bounds."
        rect.assert_isint()
        if (
                rect.left < 0 or rect.right > self.width
            or
                rect.top < 0 or rect.bottom > self.height
        ) :
            raise IndexError("rect exceeds available bounds")
        #end if
        return \
            Channel \
              (
                baseaddr = self.baseaddr + rect.top * self.stride + rect.left * self.depth,
                width = rect.width,
                height = rect.height,
                depth = self.depth,
                stride = self.stride,
                shiftoffset = self.shiftoffset,
                bitwidth = self.bitwidth
              )
    #end subrect

#end Channel

class CAIRO_PIX(enum.IntEnum) :
    "identifying components of an CAIRO.FORMAT_RGB24 or CAIRO.FORMAT_ARGB32 ImageSurface." \
    " Note this is endian-independent."
    A, R, G, B = 3, 2, 1, 0
#end CAIRO_PIX

class CairoARGBChannel(Channel) :
    "describes a pixel component for a CAIRO.FORMAT_RGB24 or CAIRO.FORMAT_ARGB32" \
    " ImageSurface. baseaddr, width, height and stride have the same meanings" \
    " as in the superclass; component is a CAIRO_PIX.xxx value indicating which" \
    " component to access."

    def __init__(self, baseaddr, width, height, stride, component) :
        super().__init__(baseaddr, width, height, 4, stride, int(component) * 8, 8)
    #end __init_

#end CairoARGBChannel

def cairo_component(pix, component) :
    "pix must be a qahirah.ImageSurface instance, and component must be a CAIRO_PIX.xxx" \
    " value. returns a Channel instance for acccessing the specified component of the" \
    " image pixmap."
    if not isinstance(pix, qahirah.ImageSurface) :
        raise TypeError("pix must be a qahirah.ImageSurface")
    #end if
    args = \
        {
            "baseaddr" : pix.data,
            "width" : pix.width,
            "height" : pix.height,
            "stride" : pix.stride,
        }
    construct = CairoARGBChannel
    if pix.format == CAIRO.FORMAT_ARGB32 :
        pass # fine
    elif pix.format == CAIRO.FORMAT_RGB24 :
        if component == CAIRO_PIX.A :
            raise ValueError("no alpha component for FORMAT_RGB24")
        #end if
    elif pix.format == CAIRO.FORMAT_A8 :
        if component != CAIRO_PIX.A :
            raise ValueError("alpha component only for FORMAT_A8")
        #end if
        component = None
        args["depth"] = 1
        args["shiftoffset"] = 0
        args["bitwidth"] = 8
        construct = Channel
    else :
        raise TypeError("unsupported format for pix")
    #end if
    if component != None :
        args["component"] = component
    #end if
    return \
        construct(**args)
#end cairo_component

class DitherMatrix :
    "base class for a general ordered-dither matrix of dimension order * order."
    " coeffs should be a tuple of the non-negative integer matrix elements in row" \
    " order, and bits should be the integer number of bits needed to hold one matrix" \
    " element value."

    __slots__ = ("bits", "coeffs", "order")

    def __init__(self, order, coeffs, bits) :
        assert isinstance(coeffs, tuple) and len(coeffs) == order * order
        self.order = order
        self.coeffs = coeffs
        self.bits = bits
    #end __init__

    def flip(self, x, y) :
        "returns a copy of this DitherMatrix with the coefficients flipped along one or both" \
        " axes. The x and y args are booleans indicating whether to flip along the corresponding" \
        " axis."
        return \
            DitherMatrix \
              (
                order = self.order,
                coeffs = tuple
                    (
                    self.coeffs[(row, self.order - row - 1)[y] * self.order + (col, self.order - col - 1)[x]]
                    for row in range(self.order)
                    for col in range(self.order)
                    ),
                bits = self.bits
              )
    #end flip

    def rotate(self, x, y) :
        "returns a copy of this DitherMatrix with the coefficients rotated by the specified" \
        " number of steps along each axis. The x and y args are integers indicating the" \
        " number of steps by which to rotate the coefficients along each axis; negative" \
        " for left/up, positive for right/down."
        return \
            DitherMatrix \
              (
                order = self.order,
                coeffs = tuple
                    (
                    self.coeffs[(row + y) % self.order * self.order + (col + x) % self.order]
                    for row in range(self.order)
                    for col in range(self.order)
                    ),
                bits = self.bits
              )
    #end rotate

#end DitherMatrix

class BayerMatrix(DitherMatrix) :
    "constructs a Bayer ordered-dither matrix. order must be a positive integer" \
    " which is a power of 2. Useful values are 2, 4, 8 or 16."

    def __init__(self, order) :
        if not isinstance(order, int) or order < 1 or order - 1 & order != 0 :
            raise ValueError("order must be a positive integer, being a power of 2")
        #end if
        if order > 2 :
            suborder = order // 2
            submat = BayerMatrix(suborder)
            mult = BayerMatrix(2).coeffs
            coeffs = tuple \
                (
                    4 * submat.coeffs[row % suborder * suborder + col % suborder]
                +
                    mult[row // suborder * 2 + col // suborder]
                for row in range(order)
                for col in range(order)
                )
            bits = submat.bits + 2
        elif order > 1 :
            coeffs = (3, 1, 0, 2)
            bits = 2
        else :
            coeffs = (1,)
            bits = 1
        #end if
        super().__init__(order, coeffs, bits)
    #end __init__

#end BayerMatrix

copy_channel = grainyx.copy_channel
ordered_dither = grainyx.ordered_dither

def copy_image_channel(src_img, src_component, dst_img, dst_component) :
    "copies a specified component from a source image to a component of a destination" \
    " image. src_image and dst_img must be qahirah.ImageSurface instances, while" \
    " src_component and dst_component must be CAIRO_PIX.xxx values selecting" \
    " the corresponding components."
    if (
            not isinstance(src_img, qahirah.ImageSurface)
        or
            not isinstance(dst_img, qahirah.ImageSurface)
    ) :
        raise TypeError("src_img and dst_img must be qahirah.ImageSurface instances")
    #end if
    src_img.flush()
    dst_img.flush()
    copy_channel \
      (
        cairo_component(src_img, src_component), # src
        cairo_component(dst_img, dst_component), # dst
      )
    dst_img.mark_dirty()
#end copy_image_channel

def ordered_dither_image(src_img, dst_img, depth, matrix, src_bounds, dst_bounds, do_a, do_r, do_g, do_b) :
    "dithers src_image into the corresponding components of dst_img using the specified" \
    " DitherMatrix, according to the booleans do_a, do_r, do_g and do_b. src_img and dst_img" \
    " may be the same image."
    if (
            not isinstance(src_img, qahirah.ImageSurface)
        or
            not isinstance(dst_img, qahirah.ImageSurface)
    ) :
        raise TypeError("src_img and dst_img must be qahirah.ImageSurface instances")
    #end if
    src_img.flush()
    dst_img.flush()
    srcchan = [None] * 4
    dstchan = [None] * 4
    nr_components = 0
    for \
        doit, component \
    in \
        (
            (do_a, CAIRO_PIX.A),
            (do_r, CAIRO_PIX.R),
            (do_g, CAIRO_PIX.G),
            (do_b, CAIRO_PIX.B),
        ) \
    :
        if doit :
            srcchan[nr_components] = cairo_component(src_img, component)
            if src_bounds != None :
                srcchan[nr_components] = srcchan[nr_components].subrect(src_bounds)
            #end if
            dstchan[nr_components] = cairo_component(dst_img, component)
            if dst_bounds != None :
                dstchan[nr_components] = dstchan[nr_components].subrect(dst_bounds)
            #end if
            nr_components += 1
        #end if
    #end for
    ordered_dither \
      (
        matrix,
        depth,
        srcchan[0], dstchan[0],
        srcchan[1], dstchan[1],
        srcchan[2], dstchan[2],
        srcchan[3], dstchan[3],
      )
    dst_img.mark_dirty()
#end ordered_dither_image
