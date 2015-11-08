#+
# Higher-level wrapper around lower-level extension module code.
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

class BayerMatrix :
    "constructs a Bayer ordered-dither matrix. order must be a positive integer" \
    " which is a power of 2. Useful values are 2, 4, 8 or 16. coeffs will be a tuple" \
    " of all the constructed matrix elements, and bits will be the integer number of" \
    " bits needed to hold one matrix element value."

    __slots__ = ("bits", "coeffs", "order")

    def __init__(self, order) :
        # setting order = None is for internal use only
        if order != None :
            if not isinstance(order, int) or order < 1 or order - 1 & order != 0 :
                raise ValueError("order must be a positive integer, being a power of 2")
            #end if
            if order > 2 :
                suborder = order // 2
                submat = BayerMatrix(suborder)
                mult = BayerMatrix(2).coeffs
                self.coeffs = tuple \
                    (
                        4 * submat.coeffs[row % suborder * suborder + col % suborder]
                    +
                        mult[row // suborder * 2 + col // suborder]
                    for row in range(order)
                    for col in range(order)
                    )
                self.bits = submat.bits + 2
            elif order > 1 :
                self.coeffs = (3, 1, 0, 2)
                self.bits = 2
            else :
                self.coeffs = (1,)
                self.bits = 1
            #end if
            self.order = order
        #end if
    #end __init__

    def flip(self, x, y) :
        "returns a copy of this BayerMatrix with the coefficients flipped along one or both" \
        " axes. The x and y args are booleans indicating whether to flip along the corresponding" \
        " axis."
        result = BayerMatrix(None)
        result.coeffs = tuple \
            (
            self.coeffs[(row, self.order - row - 1)[y] * self.order + (col, self.order - col - 1)[x]]
            for row in range(self.order)
            for col in range(self.order)
            )
        result.bits = self.bits
        result.order = self.order
        return \
            result
    #end flip

    def rotate(self, x, y) :
        "returns a copy of this BayerMatrix with the coefficients rotated by the specified" \
        " number of steps along each axis. The x and y args are integers indicating the" \
        " number of steps by which to rotate the coefficients along each axis; negative" \
        " for left/up, positive for right/down."
        result = BayerMatrix(None)
        result.coeffs = tuple \
            (
            self.coeffs[(row + y) % self.order * self.order + (col + x) % self.order]
            for row in range(self.order)
            for col in range(self.order)
            )
        result.bits = self.bits
        result.order = self.order
        return \
            result
    #end rotate

#end BayerMatrix

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

copy_channel = grainyx.copy_channel

def copy_image_channel(src_img, src_component, dst_img, dst_component) :
    "copies a specified component from a source image to a component of a destination" \
    " image. src_image and dst_img must be qahirah.ImageSurface instances, while" \
    " src_component and dst_component must be CAIRO_PIX.xxx values selecting" \
    " the corresponding components."
    if not isinstance(src_img, qahirah.ImageSurface) or not isinstance(dst_img, qahirah.ImageSurface) :
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
