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
