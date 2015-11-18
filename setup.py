#+
# Distutils script to build and install the GrainyPy modules. Invoke from the command
# line in this directory as follows:
#
#     python3 setup.py install
#
# Written by Lawrence D'Oliveiro <ldo@geek-central.gen.nz>.
#-

import distutils.core

distutils.core.setup \
  (
    name = "GrainyPy",
    version = "0.5",
    description = "Make High-Quality Images Look Grainy",
    author = "Lawrence D’Oliveiro",
    author_email = "ldo@geek-central.gen.nz",
    url = "https://github.com/ldo/grainypy",
    py_modules = ["grainy"],
    ext_modules =
        [
            distutils.core.Extension
              (
                name = "grainyx",
                sources = ["grainyx.c"],
                libraries = ["m"],
                extra_compile_args = ["-Wno-parentheses", "-Wno-maybe-uninitialized"],
              ),
        ],
  )
