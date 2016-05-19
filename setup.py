#!/usr/bin/env python

from __future__ import division, absolute_import, print_function

from numpy.distutils.core import Extension

ext = Extension('libslater.ext',
                 sources = ['src/libslater.f90'])

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(name = 'libslater',
          description       = "LibSlater",
          author            = "Steven Vandenbrande",
          author_email      = "Steven.Vandenbrande@UGent.Be",
          package_dir={'libslater': 'src'},
          packages=['libslater'],
          ext_modules = [ext]
          )
