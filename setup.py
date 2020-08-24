from setuptools import setup, Extension
import numpy as np

wavefunc_ext = Extension(
    name="wavefunc",
    sources=["ode.c", "brent.c", "capsule.c", "pywavefunc.c"],
    include_dirs=[np.get_include()])

def main():
    setup(name="wavefunc",
          version="1.0.0",
          description="Python interface for routines which solve for properties of wavefunctions",
          author="harry knight",
          author_email="harry.knight@hotmail.co.uk",
          ext_modules=[wavefunc_ext],
          install_requires=[
              'dash>=1.13.4',
              'dash-bootstrap-components>=0.10.3',
              'dash-defer-js-import',
              'pandas>=1.0.5',
              'numpy>=1.19.0'
          ])

if __name__ == "__main__":
    main()
