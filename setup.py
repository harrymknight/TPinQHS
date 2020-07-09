from setuptools import setup, Extension
import numpy as np

wavefunc_ext = Extension(
    name="wavefunc",
    sources=["ode.c", "brent.c", "pywavefunc.c"],
    include_dirs=[np.get_include()] 
)

def main():
    setup(name="wavefunc",
          version="1.0.0",
          description="Python interface for routines which solve for properties of wavefunctions",
          author="harry knight",
          author_email="harry.knight@hotmail.co.uk",
          ext_modules=[wavefunc_ext])

if __name__ == "__main__":
    main()
