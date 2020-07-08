from distutils.core import setup, Extension

def main():
    setup(name="wavefunc",
          version="1.0.0",
          description="Python interface for routines which solve for properties of wavefunctions",
          author="harry knight",
          author_email="harry.knight@hotmail.co.uk",
          ext_modules=[Extension("wavefunc", ["ode.c", "brent.c", "pwavefuncv2.c"],
          include_dirs=['/home/harry/OneDrivePersonal/Interactive_QHE', '/home/harry/OneDrivePersonal/Interactive_QHE/.qhevenv/lib/python3.7/site-packages/numpy/core/include/numpy'])])

if __name__ == "__main__":
    main()
