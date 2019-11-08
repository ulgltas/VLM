# VLM
This **V**ortex **L**attice **M**ethod code is a 3D unsteady aerodynamics solver.
# Code compilation procedure
The code requires CMake and a C compiler. It works with `gcc`, `msvc` and `clang`. For the python wrapper, geoGen, SWIG, Python 2.7 and NumPy are required as well.

geoGen can be obtained from [this repository](https://github.com/acrovato/geoGen), currently it should be in a folder in the same folder as the VLM one.

## Linux-based systems
```bash
sudo apt-get install cmake
```
### Python wrapper
```bash
sudo apt-get install python2.7 python2.7-dev libpython2.7-dev
sudo apt-get install python-numpy python-scipy
sudo apt-get install swig
```

## OS X
You need both XCode and having the CMake Command Line Interface. The latter can be installed with Homebrew:
```bash
brew install cmake
```

### Python wrapper
Python is usually already present in OS X environments.
```bash
brew install swig
```

## Windows
Using Visual Studio 2017, code compiling **without the Python wrapper** is relatively straightforward. CMake is included by default. Make sure to disable the Python wrapper before running CMake, with option `-DPYTHON_WRAPPER=OFF`.
For both versions, tests can be run from the CMake menu.

### Python wrapper
You need:
1. A working Python installation with NumPy
2. To download [SWIG](http://www.swig.org/)
3. To provide CMake with the paths to these programs
4. To change the configuration to RELEASE

In order to provide the appropriate paths to Swig and Python you can modify the configurations, including in one of them the following line:
```json
"cmakeCommandArgs": "-DSWIG_EXECUTABLE=C:\\PATH\\TO\\swig.exe -DPYTHON_EXECUTABLE=C:\\PATH\\TO\\python.exe"
```

## Compilation

```bash
mkdir build && cd build
cmake ..
make
make install
```
# Running and testing
Once compiled, in the `build` directory, it is possible to run all the tests with `ctest`.
By default the program is installed in the `bin` subdirectory and can be called with
```bash
bin/VLM path/to/input/file
```

If using the Python wrapper
```bash
python run.py path/to/input/file
```
This will automatically create its own workspace folder. If you are not using the Python wrapper this command will fail.