# VLM
This **V**ortex **L**attice **M**ethod code is a 3D unsteady aerodynamics solver.
# Code compilation procedure
The code requires CMake and a C compiler. It works both with `gcc` and `clang`.

## Linux-based systems
```bash
sudo apt-get install cmake
```
## OS X
You need both XCode and to have the CMake Command Line Interface. The latter can be installed with Homebrew:
```bash
brew install cmake
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