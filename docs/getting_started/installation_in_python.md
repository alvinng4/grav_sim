For basic install, follow the quick guide at [Getting Started](getting_started.md).
Here we provide a more detailed installation guide for Python.

Below are the dependencies for the package:

- C compiler, preferably GCC or Clang
- CMake
- Python version >= 3.9

Optional dependencies:

- HDF5 (For data storage. Optional but recommended)
- FFTW3 (Only required for cosmological simulations)
- OpenMP

Python packages: (Check `requirements.txt` for the latest version)

- NumPy
- Matplotlib
- h5py

There are multiple ways to install the package, depending on your use case.

/// tab | Installing with pip and local compilation (MacOS / Linux)
If you want to compile the package locally, you could run
```
// Choose the options you want
export CMAKE_ARGS="-DUSE_HDF5=ON -DUSE_OPENMP=ON -DUSE_FFTW3=ON -DCMAKE_C_COMPILER=gcc .."

pip install grav_sim --no-binary grav_sim --no-cache-dir
```

  * `-DUSE_OPENMP=ON -DUSE_FFTW3=ON -DUSE_HDF5=ON`: optional flags for the dependencies.
  * `-DCMAKE_C_COMPILER=gcc`: optional flag to specify the C compiler.

To check whether it is successfully installed, run
```
python -m grav_sim [--path=c_lib_path]
```
where `--path=c_lib_path` is optional argument to specify the path
to the compiled C library. You should see the compilation information and the path to the
compiled library.
```
-----------------------------------------------------------------
                                              __                   
    __   _ __    __     __  __           ____/\_\    ___ ___       
  /'_ `\/\`'__\/'__`\  /\ \/\ \         /',__\/\ \ /' __` __`\     
 /\ \L\ \ \ \//\ \L\.\_\ \ \_/ |       /\__, `\ \ \/\ \/\ \/\ \    
 \ \____ \ \_\\ \__/.\_\\ \___/        \/\____/\ \_\ \_\ \_\ \_\   
  \/___L\ \/_/ \/__/\/_/ \/__/   _______\/___/  \/_/\/_/\/_/\/_/   
    /\____/                     /\______\                          
    \_/__/                      \/______/                          


grav_sim version 0.0.4

Operating System: MacOS
Compilation Info:
Compiled with OpenMP: true
Compiled with HDF5: true
    Version: 1.14.6
Compiled with FFTW3: false

Build time: Apr 16 2025 13:12:32
Compiler: GCC (version: 14)
-----------------------------------------------------------------
C library location: /Library/Frameworks/Python.framework/Versions/3.11/lib/python3.11/site-packages/libgrav_sim.dylib
```
///

/// tab | Compiling the C library directly from source

This way should work on all platforms as long as you could get the
C library compiled on your system.

1. Clone the repository to your local machine and navigate to the directory:
```
git clone https://github.com/alvinng4/grav_sim
cd grav_sim
```
2. Create a build directory inside the repository:
```
mkdir build
cd build
```
3. Compile the library with CMake
```
cmake [-DUSE_OPENMP=ON -DUSE_FFTW3=ON -DUSE_HDF5=ON] [-DCMAKE_C_COMPILER=gcc] ..
cmake --build .
```

    * `-DUSE_OPENMP=ON -DUSE_FFTW3=ON -DUSE_HDF5=ON`: optional flags for the dependencies.
    * `-DCMAKE_C_COMPILER=gcc`: optional flag to specify the C compiler.

4. Check the compilation. You should see one of the following files in the `build` directory.
```
libgrav_sim.dylib
libgrav_sim.so
libgrav_sim.dll
```
This is the compiled C library. You will need to keep track of the path to this file.

5. Now you can install the Python wrapper either from PyPI or from source.

    /// tab | MacOS / Linux
      - To install from PyPI, run
      ```
      SKIP_CMAKE=1 pip install grav_sim
      ```
      - To install from source, navigate to the parent directory and run
      ```
      SKIP_CMAKE=1 pip install .
      ```
    ///

    /// tab | Windows
      - To install from PyPI, run
      ```
      set SKIP_CMAKE=1
      pip install grav_sim
      ```
      - To install from source, navigate to the parent directory and run
      ```
      set SKIP_CMAKE=1
      pip install .
      ```
    ///

6. To check whether it is successfully installed, run
  ```
  python -m grav_sim --path=/path/to/compiled/c_lib
  ```
  where `--path=/path/to/compiled/c_lib` is the path to the compiled C library.
  You should see the compilation information and the path to the compiled library.

  ```
  -----------------------------------------------------------------
                                              __                   
    __   _ __    __     __  __           ____/\_\    ___ ___       
  /'_ `\/\`'__\/'__`\  /\ \/\ \         /',__\/\ \ /' __` __`\     
  /\ \L\ \ \ \//\ \L\.\_\ \ \_/ |       /\__, `\ \ \/\ \/\ \/\ \    
  \ \____ \ \_\\ \__/.\_\\ \___/        \/\____/\ \_\ \_\ \_\ \_\   
    \/___L\ \/_/ \/__/\/_/ \/__/   _______\/___/  \/_/\/_/\/_/\/_/   
      /\____/                     /\______\                          
      \_/__/                      \/______/                          


  grav_sim version 0.0.5

  Operating System: MacOS
  Compilation Info:
    Compiled with OpenMP: false
    Compiled with HDF5: false
    Compiled with FFTW3: false

  Build time: Apr 28 2025 17:06:18
  Compiler: Clang (version: 15)
  -----------------------------------------------------------------
  C library location: /Users/alvinng/Desktop/cuhk/gravity-simulator/Final Project/grav_sim/build/lib.macosx-10.9-universal2-cpython-311/libgrav_sim.dylib
  ```
///