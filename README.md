## License

BSD-type license, see legal.txt and license.txt.


## Compilation

Cadmus requires MPI and DIY2.
DIY2 is used as a header-only library.


```
git clone git@gitlab.kitware.com:diatomic/diy.git
git@github.com:anigmetov/cadmus.git
cd cadmus
mkdir build
cd build
cmake .. -DDIY_INCLUDE_DIR=../diy/include
```

## Usage

Cadmus reads data from .npy files and writes the computed persistence diagrams
to disk.  Since it is now templated on the cube dimension, there are 3 executables,
`cadmus_1`, `cadmus_2` and `cadmus_3` in the `build/src` directory.

To run it, use
``mpirun -n 64 ./cadmus_3 [-n] [-c] [-f] input_3D_array.npy output_dgm_filename``
or an equivalent `srun` command.

Flags: 
- `-n` for `negate` (if you want upper-star persistence)
- `-c` for `clearing` (clearing optimization is highly recommended)
- `-f` for `Freudenthal` (use Freudenthal triangulation instead of cubical
        complexes)

