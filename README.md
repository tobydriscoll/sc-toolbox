sc-toolbox
==========

Schwarz-Christoffel Toolbox for conformal mapping in MATLAB

The SC Toolbox contains numerical routines and graphical interfaces to work with Schwarz-Christoffel conformal maps--those to regions bounded by polgons in the complex plane. Many map variations are present. The software has no requirements other than core MATLAB.

You might prefer to view the [page at the File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/1316-schwarz-christoffel-toolbox), where you can try the package out online without downloading and installing it.

For more details on the maps, see _Schwarz Christoffel Mapping_, by Driscoll and Trefethen. For a user's guide, visit https://tobydriscoll.net/project/sc-toolbox/.

## Ports

- **C++** — a function-to-function C++17 port of the numerical core lives in [`cpp/`](cpp), validated against golden values generated from the MATLAB reference. See [`CPP_PLAN.md`](CPP_PLAN.md).
- **Python** — NumPy-native bindings over the C++ port (built with nanobind) live in [`python/`](python). See [`python/README.md`](python/README.md) for install instructions, examples, and a rendered gallery.
