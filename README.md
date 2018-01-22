# R-tree (CUDA)
ala_tree (C/C++ + CUDA)

Implementation of [Rtree](https://en.wikipedia.org/wiki/R-tree). Works under CPU and NVidia GPU (CUDA). Depending on the load, you can do calculations on the CPU or GPU.
The code is cross-platform for Windows and Linux.
To build under Windows, there is a solution for Visual Studio 2015+, for Linux - Makefile. In order for everything to work, you need a cuda-compatible device and a CUDA Toolkit.

At the moment, preprocessor definitions can be used in the assembly:
- USE_CUDA - compile the code using CUDA;
- CALC_POINT - compile the code that searches for the nearest object (only for CPU). Mutually excludes CALC_CIRCLE;
- CALC_CIRCLE - compile code that can return objects not only from a rectangle, but also from a circle (only for CPU). Mutually excludes CALC_POINT.

Preprocessor definitions (#define) used in the code:
- MAX_ITEMS_IN_NODE (first.h) - the maximum number of objects (leaves) in node;
- MAX_NODES (first.h) - the maximum number of nodes at any level;
- MAX_RESULTS (search2.cu) - the maximum number of results when calculating on the GPU;
- PACK_RESULTS (search2.cu) - flag that allows you to remove duplicate results when calculating on the GPU.

The code is embedded in any application in c / c ++. To demonstrate the work there is a file test.cpp, in which an example of usage is given.
