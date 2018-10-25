Partitions of Unity Extension (PUX)
Extends function, given on 2D domain, to the plane. Extension has compact support and high regularity. Based on the paper https://arxiv.org/abs/1712.08461.

Tests are located in the test directory. DEMOextension shows how to use PUX. In the DEMOextension.m file the user can choose between different cases: different geometries and functions to extend. The demos example2 and example3 are example 2 and example 3 from paper https://arxiv.org/abs/1712.08461. They show how to incorporate PUX in a Poisson solver. It is recommended to save BIMstruct, xe and idxbO if the domain is to be reused.

Prerequisites
RBF-QR for the 2D case: http://www.it.uu.se/research/scientific_computing/software/rbf_qr
Non-uniform Fast Fourier Transform for 2D (NUFFT): https://cims.nyu.edu/cmcl/nufft/nufft.html

Store both in src directory.
The  NUFFT is called from the function callNUFFT.m and RBF-QR is called from setupPUX.

Tests:
The tests and demos are run by executing DEMO_extension.m, exmaple2.m and example3.m.

Remark: This code is intended to show the structure of PUX and the Poisson solver from https://arxiv.org/abs/1712.08461. Optimization has not been the main focus.


Authors

Fredrik Fryklund
Erik Lehto.

License

This project is licensed under the MIT License - see the LICENSE.md file for details
