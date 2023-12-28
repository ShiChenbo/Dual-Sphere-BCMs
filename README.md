# Characteristic Modes amidst Arbitrary Background using Scattering Matrices

We provide open source code for computing **Characteristic Modes amidst Arbitrary Backgrounds** (BCMs)  of two-sphere systems based on scattering matrices.

The three MATLAB code files below correspond to the examples in the paper, and they can be \textcolor{red}{downloaded and run directly}. Note that the author's MATLAB version is R2020b, and if there are problems running them, the user may want to check the MATLAB version for consistency: 

1.  Case1_PEC.m; Calculate the dual PEC sphere system with sphere 1 of radius $a$ and separation distance $d=3a$.
2. Case2_Delectric.m; Calculate the dual dielectric sphere system with sphere 1 of radius $a$ and sphere 2 of radius $0.75a$. The separation distance $d=2.25a$. The relative permittivities of the two spheres are $\epsilon_1=2,\epsilon_2=8$, resp. Note that sphere 1 is set as background structure.
3. Case3_Composite; The sphere 1 is a PEC sphere with dielectric coating. the radius of this sphere is $0.8a$, the thickness of the coat is $\Phi=0.2a$, with relative permittivity $\epsilon_r=15$. Sphere 2 is a PEC sphere with radius $a$. The separation distance is  $d=3a$

### Notes:

The translation matrices of the vector spherical wave functions are calculated using C++. Note that we use C++17 language style. The source files are in the folder "cpp Files for Calculating Translation Matrix". We packaged them into the executable file "Trans_Solver.exe".
