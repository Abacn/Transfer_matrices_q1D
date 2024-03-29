# Transfer matrix method for quasi-one-dimensional hard sphere systems


Transfer matrix is a technique to solve the partition function of such simple systems in statistical mechanics where one unit (such as spin, particle) only interacts with its limited neighboring units.

The code here calculates the partition function of hard disks confined between two parallel lines (2D system), and hard spheres confined within cylindrical pores (3D systems). The algorithm here is exact if up to next-nearest-neighbor along the axial direction can interact.

## Summary

- Language: Matlab (work with Matlab R2016a, R2017a)

- Publications

  - Yi Hu, Lin Fu and Patrick Charbonneau. "Correlation lengths in quasi-one-dimensional systems via transfer matrices." Molecular Physics, 21, 3345 (2018). http://doi.org/10.1080/00268976.2018.1479543 (arXiv:1804.00693)

  - Yi Hu and Patrick Charbonneau. "Comment on 'Kosterlitz-Thouless-type caging-uncaging transition in a quasi-one-dimensional hard disk system'." Physical Review Research, 3, 038001 (2021). https://doi.org/10.1103/PhysRevResearch.3.038001 (arXiv:2009.11194)

## File Description

### prefixes description
- NN : nearest-neighbor cases
- NNN : next-nearest-neighbor cases
- 2D : 2D system, hard disks confined between two parallel lines
- 3D : 3D system, hard spheres confined within cylindrical pores

### File specification
- \*Calc.m : run scripts, need param.dat (2D) or param3D.dat (3D) as input
- \*Eigs.m : calculating eigenvalue and eigenvectors
- \*Matrix.m : return M\*x for given x, where M is the transfer matrix
- \*TMatrix.m : return M^T\*x for given x, where M^T is the transpose of the transfer matrix
- \*Density.m : get the axial density rho=N/L 
- \*Cor.m : calculating correlation functions g(|i-j|)
- \*gofr : Planting equilibrium configuration and generating g(r) statistics
- NN3DReduce.m : obtaining eigenvalue and eigenvector from reduced transfer matrix. The angular integral has been calculated separately (see reference).
- redoxis.m : Calculate xi from g(|i-j|), need gs.dat, which is generated by \*Cor.m

### Input files
param.dat : for 2D systems. 4 lines in total.

1. ny (integer): Number of grid
2. ymax (float): the disk center are confined between (-ymax/2, ymax/2). The line seperation is ymax+1
3. betaF (float): the reduced pressure betaF=beta*P*ymax
4. c (float): the parameter used by change-of-variable scheme. More samples distributed in the regime close to the boundary in larger c

param3D.dat : for 3D systems. 4 lines in total.

1. nr (integer): Number of grids for radial part
2. nt (integer): Number of grid for angular part
3. dmax (float): the radial coordinates of sphere center are confined between (0, dmax/2). The cylinder diameter seperation is dmax+1
4. betaF (float): the reduced pressure betaF=beta\*P\*pi\*(ymax/2)^2

