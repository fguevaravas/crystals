# Periodic particle arrangements using standing acoustic waves
Code accompanying the paper "Periodic particle arrangements using standing acoustic waves", by [Fernando Guevara Vasquez](https://www.math.utah.edu/~fguevara) and China Mauck. For a preprint see: [arXiv:1908.08664](https://arxiv.org/abs/1908.08664). The published version is in xxxx.

This code was developed and tested with Matlab version R2018b. To generate the figures:
* `fig2.m` - example of acoustic radiation potential (ARP) corresponding to a tetragonal lattice arrangement (see example 2.1)
* `fig3.m` - example of ARP corresponding to lines (see example 2.1)
* `fig4.m` - another ARP showing lines (see example 2.2)
* `fig5.m` - position of minima of ARP in 3D, showing planes (see example 2.3)
* `fig6.m` - position of minima of ARP in 3D, showing lines  (see example 2.4)
* `fig7.m` - representatives of the 3 achievable Bravais lattices in 2D (orthorhombic centered, hexagonal and tetragonal). Please see section 3 for what we mean by "achievable"
* `fig8.m` - representatives of the 6 achievable Bravais lattices in 3D (triclinic primitive, orthorombic face-centred, trigonal primitive, cubic primitive, cubic face-centred and cubic body-centred). Please see section 3 for what we mean by "achievable"

The following Matlab packages are needed for display purposes. The two arrow packages (arrow and Arrow3) are used to include plane wave directions in each figure. The Inhull package is used to find points in a unit cell for figure 8.
https://www.mathworks.com/matlabcentral/fileexchange/278-arrow
https://www.mathworks.com/matlabcentral/fileexchange/14056-arrow3
https://www.mathworks.com/matlabcentral/fileexchange/10226-inhull
