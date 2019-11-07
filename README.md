# Periodic particle arrangements using standing acoustic waves
Code accompanying the paper "Periodic particle arrangements using standing acoustic waves", by [Fernando Guevara Vasquez](https://www.math.utah.edu/~fguevara) and China Mauck. For a preprint see: [arXiv:1908.08664](https://arxiv.org/abs/1908.08664). The published version is in (TBA).

This code was developed and tested with Matlab version R2018b. To re-generate the figures in the manuscript above:
* `code/fig2.m` - example of acoustic radiation potential (ARP) corresponding to a tetragonal lattice arrangement (see example 2.1)
* `code/fig3.m` - example of ARP corresponding to lines (see example 2.1)
* `code/fig4.m` - another ARP showing lines (see example 2.2)
* `code/fig5.m` - position of minima of ARP in 3D, showing planes (see example 2.3)
* `code/fig6.m` - position of minima of ARP in 3D, showing lines  (see example 2.4)
* `code/fig7.m` - representatives of the 3 achievable Bravais lattices in 2D (orthorhombic centered, hexagonal and tetragonal). Please see section 3 for what we mean by "achievable"
* `code/fig8.m` - representatives of the 6 achievable Bravais lattices in 3D (triclinic primitive, orthorombic face-centred, trigonal primitive, cubic primitive, cubic face-centred and cubic body-centred). Please see section 3 for what we mean by "achievable"

The following Matlab packages are dependencies. The two arrow packages (`arrow` and `Arrow3`) are used to display plane wave directions in each figure. The `Inhull` package is used to find points in a unit cell for figure 8. The packages are available at:
* [arrow](https://www.mathworks.com/matlabcentral/fileexchange/278-arrow) by Erik Johnson
* [Arrow3](https://www.mathworks.com/matlabcentral/fileexchange/14056-arrow3) by Tom Davis
* [Inhull](https://www.mathworks.com/matlabcentral/fileexchange/10226-inhull) by John D'Errico

Each of the m-files above generates a `PNG` file that is cropped to remove whitespace with a `system` call to `mogrify`. These `system` calls can be safely commented out. To remove whitespace we use `mogrify` command, part of the open source image manipuation suite [ImageMagick](https://imagemagick.org).
