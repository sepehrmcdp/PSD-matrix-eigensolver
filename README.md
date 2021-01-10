# PSD-matrix-eigensolver

C code for calculation of eigenvalues and eigenvectors of a n*n real positive semi-definite matrix.
Implementation based on "Numerical Recipes in C" (Cambridge Univ. Press) by W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery.  

I was working on a project which involved getting eigenvectors and eigenvalues of a psd matrix (in c code) and I searched the web and couldn't find appropriate independent code so I decided to make my own version for better understanding.  
This algorithm uses householder reduction and then QL algorithm. This is algorithm is very efficient for calculating all values especially for middle sized matrices.  
functions work with linearized (1D) input and output arrays.  
There are two top function eig_dec_full (output eigenvalues and eigenvectors) and eig_dec_lite (output only eigenvalues, less runtime) which can be simply used by giving input matrix and sufficient allocated matrices.
Note that this implementation is pure single core and doesn't utilize any vector processing. For high performance parallel implementations check out MKL or similar libraries.
