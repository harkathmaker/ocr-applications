# This script generates a symmetric, positive-definite matrix.
# This script can be modified with the file name after 'save -ascii',
# the amount of rows/columns by n, and the matrix density by the
# density value. Configure those values to your liking, then run
# this script in matlab or octave to export a file containing the matrix
# in ASCII format. Use this matrix as the file name parameter in the
# HPCG program.

# Example:   octave ./cgmatrix.m 

# Set matrix size here (in # of rows/columns)
n = 1000;
# Set matrix density here (how dense the matrix is)
density = 0.2;

# Matrix generator for Conjugate Gradient algorithm
function A = cgmatrix(n,density)
A = sprandsym(n,density);
A = A + n*speye(n);
end

# Generate matrix and output to file
M = cgmatrix(n,density);

# Configure file name here after 'save -ascii'
save -ascii cg1000_s.mat M

