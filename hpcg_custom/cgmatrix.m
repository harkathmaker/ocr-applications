# Set matrix size here (in # of rows/columns)
n = 950;
density = 0.9;

# Matrix generator for Conjugate Gradient algorithm
function A = cgmatrix(n,density)
A = sprandsym(n,density);
A = A + n*speye(n);
end

# Generate matrix and output to file
M = cgmatrix(n,density);
save -ascii cg950.mat M

