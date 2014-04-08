# Set matrix size here (in # of rows/columns)
n = 100;
density = 0.7;

# Matrix generator for Conjugate Gradient algorithm
function A = cgmatrix(n,density)
A = sprandsym(n,density);
A = A + n*speye(n);
end

# Generate matrix and output to file
M = cgmatrix(n,density);
save -ascii cg.mat M

