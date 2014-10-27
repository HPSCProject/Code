mtx_cmp assumes the following:
1) Both the files are in the directory containing mtx_cmp and hence named differently
2) Dimensions of the matrix in both the files are same

Usage:

./mtx_cmp.exe <col(which is timestemp in our case)> <row> <File1> <File2>

Example:for 250 timesteps

./mtx_cmp.exe 250 2 X_x.dat Y_y.dat



