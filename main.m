% AM900 A3
% This is the main script for computing the L2 error and order of
% convergence for Space-time Discontinuous Galerkin for advection

% TEST PARAMETERS
% See advSTDG.m for different options
movingmesh = 'full'
testcase = 'full'
% LEVELS OF REFINEMENT
N = [2 4 8 16 32 64 128 256 512 1024];
% COMPUTE ERRORS FOR DIFFERENT MESHES
err = zeros(size(N));
for i = 1:length(N)
    [u,approx_err] = advSTDG(N(i),N(i),1,movingmesh,testcase,false);
    err(i) = approx_err;
end
% COMPUTE ORDER OF CONVERGENCE
order = zeros(size(N));
for i = 2:length(err)
    order(i) = log2(err(i-1)/err(i));
end
disp('order of convergence in L2 norm')
disp(order')
% DEMONSTRATE SOLUTION
[u,approx_err] = advSTDG(64,64,1,movingmesh,testcase,true);
