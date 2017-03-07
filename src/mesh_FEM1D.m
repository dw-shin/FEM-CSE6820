function [c4n, n4e, n4db, ind4e]=mesh_FEM1D(a, b, M, k)
%% 
% mesh_1D    Mesh geometry on 1 dimensional domain
%    mesh_1D(a, b, M, k) generates an uniform mesh on the domain [a,b] in 
%    1D with mesh size h = 1/M. Also this code returns an index matrix for
%    continuous k-th order polynomial approximations.
%
%    - Input
%      a      left-end point of the domain
%      b      right-end point of the domain
%      M      the number of elements
%      k      polynomial order for the approximate solution
%
%    - Output
%      c4n    coordinates for nodes.
%             c4n is a N dimensional vector and it contains all
%             coordinates for nodes of the approximate solution.
%      n4e    nodes for elements.
%             n4e is a 2-by-M matrix. Each column of n4e contains indices 
%             into the left-end and the right-end points of the
%             corresponding element.
%      n4db   nodes for Dirichlet boundary.
%             n4db is a 2 dimensional vector and it contains node number
%             for Dirichlet boundary.
%      ind4e  indices for elements
%             ind4e is a (k+1)-by-M matrix. Each column of ind4e contains 
%             indices into all nodes in the corresponding element from 
%             the left-end point to the right-end point.

%% 
nrNodes = k*M + 1; % the number of nodes on the mesh in terms of k and N
c4n = linspace(a, b, nrNodes); % or c4n = a:(b-a)/(k*N):b
n4e = [1:k:(nrNodes-1); (k+1):k:nrNodes];
n4db = [1, nrNodes];
ind4e = repmat(n4e(1,:), k+1, 1)+repmat((0:k)', 1, M);
end