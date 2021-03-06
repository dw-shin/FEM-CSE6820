function [c4n, n4e, ind4e, inddb] = mesh_FEM2D_Tri_rectangle(xl, xr, yl, yr, Mx, My, k)
%% 
% mesh_FEM2D_Tri_rectangle    Mesh geometry on 2D triangular domain
%    mesh_FEM2D_Tri_rectangle(xl, xr, yl, yr, Mx, My, k) generates an 
%    uniform triangular mesh on the domain [xl,xr]x[yl,yr] in 2D with Mx 
%    elements along x-direction and My elements along y-direction. Also  
%    this code returns an index matrix for continuous k-th order polynomial 
%    approximations.
%
%    - Input
%      xl     x-coordinate of bottom-left vertex of the domain
%      xr     x-coordinate of top-right vertex of the domain
%      yl     y-coordinate of bottom-left vertex of the domain
%      yr     y-coordinate of top-right vertex of the domain
%      Mx     the number of elements along x-direction
%      My     the number of elements along y-direction
%      k      polynomial order for the approximate solution
%
%    - Output
%      c4n    coordinates for nodes.
%             c4n is a (k*Mx+1)*(k*My+1) dimensional vector and it contains
%             all coordinates for nodes of the approximate solution.
%      n4e    nodes for elements.
%             n4e is a 3-by-2*Mx*My matrix. Each column of n4e contains 
%             3 vertices of the corresponding element in a counterclockwise
%             orientation. The first two vertices are the end-points of the
%             longest edge.
%      n4db   nodes for Dirichlet boundary.
%             n4db is a 2*k*(Mx+My) dimensional vector and it contains node
%             number for Dirichlet boundary.
%      ind4e  indices for elements
%             ind4e is a (k+1)*(k+2)/2-by-2*Mx*My matrix. Each column of
%             ind4e contains indices into all nodes in the corresponding  
%             element. These nodes are ordered from left to right and from 
%             bottom to top, and the starting node is in the third row of n4e.

%% index
ind4e = zeros((k + 1)*(k + 2)/2, 2*Mx*My);
tmp = repmat((1:k:k*Mx)', 1, My) + repmat((0:k*(k*Mx + 1):((k*Mx + 1)*((My - 1)*k + 1) - 1)), Mx, 1);
tmp2 = repmat(((k + 1):k:(k*Mx + 1))', 1, My) + repmat((k*(k*Mx + 1):k*(k*Mx + 1):((k*Mx + 1)*(k*My))), Mx, 1);
tmp = tmp(:)';
tmp2 = tmp2(:)';
for j=1:k+1
    ind4e((1 + (j - 1)*(k + 2 - j/2)) + (0:(k + 1 - j)), 1:2:2*Mx*My) = ...
        repmat(tmp,k + 2 - j, 1) + repmat(((j - 1)*(Mx*k + 1)+(0:(k + 1 - j)))', 1, Mx*My);
    ind4e((1 + (j - 1)*(k + 2 - j/2)) + (0:(k + 1 - j)), 2:2:2*Mx*My) = ...
        repmat(tmp2, k + 2 - j, 1) + repmat((-(j - 1)*(Mx*k + 1) - (0:(k + 1 - j)))', 1, Mx*My);
end

%% n4e
n4e = ind4e([k + 1, (k + 1)*(k + 2)/2, 1],:);

%% indDb
inddb = [1:(k*Mx + 1), 2*(k*Mx + 1):(k*Mx + 1):(k*Mx + 1)*(k*My + 1), ...
                ((k*Mx + 1)*(k*My + 1) - 1):-1:(k*My*(k*Mx + 1) + 1), ((k*My - 1)*(k*Mx+1) + 1):-(k*Mx + 1):(k*Mx + 2)];

%% c4n
x = linspace(xl, xr, k*Mx + 1);
y = linspace(yl, yr, k*My + 1);
y = repmat(y, k*Mx + 1, 1);
x = repmat(x, k*My + 1, 1)';
c4n = [x(:), y(:)]';
end