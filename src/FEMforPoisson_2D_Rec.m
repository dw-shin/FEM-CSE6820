function [u, A, b, fns] = ...
    FEMforPoisson_2D_Rec(c4n, n4e, n4db, ind4e, M_R, Srr_R, Sss_R, f, u_D)
%% 
% FEMforPoisson_2D_Rec    FEM solver for Poisson problem in 2D with
%                         rectangular elements
%    FEMforPoisson_2D_Rec(c4n,n4e,n4db,ind4e,M_R,Srr_R,Sss_R,f,u_D) solves 
%    the Poisson problem. In order to use this code, mesh information (c4n,
%    n4e, n4db, ind4e), matrices (M_R, Srr_R, Sss_R), the source f, and the
%    boundary condition u_D. Then the results of this code are the  
%    numerical solution u, the global stiffness matrix A, the global load 
%    vector b and the freenodes.
%
%    - Input
%      c4n    coordinates for nodes.
%             c4n is a N dimensional vector and it contains all coordinates
%             for nodes of the approximate solution.
%      n4e    nodes for elements.
%             n4e is a matrix with 4 rows. Each column of n4e contains 
%             4 vertices of the corresponding element in a counterclockwise
%             orientation. The first vertex is the bottom-left node.
%      n4db   nodes for Dirichlet boundary.
%             n4db is a vector and it contains all node number on the 
%             Dirichlet boundary.
%      ind4e  indices for elements
%             ind4e is a matrix with (k+1)^2 rows. Each column of ind4e 
%             contains indices into all nodes in the corresponding element 
%             from left to right and from bottom to top.
%      M_R    Mass matrix on the reference interval
%      Srr_R  Stiffness matrix on the reference interval
%      Sss_R  Stiffness matrix on the reference interval
%      f      RHS in the Poisson problem
%      u_D    Dirichlet boundary condition for the solution u
%
%    - Output
%      u      numerical solution
%             u is a N dimensional vector and it contains all values
%             corresponding to the nodes in c4n.
%      A      Global stiffness matrix
%             A is a N-by-N matrix.
%      b      Global right-hand side
%             b is a N dimensional vector.
%      fns    free nodes
%             fns is a set of all nodes except boundary nodes.

%%
A = sparse(length(c4n), length(c4n));
b = zeros(length(c4n), 1);
u = b;
for j = 1:length(n4e)
    xr = (c4n(1,n4e(2,j))-c4n(1,n4e(1,j)))/2;
    ys = (c4n(2,n4e(4,j))-c4n(2,n4e(1,j)))/2;
    J = xr*ys;
    rx=ys/J; sy=xr/J;

    A(ind4e(:,j), ind4e(:,j)) = A(ind4e(:,j),ind4e(:,j)) ...
        + J*(rx^2*Srr_R + sy^2*Sss_R);
    b(ind4e(:,j)) = b(ind4e(:,j)) + J*M_R*f(c4n(:,ind4e(:,j))');
end
fns = setdiff(1:length(c4n), n4db);
u(ind4e) = u_D(c4n(ind4e(:)));
b = b - A*u;
u(fns) = A(fns, fns)\b(fns);
end