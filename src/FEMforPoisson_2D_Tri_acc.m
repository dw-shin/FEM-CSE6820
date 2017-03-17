function [u, A, b, fns] = FEMforPoisson_2D_Tri_acc(c4n, n4e, n4db, ind4e, M_R, Srr_R, Srs_R, Ssr_R, Sss_R, f, u_D)
%% 
% FEMforPoisson_2D_Tri_acc    accelerated FEM solver for Poisson problem 
%                             in 2D with triangular elements
%    FEMforPoisson_2D_Tri_acc(c4n,n4e,n4db,ind4e,M_R,Srr_R,Srs_R,Ssr_R,Sss_R,f,u_D) 
%    solves the Poisson problem. This code is much faster than 
%    FEMforPoisson_2D_Tri.
%    In order to use this code, mesh information (c4n, n4e, n4db, ind4e),
%    matrices (M_R, Srr_R, Srs_R, Ssr_R, Sss_R), the source f, and the  
%    boundary condition u_D. Then the results of this code are the numerical  
%    solution u, the global stiffness matrix A, the global load vector b 
%    and the freenodes.
%
%    - Input
%      c4n    coordinates for nodes.
%             c4n is a N dimensional vector and it contains all
%             coordinates for nodes of the approximate solution.
%      n4e    nodes for elements.
%             n4e is a matrix with 3 rows. Each column of n4e contains 
%             3 vertices of the corresponding element in a counterclockwise
%             orientation. The first two vertices are the end-points of the
%             longest edge.
%      n4db   nodes for Dirichlet boundary.
%             n4db is a vector and it contains all node number on the 
%             Dirichlet boundary.
%      ind4e  indices for elements
%             ind4e is a matrix with (k+1)(k+2)/2 rows. Each column of ind4e 
%             contains indices into all nodes in the corresponding element. 
%             These nodes are ordered from left to right and from bottom to
%             top, and the starting node is in the third row of n4e.
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
[u, A, b, fns] = FEMforPoisson_2D_Tri_p(c4n, n4e, n4db, ind4e, M_R, Srr_R, Srs_R, Ssr_R, Sss_R, f, u_D);
end