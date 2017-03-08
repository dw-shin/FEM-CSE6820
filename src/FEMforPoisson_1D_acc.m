function [u, A, b, fns] = FEMforPoisson_1D_acc(c4n, n4e, n4db, ind4e, M_R, S_R, f, u_D)
%% 
% FEMforPoisson_1D_acc    accelerated FEM solver for Poisson problem in 1D
%    FEMforPoisson1D_acc(c4n,n4e,n4db,ind4e,M_R,S_R,f,u_D) solves the 
%    Poisson problem. This code is much faster than FEMforPoisson_1D.
%    In order to use this code, mesh information (c4n, n4e, n4db, ind4e),
%    matrices (M_R, S_R), the source f, and the boundary condition u_D.
%    Then the results of this code are the numerical solution u, the
%    global stiffness matrix A, the global load vector b and the freenodes.
%
%    - Input
%      c4n    coordinates for nodes.
%             c4n is a N dimensional vector and it contains all coordinates
%             for nodes of the approximate solution.
%      n4e    nodes for elements.
%             n4e is a 2-by-M matrix. Each column of n4e contains indices 
%             into the left-end and the right-end points of the
%             corresponding element.
%      n4db   nodes for Dirichlet boundary.
%             n4db is a 2 dimensional vector and it contains node number
%             for Dirichlet boundary.
%      ind4e  indices for elements
%             ind4e is a (k+1)-by-M matrix. Each column of n4db contains 
%             indices into all nodes in the corresponding element from 
%             the left-end point to the right-end point.
%      M_R    Mass matrix on the reference interval
%      S_R    Stiffness matrix on the reference interval
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
[u, A, b, fns] = FEMforPoisson_1D_p(c4n, n4e, n4db, ind4e, M_R, S_R, f, u_D);
end