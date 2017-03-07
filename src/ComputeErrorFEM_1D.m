function error = ComputeErrorFEM_1D(c4n, ind4e, M_R, D_R, u, Du)
%% 
% ComputeErrorFEM_1D    Semi H1 error 
%    ComputeErrorFEM_1D(c4n,ind4e,M_R,D_R,u,Du) computes the semi H1 error 
%    between the exact solution and the FE solution.
%
%    - Input
%      c4n    coordinates for nodes.
%             c4n is a N dimensional vector and it contains all
%             coordinates for nodes of the approximate solution.
%      ind4e  indices for elements
%             ind4e is a (k+1)-by-M matrix. Each column of n4db contains 
%             indices into all nodes in the corresponding element from 
%             the left-end point to the right-end point.
%      M_R    Mass matrix on the reference interval
%      D_R    Differentiation matrix on the reference interval
%      u      numerical solution
%             u is a N dimensional vector and it contains all values
%             corresponding to the nodes in c4n.
%      Du     Derivative of the exact solution for the model problem
%
%    - Output
%      error  Semi H1 error between the exact solution and the FE solution.

%%
error = 0;
for j = 1:size(ind4e,2)
    J = (c4n(ind4e(end,j)) - c4n(ind4e(1,j))) / 2;
    De = Du(c4n(ind4e(:,j))') - D_R*u(ind4e(:,j)) / J;
    error = error + J*De'*M_R*De;
end
error = sqrt(error);
end