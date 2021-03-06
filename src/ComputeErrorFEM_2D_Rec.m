function error = ComputeErrorFEM_2D_Rec(c4n, n4e, ind4e, M_R, Dr_R, Ds_R, u, ux, uy)
%% 
% ComputeErrorFEM_2D_Rec    Semi H1 error (2D rectangular element)
%    ComputeErrorFEM_2D_Rec(c4n,n4e,ind4e,M_R,Dr_R,Ds_R,u,ux,uy) computes 
%    the semi H1 error between the exact solution and the FE solution.
%
%    - Input
%      c4n    coordinates for nodes.
%             c4n is a N dimensional vector and it contains all
%             coordinates for nodes of the approximate solution.
%      n4e    nodes for elements.
%             n4e is a 4-by-M^2 matrix. Each column of n4e contains indices 
%             into the left-end and the right-end points of the
%             corresponding element.
%      ind4e  indices for elements
%             ind4e is a (k+1)-by-M^2 matrix. Each column of n4db contains 
%             indices into all nodes in the corresponding element from 
%             the left-end point to the right-end point.
%      M_R    Mass matrix on the reference interval
%      Dr_R   Differentiation matrix along r-direction on the reference 
%             interval
%      Ds_R   Differentiation matrix along s-direction on the reference 
%             interval
%      u      numerical solution
%             u is a N dimensional vector and it contains all values
%             corresponding to the nodes in c4n.
%      ux     Derivative of the exact solution along x-direction for the 
%             model problem
%      uy     Derivative of the exact solution along y-direction for the 
%             model problem
%
%    - Output
%      error  Semi H1 error between the exact solution and the FE solution.

%%
error = 0;
for j=1:size(ind4e,2)
    xr = (c4n(1,n4e(2,j))-c4n(1,n4e(1,j)))/2;
    ys = (c4n(2,n4e(4,j))-c4n(2,n4e(1,j)))/2;
    J = xr*ys;
    rx = ys/J; sy = xr/J;
    
    Dex = ux(c4n(:,ind4e(:,j))') - rx*Dr_R*u(ind4e(:,j));
    Dey = uy(c4n(:,ind4e(:,j))') - sy*Ds_R*u(ind4e(:,j));
    error = error + J*(Dex'*M_R*Dex + Dey'*M_R*Dey);
end
error = sqrt(error);
end