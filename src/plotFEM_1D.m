function plotFEM_1D(c4n, ind4e, u, k, pts)
%%
% plotFEM_1D    1 dimensional FE solution plot
%    plotFEM_1D(c4n,n4e,u,index,k,pts) displays finite element solution u.
%    In each element, the piecewise linear intepolation Iu is obtained from
%    the solution corresponding to (k+pts+1) uniform nodes. In order to 
%    obtain a smooth solution for the high-order case, pts should be large 
%    enough.
%
%    plotFEM_1D(c4n,n4e,u,index,k) displays finite element solution u
%    without a interpolation. In this case, the solution may not be smooth.
%
%     - Input
%      c4n    coordinates for nodes.
%             c4n is a N dimensional vector and it contains all
%             coordinates for nodes of the approximate solution.
%      ind4e  indices for elements
%             ind4e is a (k+1)-by-M matrix. Each column of n4db contains 
%             indices into all nodes in the corresponding element from 
%             the left-end point to the right-end point.
%      u      the FE solution
%             u is a N dimensional vector and it contains all values
%             corresponding to c4n.
%      k      polynomial order for the approximate solution
%      pts    the number of additional points in each element
%             If pts > 0, (k+pts+1) uniform nodes are generated in each
%             element. In this case, the piecewise linear interpolation is
%             computed in each element.

%%
if nargin > 4
    plotFEM_1D_p(c4n, ind4e, u, k, pts)
else
    plotFEM_1D_p(c4n, ind4e, u, k)
end
end