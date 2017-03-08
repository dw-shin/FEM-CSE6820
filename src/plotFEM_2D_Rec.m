function plotFEM_2D_Rec(c4n, ind4e, u, k, order)
%%
% plotFEM_2D_Rec    2D FE solution plot with rectangular elements
%    plotFEM_2D_Rec(c4n, ind4e, u, k, order) displays finite element 
%    solution u. In each element, the piecewise linear intepolation Iu is 
%    obtained from the solution corresponding to (k+order)-th order uniform
%    nodes. In order to obtain a smooth solution for the high-order case,
%    order should be large  enough.
%
%    plotFEM_2D_Rec(c4n, ind4e, u, k) displays finite element solution u
%    without a interpolation. In this case, the solution may not be smooth.
%
%     - Input
%      c4n    coordinates for nodes.
%             c4n is a N dimensional vector and it contains all
%             coordinates for nodes of the approximate solution.
%      ind4e  indices for elements
%             ind4e is a matrix with (k+1)^2 rows. Each column of ind4e 
%             contains indices into all nodes in the corresponding element 
%             from left to right and from bottom to top.
%      u      numerical solution
%             u is a N dimensional vector and it contains all values
%             corresponding to the nodes in c4n.
%      k      polynomial order for the approximate solution
%      order  the additional polynomial order in each element
%             If order>0, (k+order)-th uniform nodes are generated in each
%             element. In this case, the piecewise linear interpolation is
%             computed in each element.

%%
if nargin > 4
    plotFEM_2D_Rec_p(c4n, ind4e, u, k, order)
else
    plotFEM_2D_Rec_p(c4n, ind4e, u, k)
end
end