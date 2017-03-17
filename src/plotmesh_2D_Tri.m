function plotmesh_2D_Tri(c4n, n4e, ~)
%%
% plotmesh_2D_Tri    Mesh plot in 2D with rectangular elements
%    plotmesh_2D_Tri(c4n,n4e) displays the 2D mesh geometry with 
%    rectangular defined in the 3-by-N matrix n4e as a mesh. Each column of 
%    n4e contains indices into the 3 vertices of the corresponding element.
%
%    plotmesh_2D_Tri(c4n,n4e,onoff) displays the 2D mesh geometry with 
%    rectangular defined in the 3-by-N matrix n4e as a mesh and all nodes 
%    for the approximate solution. 
%    In this case, onoff is an arbitrary number.
%
%     - Input
%      c4n    coordinates for nodes.
%             c4n is a N dimensional vector and it contains all
%             coordinates for nodes of the approximate solution.
%      n4e    nodes for elements.
%             n4e is a matrix with 3 rows. Each column of n4e contains 
%             3 vertices of the corresponding element in a counterclockwise
%             orientation. The first two vertices are the end-points of the
%             longest edge.

%%
nElems = size(n4e, 2);
X = reshape(c4n(1, n4e([1 2 3 1], :)), 4 ,nElems);
Y = reshape(c4n(2, n4e([1 2 3 1], :)), 4, nElems);
EU = zeros(size(X));
set(gcf, 'Color', 'w')
patch(X, Y, EU, 1, 'facecolor', 'none')
view(2)

if nargin > 2
    hold on
    plot(c4n(1, :), c4n(2, :), 'ro', 'markersize', 7)
    
    % node number
    for k = 1:size(c4n, 2)
        text(c4n(1, k), c4n(2, k), num2str(k, 'v_{%d}'),...
            'VerticalAlignment', 'bottom');
    end
    
    % element number
    for k = 1:size(n4e, 2)
        mid = (c4n(:, n4e(1, k)) + c4n(:, n4e(2, k)) + c4n(:, n4e(3, k)))/3;
        text(mid(1), mid(2), num2str(k,'T_{%d}'), 'Color', 'b');
    end
    hold off
end
end