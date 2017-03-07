function [M_R, S_R, D_R] = MatrixforPoisson_1D(k)
%% 
% MatrixforPoisson_1D    Matrices for Poisson using FEM in 1D
%    MatrixforPoisson_1D(k) generates the mass matrix M_R, the stiffness
%    matrix S_R and the differentiation matrix D_R for continuous k-th 
%    order polynomial approximations on the reference interval.
%
%    - Input
%      k      polynomial order for the approximate solution
%
%    - Output
%      M_R    Mass matrix on the reference interval
%      S_R    Stiffness matrix on the reference interval
%      D_R    Differentiation matrix on the reference interval

%%
if k==1
    M_R = [2 1; 1 2]/3;
    S_R = [1 -1; -1 1]/2;
    D_R = [-1 1; -1 1]/2;
elseif k==2
    M_R = [4 2 -1; 2 16 2; -1 2 4]/15;
    S_R = [7 -8 1; -8 16 -8; 1 -8 7]/6;
    D_R = [-3 4 -1; -1 0 1; 1 -4 3]/2;
elseif k==3
    M_R = [128 99 -36 19; 99 648 -81 -36;
        -36 -81 648 99; 19 -36 99 128]/840;
    S_R = [148 -189 54 -13; -189 432 -297 54;
        54 -297 432 -189; -13 54 -189 148]/80;
    D_R = [-11 18 -9 2; -2 -3 6 -1; 1 -6 3 2; -2 9 -18 11]/4;
else
    M_R = 0; S_R = 0; D_R = 0;
end
end