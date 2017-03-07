clear
tic
iter = 10;
a = 0;      % left-end point of the domain
b = 1;      % right-end point of the domain
k = 3;      % polynomial order for the approximate solution
M = 2.^(1:iter);        % the number of elements
f = @(x) pi^2*sin(pi*x);        % RHS in the Poisson problem
u_D = @(x) x*0;     % Dirichlet boundary condition for the solution u
Du = @(x) pi*cos(pi*x);     % Derivative of the exact solution for the model problem

error = zeros(1,iter);
time = error;
h = 1./M;
for j=1:iter
    tic;
    [c4n, n4e, n4db, ind4e] = mesh_FEM1D(a, b, M(j), k);
    
    [M_R, S_R, D_R] = MatrixforPoisson_1D(k);
    
%     u = FEMforPoisson_1D_acc(c4n, n4e, n4db, ind4e, M_R, S_R, f, u_D);
    u = FEMforPoisson_1D(c4n, n4e, n4db, ind4e, M_R, S_R, f, u_D);
    time(j) = toc(tic);
%     figure; plotFEM_1D(c4n, ind4e, u, k, ceil(20/j));
    error(j) = ComputeErrorFEM_1D(c4n, ind4e, M_R, D_R, u, Du);
end
rateE = (log(error(2:end)) - log(error(1:end-1))) ./ (log(h(2:end)) - log(h(1:end-1)));
disp(rateE)
toc
