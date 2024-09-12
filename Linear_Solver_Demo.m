set(0,'defaulttextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex'); 
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultLineLineWidth',3);
set(0,'defaultAxesFontSize',35)

clf
N = 1000;
try_dense = true;
x = linspace(0,2*pi,N+2)';
x(1) = [];
x(end) = [];
dx = x(2)-x(1);

%%% exact solution and RHS
u_exact = sin(3*x);
f = -9*u_exact;

if try_dense
%%% Dense, 'exact' linear algebra
D2_dense = diag((-2/dx/dx)*ones(1,N)) + diag((1/dx/dx)*ones(1,N-1),1) + diag((1/dx/dx)*ones(1,N-1),-1);
% D2_dense(1,end) = (1/dx/dx); % for periodic BCs
% D2_dense(end,1) = (1/dx/dx); % for periodic BCs
tic
u_dense_slash = D2_dense\f;
disp('time for dense solve')
toc
else
u_dense_slash = u_exact;
end

%%% Sparse, 'exact' linear algebra
D2_sparse = spdiags(ones(N,1)*[1 -2 1]./(dx^2),-1:1,N,N);
% D2_sparse(1,end) = (1/dx/dx); % for periodic BCs
% D2_sparse(end,1) = (1/dx/dx); % for periodic BCs
tic
u_sparse_slash = D2_sparse\f;
disp('time for sparse solve')
toc

%%% Sparse itterative linear algebra
tic
max_its = N;
tol = 0.1*dx*dx; %try = 1e-14; 
[u_sparse_gmres,FLAG,RELRES,ITER] = gmres(D2_sparse,f,max_its,tol);
disp(['GMRES converged in ' num2str(ITER(2)) ' itterations and took '])
toc

%%% Matrix-free itterative linear algebra
tic
max_its = N;
tol = 0.1*dx*dx; %try = 1e-14; 
D2_times = @(u) Apply_D2(u,dx);
[u_MF_gmres,FLAG,RELRES,ITER] = gmres(D2_times,f,max_its,tol);
disp(['GMRES converged in ' num2str(ITER(2)) ' itterations and took '])
toc


plot(x,(u_exact-u_dense_slash),'-','linewidth',2)
hold all
plot(x,(u_exact-u_sparse_slash),':','linewidth',4)
hold all
plot(x,(u_exact-u_sparse_gmres),'o','markersize',12)
hold all
plot(x,(u_exact-u_MF_gmres),'.','markersize',25)
hold all
legend('error: dense','error: sparse','error: sparse, GMRES','error: matrix-free')
xlabel('x')
ylabel('Error')


function y = Apply_D2(u,dx)
    y = 0*u; % initialize output
    N = length(y);
    %%%%%%%%%%%%%%%%%%
    i = 1;
    i_next = i+1;
    y(i) = (1/(dx^2))*(-2*u(i) + u(i_next));
    %%%%%%%%%%%%%%%%%%
    i = N;
    i_prev = i-1;
    y(i) = (1/(dx^2))*(-2*u(i) + u(i_prev));
    %%%%%%%%%%%%%%%%%%
    for i = 2:N-1
        i_next = i+1;
        i_prev = i-1;
        y(i) = (1/(dx^2))*(-2*u(i) + u(i_next) + u(i_prev));
    end
end


function y = Apply_D2_periodic(u,dx)
    y = 0*u; % initialize output
    N = length(y);
    %%%%%%%%%%%%%%%%%%
    i = 1;
    i_next = i+1;
    i_prev = N;
    y(i) = (1/(dx^2))*(-2*u(i) + u(i_next) + u(i_prev));
    %%%%%%%%%%%%%%%%%%
    i = N;
    i_next = 1;
    i_prev = i-1;
    y(i) = (1/(dx^2))*(-2*u(i) + u(i_next) + u(i_prev));
    %%%%%%%%%%%%%%%%%%
    for i = 2:N-1
        i_next = i+1;
        i_prev = i-1;
        y(i) = (1/(dx^2))*(-2*u(i) + u(i_next) + u(i_prev));
    end
end