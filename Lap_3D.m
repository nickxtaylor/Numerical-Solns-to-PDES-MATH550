N = 128;
x = linspace(0,2*pi,N+2)';
x(1) = [];
x(end) = [];
dx = x(2)-x(1);

D2 = spdiags(ones(N,1)*[1 -2 1]./(dx^2),-1:1,N,N);

In = speye(N);
tic 
Lap = kron(In,kron(In,D2)) + kron(In,kron(D2,In)) + kron(kron(D2,In),In);
toc

tic
[X,Y,Z] = meshgrid(x,x,x);
toc

tic
f = -3.0*sin(X).*sin(Y).*sin(Z);
u_exact = sin(X).*sin(Y).*sin(Z);
toc

%tic 
%sol = Lap\f(:);
%toc

tic
max_its = N;
tol = 0.1*dx*dx; %try = 1e-14; 
[u_sparse_gmres,FLAG,RELRES,ITER] = gmres(Lap,f(:),[],tol);
disp(['GMRES converged in ' num2str(ITER(2)) ' itterations and took '])
toc

u_approx = reshape(u_sparse_gmres,N,N,N);
Rel_L2_err = sqrt(mean((u_exact(:) - u_approx(:)).^2))./sqrt(mean(u_exact(:).^2));
disp(['sparse relative error: ' num2str(Rel_L2_err)])

tic
max_its = N;
tol = 0.1*dx*dx; %try = 1e-14; 
MatVec = @(x) Apply_Lap(x,dx,N);
[u_matVec_gmres,FLAG,RELRES,ITER] = gmres(MatVec,f(:),[],tol);
disp(['mat vec GMRES converged in ' num2str(ITER(2)) ' itterations and took '])
toc

u_approx_MV = reshape(u_matVec_gmres,N,N,N);
Rel_L2_err_matVec = sqrt(mean((u_exact(:) - u_approx_MV(:)).^2))./sqrt(mean(u_exact(:).^2));
disp(['Mat. vec. relative error: ' num2str(Rel_L2_err_matVec)])

function y = Apply_Lap(u,dx,N)
    u = reshape(u,N,N,N); % reshape from vector to scalar
    y = 0*u; % initialize output
    %%%%%%%%%%%%%%%%%%
    for i = 1:N
        for j = 1:N
            for k = 1:N
                u_E = 0.0;
                u_W = 0.0;
                u_N = 0.0;
                u_S = 0.0;
                u_T = 0.0;
                u_B = 0.0;
                u_C = u(i,j,k);
                if i > 1; u_W = u(i-1,j,k); end
                if i < N; u_E = u(i+1,j,k); end
                if j > 1; u_S = u(i,j-1,k); end
                if j < N; u_N = u(i,j+1,k); end
                if k > 1; u_B = u(i,j,k-1); end
                if k < N; u_T = u(i,j,k+1); end
                y(i,j,k) = (1/(dx^2))*(u_W + u_E + ...
                                   u_S + u_N + ...
                                   u_T + u_B + ...
                                   -6.0*u_C);
            end
        end
    end
    y = y(:);
end