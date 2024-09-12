set(0,'defaulttextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex'); 
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultLineLineWidth',3);
set(0,'defaultAxesFontSize',35)

clc
clear
clf

% solve u''(x) = cos(x) on domain 0 <= x <= pi
% BCs: u'(0) = 0, u(pi) = 1. 
% note neumann condition on left end point, dirchilet condition on right

N = 16;
x = linspace(0, pi, N);
dx = x(2)-x(1);

x_interior = x(1:end-1);

% transpose on x since we want a column vector
f = cos(x_interior');

% construct D2 matrix
e = ones(N-1, 1);
D2 = (1/dx^2)*spdiags([e, -2*e, e], -1:1, N-1, N-1);

% neumann boundary condition on left end point
% ghost points: do centered difference approx of u'(x_i) for i = 0 (boundary points)
% gives u'(0) = (u_1 - u_{-1})/(2dx).
% solve for u_{-1} and use BC to get u_{-1} = u_1
% substitute into D^2[u_0] = 2u_1 - 2u_0.
D2(1,2) = 2/dx^2; % set u_1 coeff to 2

% modify right end point with dirchlet condition as before
BC = zeros(N-1, 1);
BC(end) = 1;

RHS = f - (1/dx^2)*BC;

% solve
u = D2\RHS;
u_full = [u; 1]; % add on BC

% compare
sol = -cos(x); % analytic solution

hold on
plot(x, sol);
plot(x, u_full, '--p','linewidth',3,'markersize',10);
legend("Exact solution", "Numerical solution", Location="northwest")
xlim([0, pi])