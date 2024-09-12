set(0,'defaulttextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex'); 
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultLineLineWidth',3);
set(0,'defaultAxesFontSize',35)

clc
clear

% solve u''(x) = cos(x) on domain 0 <= x <= pi
% BCs: u(0) = -1, u(pi) = 1. note dirichlet BCs

N = 64;
x = linspace(0, pi, N);
dx = x(2)-x(1);

% dirichlet BCs: ignore endpoints for now
x_interior = x(2:end-1);
f = cos(x_interior');

% BCs
BC = zeros(N-2, 1);
BC(1) = -1;
BC(end) = 1;

% construct D2 matrix
e = ones(N-2, 1);
D2 = (1/dx^2)*spdiags([e, -2*e, e], -1:1, N-2, N-2);

% modify with BCs
RHS = f - (1/dx^2)*BC;

% solve
u = D2\RHS;

% add BCs of solution back in
u_full = [-1; u; 1];

% compare
sol = -cos(x); % analytic solution

plot(x, sol);
hold on
plot(x, u_full, '--p','linewidth',3,'markersize',10);
legend("Exact solution", "$D^2 \symbol{92} f$", Location="northwest")
xlim([0, pi])