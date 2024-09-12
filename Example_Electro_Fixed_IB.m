set(0,'defaulttextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex'); 
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultLineLineWidth',3);
set(0,'defaultAxesFontSize',35)


clf
clear all
%%% Setup grid in x-y
N = 128; % number of grid point along one direction
L = 1.0;
x = linspace(0,L,N+2); % type 1 grid
dx = x(2)-x(1);

y = x;
dy = dx;


[X,Y] = meshgrid(x,y); % make 2D grid
[Xint,Yint] = meshgrid(x(2:end-1),y(2:end-1)); % make 2D grid of interior points

%%% make finite difference laplacian

% dirichelet boundary conditions
e = (1/dx^2)*ones(N,1);
D2 = spdiags([e -2*e e], -1:1, N, N);

colormap('turbo')

I_n = speye(N);

Lap = kron(I_n, D2) + kron(D2, I_n);

[sx,sy,sz] = sphere(50);

dLap = decomposition(Lap);
%%%%%%%%%%%%
% Make immersed boundary mats
xib = [0.5*L 0.25*L 0.3*L];
yib = [0.75*L 0.5*L 0.3*L];
%%%%% Plot IB stuff
figure(1)
q = [-3 1 2]; %sum(Q(:))
delta_a = @(r,a) (1/(2*pi*a^2))*exp(-0.5*(r/a).^2); 
delta = @(r) delta_a(r,1.2*dx);

Q = spreadQ(Xint,Yint,xib,yib,q,delta);

meanQ = sum(Q(:))*dx*dy; %;
Phi = reshape(dLap\Q(:),N,N);
Phi_ib = interpPhi(Xint,Yint,xib,yib,Phi,delta);
subplot(1,2,1)
surf(Xint,Yint,Q)
xlabel('x')
ylabel('y')
title('charge  $$\displaystyle \sum_i q_i \delta\left(r_i\right)$$')
subplot(1,2,2)
sf = surf(Xint,Yint,Phi);
hold all
plot3(xib,yib,Phi_ib,'ob','markerfacecolor',[0.6,0.65,1.0],'markersize',25)
% a = 3*dx;
% for k = 1:3
%     sh = surface(a*sx+Xib(k,1),a*sy+Xib(k,2),a*sz+Phi_ib(k));
%     set(sh,'edgecolor','none','facecolor',[0.6,0.65,1.0])
% end
view([0 90])
xlabel('x')
ylabel('y')
title('potential  $$\phi$$')

%%


function [Sq] = spreadQ(X,Y,xq,yq,q,delta)
    Sq = 0*X;
    Nq = length(q);
    for k = 1:Nq
        Rk = sqrt((X-xq(k)).^2 + (Y-yq(k)).^2);
        Sq = Sq + q(k)*delta(Rk);
    end
end

function [Jphi] = interpPhi(X,Y,xq,yq,Phi,delta)
    Jphi = 0*xq;
    Nq = length(xq);
    dx = X(1,2)-X(1,1);
    dy = Y(2,1)-Y(1,1);
    for k = 1:Nq
        Rk = sqrt((X-xq(k)).^2 + (Y-yq(k)).^2);
        Jphi(k) = dx*dy*sum(sum(Phi.*delta(Rk)));
    end
end