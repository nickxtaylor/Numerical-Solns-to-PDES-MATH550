set(0,'defaulttextInterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultLineLineWidth',3);
set(0,'defaultAxesFontSize',35)

clc
clear all

N = 64;
L = 10;
h = L/N;
x = h*(1:N)'-L/2;

c = 100;
H_0 = @(x) -0.5*c*(sech(0.5*sqrt(c).*(x))).^2;
%H_0 = @(x) sech(10*x).^2;

%%%%%%%%%%%
kk = [0:N/2-1 0 -N/2+1:-1]';
ik = ((2*pi)/L)*1i*kk;

H_hat = fft(H_0(x));
H_hat_ai = fft(H_0(x));

dt = 0.00004;

t_final = 5;
N_steps = round(t_final/dt);


A_op = ik.^3;

for n = 1:N_steps
    Integrate_Factor = (1.0./(A_op + 1e-8)).*(exp(dt*A_op)-1);
    B_H = -3*ik.*fft(ifft(H_hat).^2);
    B_H_ai = -3*ik.*AntiAlias(H_hat_ai', H_hat_ai')';
    
    
    H_hat = exp(dt*ik.^3).*H_hat ...
          + Integrate_Factor.*B_H;
      
    H_hat_ai = exp(dt*ik.^3).*H_hat_ai ...
          + Integrate_Factor.*B_H_ai;
      
   plot(x,real(ifft(H_hat)),'-o')
   hold all
   plot(x,real(ifft(H_hat_ai)),'-p')
   hold off
   leg = legend('no AI','with AI');
   set(leg,'location','southeast')
   ylim([min(H_0(x))-0.5*max(abs(H_0(x))) max(H_0(x))+0.5*max(abs(H_0(x)))])
   drawnow
end

