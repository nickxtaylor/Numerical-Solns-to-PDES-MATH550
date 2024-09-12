clear all

N = 2048;
L = 2;
h = L/N;
x = h*(1:N)'-L/2;

H_0 = @(x) max(0,1-L*abs(x));
dH_0 = @(x) -sign(x).*L.*H_0(x)./(1-L*abs(x));

%%%%%%%%%%%
kk = [0:N/2-1 0 -N/2+1:-1]';
ik = ((2*pi)/L)*1i*kk;

f_d_H = ifft(ik.*fft(H_0(x)));

H_hat = fft(H_0(x));
H_tilde = ifft(H_hat);

x_plot = linspace(-L/2,L/2,5*N);
H_sum = 0*x_plot;
for i = 1:N
   H_sum = H_sum + (1/N)*H_hat(i).*exp(1i*2*pi*kk(i)*(x_plot-L/2-h)./L); 
end
    

%%%%%%%%%%
clf
subplot(1,2,1)
plot(x,H_0(x),'-o');
hold all
plot(x,H_tilde,'--p')
hold all
plot(x_plot,H_sum,'-')
xlim([-0.1 0.1])
subplot(1,2,2)
plot(x,dH_0(x),'-o');
hold all
plot(x,f_d_H,'--p')
