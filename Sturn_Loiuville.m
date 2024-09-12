N = 32;
L = 2*pi;
h = L/N;
x = h*(1:N)';

sig = 100;
H_0 = @(x) exp(-sig*x.^2);
c_x = sin(3*pi*x./L);

kk = (2*pi/L).*[0:N/2-1 0 -N/2+1:-1]';
kk2 = (2*pi/L).*[0:N/2 -N/2+1:-1]';

H_hat = fft(H_0(x));
H_Nby2 = H_hat(N/2+1);
H_prime = ifft(1i.*kk.*H_hat);

c_X_H = c_x.*H_prime;

c_X_H_hat = fft(c_X_H);
L_H = 1i.*kk.*c_X_H_hat;
L_H(N/2+1) = mean(c_x)*(pi*N/L)^2*H_Nby2;




