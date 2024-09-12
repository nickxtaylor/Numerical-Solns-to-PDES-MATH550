N = 32;
L = 2*pi;
h = L/N;
x = h*(1:N)'-L/2;
[X,Y] = meshgrid(x,x);

sig = 100;
H = @(x,y) sin(3*x).*cos(2*y);
H_x = @(x,y) 3*cos(3*x).*cos(2*y);
H_y = @(x,y) -2*sin(3*x).*sin(2*y);

%%%%%%%%%%%
kk = [0:N/2-1 0 -N/2+1:-1]';
kk2 = [0:N/2 -N/2+1:-1]';
ik = ((2*pi)/L)*1i*kk;
ik2 = ((2*pi)/L)*1i*kk2;

[ikX,ikY] = meshgrid(ik,ik);
[ik2X,ik2Y] = meshgrid(ik2,ik2);


fft_H_x = ifft2(ikX.*fft2(H(X,Y)));
fft_H_y = ifft2(ikY.*fft2(H(X,Y)));
%%%%%%%%%%
clf
colormap(turbo)
subplot(1,3,1)
pcolor(X,Y,H(X,Y));
hold all
title('$$H(x,y)$$')
xlabel('x')
ylabel('y')
colorbar
subplot(1,3,2)
pcolor(X,Y,((H_x(X,Y)-fft_H_x)));
title('Error in $$\frac{d}{dx} H(x,y)$$')
colorbar
xlabel('x')
ylabel('y')
subplot(1,3,3)
pcolor(X,Y,((H_y(X,Y)-fft_H_y)));
title('Error in $$\frac{d}{dy} H(x,y)$$')
colorbar
xlabel('x')
ylabel('y')