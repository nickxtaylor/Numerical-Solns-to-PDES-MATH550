function w_hat = AntiAlias(u_hat, v_hat)
	N = length(u_hat);
	M = 3*N/2; % 3/2th rule
	u_hat_pad = [u_hat(1:N/2) zeros(1, M-N) u_hat(N/2+1:end)];
	v_hat_pad = [v_hat(1:N/2) zeros(1, M-N) v_hat(N/2+1:end)];
	u_pad = ifft(u_hat_pad);
	v_pad = ifft(v_hat_pad);
	w_pad = u_pad.*v_pad;
	w_pad_hat = fft(w_pad);
	w_hat = 3/2*[w_pad_hat(1:N/2) w_pad_hat(M-N/2+1:M)];
end % AntiAlias()