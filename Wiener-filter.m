%%% filter
% Creates a Wiener Filter given scene autocorrelation model parameter, 
% Gaussian system transfer model parameter, noise autocorelation paramter,
% and scene standard deviation

function f = filter(hpar,sstd,spar,npar)

hHat = zeros(256); % Memory allocation
for u = 1:256
    for v = 1:256
        hHat(u,v) = exp(-1 * (u^2 + v^2) / hpar(u,v));% Blackbox6 transfer function         
    end
end

phi_s = zeros(256); % Memory allocation
for m = 1:256
    for n = 1:256
        phi_s = spar(m,n) ^ (sqrt(m^2 + n^2)); % Create scene autocorrellation matrix
    end
end

MN = 256 ^2; % Number of elements in matrix
Phi_s = fft2(phi_s) * (1/ MN); % Inverse FFT of scene Autocorrelation (power spectrum)
fHat_w = conj(hHat) .* (1 ./ ((abs((hHat .* hHat))) + ((npar ^2) ./ ((sstd^2) * Phi_s)))); %API
f = fHat_w; % Assign fHat_w to f;
end
