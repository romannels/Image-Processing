%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------Roman Nelson---------------------------------%
%                                                                         %
% Noise filter                                                            %
%-------------------------------------------------------------------------%

function F = filter5()
flat_field = 128 * ones(256,256);


% average noisy images
% img_new = double(zeros(256,256));
% f_img_new = 0;
% for i = 1:128
%     img = double(blackbox5(flat_field));
%     img_new = double(img+img_new);
%     img = 0;
%     f_img = fft2(double(flat_field));
%     f_img_new = f_img + double(f_img_new);
% end
% img_new = (img_new / 128);
% mean_img_new = mean(img_new(:));
% std_img_new = std(img_new(:));
% f_img_new = double(f_img_new) / 128;
%disp(f_img_new);
%imshow(f_img_new)


% disp(mean_img_new);
% disp(std_img_new);

F_new = zeros(256,256);
k = 45;
%average the noisy pictures FFT's together and create filter
for i = 1:k
%convert to frequency domain
Z = abs(fft2(double(blackbox5(flat_field))));

%imshow(uint8(img_new));
%plot(2:length(Z),abs(Z(2:256)));
%figure

% initialize filter output
F = ones(256,256);
%set DC component to one (pass)
F(1,1) = 1;
% loop through transformed image to find noisy frequencies
for n = 2:length(Z)
    for m = 2:length(Z)
        
        % most noise is greater than 1e-13
    if Z(n,m) >= 1e-2;
        %destroy noisy frequencies
    F(n,m) = 0;
    else
        % pass all non-noisy frequencies
        F(n,m) = 1;
    end    
    end
end
%plot(F);
%imshow(F);

F_new = F + F_new;
end
F = F_new;
x = ifft2(F_new);
%imshow(x);
%figure
%plot(F);
%figure

% Normalize to maximum of one

% Set DC to one
for v = 1:256
    F(1,v) = 1;
    F(v,1) = 1;
end
% loop through to set values below threshold to zero
for n = 2:256
    for m = 2:256
        
        % if F has never had noise at this frequency
    if F(n,m) == k;
    F(n,m) = 1;
    else
        % pass all non-noisy frequencies
        F(n,m) = 0;
    end    
    end
end
%plot(F);
end
