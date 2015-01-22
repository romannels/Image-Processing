%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Roman Nelson                                                            %
% 5 October 2013                                                          %
%                                                                         %
% Program to simulate calibration of an imaging systemm                   %
% Implements: blackbox3()                                                 %
% INPUTS: 1024 x 1024 image matrix (p), 1024 x 1024 gain matrix, 1024 x   %
% 1024 offset matrix                                                      %
% OUTPUTS: 1024 x 1024 grayscale uint16 image matrix                      %
% Use program calibrate3 to determine gain and offset matrices            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  r = hw3correct(p,g,d)
    p = uint8(p); % convert input image to uint8
    avg = 150; % number of iterations to average scene to reduce noise
    gain = g; % assigns input "g" to variable "gain" to enhance readability
    offset = d; % assigns input "d" to variable "gain" to enhance 
                % readability
    r = double(zeros(1024,1024)); % creates zeros output matrix r to reduce
                                  % computation time
    imgnew = 0; %initializes variable "imgnew" to zero, used to average
                % images
    for i = 1:avg % for loop to average input scene
        img = double(blackbox3(p)); % creates new blackbox image
        imgnew = imgnew + img; % sums images together
    end
    
   img = imgnew / avg; % find avareage of images
   
    for m = 1:1024 % for loops for image correction (calibration)
        for n = 1:1024
           r(m,n) = (img(m,n) - offset(m,n)) / gain(m,n); % removes
           % offset and normalizes by the gain
        end
    end
    r = uint16(r); % converts output image to uint16 integer
    %disp('RMS = '); % debugging
    %r_m_s = rms1(p,r); % debugging
    %disp(r_m_s); % debugging
end

function r_m_s = rms1(img1, img2)
% function to find rms of two input images (values)
img1 = double(img1); % convert image one to type double
img2 = double(img2); % convet image two to type double
msqr = 0; % initialize variable msqr to zero
for m = 1:1024 % for loops to determine mean square
    for n = 1:1024
        sum = img1(m,n) - img2(m,n); % subtract pixel value from mean
        msqr = msqr + sum * sum; % sum previous msqr with current msqr
    end
end
msqr = msqr / (1024*1024); % divide by image size
r_m_s = sqrt(msqr); % find root msqr
end
