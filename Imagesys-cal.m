%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Roman Nelson                                                            %
% 5 October 2013                                                          %
%                                                                         %
% Program to determine gain and offset coeficient matrices of a given     %
% input image                                                             %
% Implements: blackbox3()                                                 %
% INPUTS: none                                                            %
% OUTPUTS: 1024 x 1024 gain matrix, 1024 x 1024 offset matrix             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [g,d] = hw3calibrate()

% x min = 45 (no zeros), determined via experimentation
% x max = 947(948-950), determined via experimentation
% MATHS
% darkoutput = 45*gain-offset
% brightoutput = 947*gain-offset
%
%               darkoutput =  45*gain-offset;
% -1*(45/947) * brightoutput = 45*gain + (45/947)*offset;
% darkoutput - (45/947) * brightoutput = offset((45/947) - 1);
% (darkoutput - (45/947) * brightoutput) / ((45/947) - 1); = offset;
%(darkoutput - A * brightoutput) / B; = offset;
% (darkoutput + offset) / 45 = gain


min = 45; % lowest value that gives all zeros
max = 940; % highest value that gives all 255 (max for 8 bit)
A = (min/max); % see above for derivation
B = ((min/max) - 1); % see above for derivation
dark = ones(1024,1024) * min; % darkfield image
bright = ones(1024,1024) * max; % brightfield image

dark_convert = double(blackbox3(dark)); % run darkfield through blackbox
bright_convert = double(blackbox3(bright)); % run brightfield through 
                                            % blackbox
offset = double(zeros(1024,1024)); % initialize offset matrix for speed
gain = double(zeros(1024,1024)); % intializes gain matrix for speed

for m = 1:1024 % loop to calculate calibration coeefficients
    for n = 1:1024
        offset(m,n) = double((dark_convert(m,n) - A * bright_convert(m,n)) / B);
        gain(m,n) = double((dark_convert(m,n) +offset(m,n)) / 45);
    end
end
g = gain; % assigns output variable "g" to variable gain
d = offset; % assigns output variable "d" to variable offset


% DEBUGGING
%disp('Max offset:');
%maxoff = max(d(:));
%disp(maxnoff);
%disp('Min offset');
%minoff = min(d(:));
%disp(minoff);
%disp('Max gain');
%maxgain = max(g(:));
%disp(maxgain);
%disp('Min gain');
%mingain = min(g(:));
%disp(mingain);
end
