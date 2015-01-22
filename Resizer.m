%-------------------------------------------------------------------------%
% Digital Image Aquisition and Display 
% Roman Nelson
% 1 September 2013
% 
% DESCRIPTION
% Function that takes an intensity image with unit8 data, reduces spatial
% and grayscale resolution as specified by the user and then reconstructs
% spatial resolution using nearest neighboor resampling and remaps
% grayscale resolution using requantization. 
%
% INPUT: s, grayscale, unit, input image
%       M, number of rows in simulated image
%       N, number of columns in simulated image
%       L, number of intensity levels in the simulated image
%
% OUTPUT: r, uint8, grayscale image same size as input image s
%-------------------------------------------------------------------------%
function r = hw1(s, M, N, L) 
%Digital Image Acquistion Simulation
p = imread(s); %push image into an array;
[Ms,Ns] = size(p); %size of input image array
Ls = max(p(:)) + 1; %maximum grayscale value of input image array
sx = M / Ms;
sy = N / Ns;
Lp = L;
rsim = uint8(zeros(M,N));
pnew = uint8(zeros(Ms,Ns));

%call the reduction function
rsim = reduce(M, N, Ms, Ns, Ls, Lp, sx, sy, p);

% Displays "camera" resolution
[newx, newy] = size(rsim);
disp('Y Dimension:');
disp(newx);
disp('X Dimension:');
disp(newy);
disp(Ls)

%Digital Image Display Simulation
%pnew = expand(Ms, Ns, M, N, sx, sy, rsim);
% nested for loop to rescale image
% NOTE: COULD NOT GET 'expand()' FUNCTION TO WORK
% WOULD ONLY OUTPUT BLACK IMAGE 
for ms = 1:Ms
    mp(ms)= imgmod(((ms + 0.5) * sx), M) + 1;
    for ns = 1:Ns
        np(ns) = imgmod(((ns + 0.5) * sy),N) + 1;
        rsimprime = rsim(mp(ms),np(ns));
        ls = rsimprime * 1/(((Ls / Lp)));
        pnew(ms, ns) = ls;
    end
end
%print the image

disp('Original Y Dimension:');
disp(Ns);
disp('Original X Dimension:');
disp(Ms);


x = imshow(pnew);

end


% function to reduce image size
function rsim = reduce(M, N, Ms, Ns, Ls, Lp, sx, sy, p)
  for mp = 1:M    
    ms(mp)= imgmod(((mp + 0.5) / sx),Ms) + 1; %picks out pixels and ensures integer values using imgmod function
   for np = 1:N
       ns(np) = imgmod(((np + 0.5) / sy),Ns) + 1; %picks out pixels and ensures integer values using imgmod function
       pprime = p(ms(mp),ns(np)); %assings variable to matrix p to make code easier to read
        lp = pprime/(Ls / Lp); %converts gray scale range 
        rsim(mp,np) = lp; %assigns new pixel value to the "cameras" matrix
   end
  end
end

% function to increase resolution to original resolution
function pnew = expand(Ms, Ns, M, N, sx, sy, rsim)
 for ms = 1:Ms
    mp(ms)= imgmod(((ms + 0.5) * sx), M) + 1; %picks out pixels and ensures integer values using imgmod function
    for ns = 1:Ns
        np(ns) = imgmod(((ns + 0.5) * sy),N) + 1;%picks out pixels and ensures integer values using imgmod function
        rsimprime = rsim(mp(ms),np(ns)); %assings variable to p to make code easier to read
        ls = rsimprime * 1/(((Ls / Lp)));%converts gray scale range back to original range
        pnew(ms, ns) = ls;%assigns new pixel value to the scene matrix
    end
 end
end

% function to enable periodic sampling and round non-integer values down (nearest neighbor sampling)
function x = imgmod(scaledinput,dim)
x = mod(floor(scaledinput),dim);
end
