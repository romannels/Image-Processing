%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------Joint Image Statistics----------------------------%
%                                                                         %
% Program that finds subimages within a larger image using autocorrelation%
% Part of program taken from Dr. Stephen Reichenbach's book:              %
% Digital Image Processing:Foundations, Models, and Programs              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [ r ] = hw4(p)
t = 0.9;
e = p(46:60,43:55); % "e" subimage extracted from "p"
f = p(40:60,202:210); % "f" subimage extracted from "p"
l = p(42:59,395:398); % "l" subimage extracted from "p"
[ey,ex] = size(e); %size of "e" image
e =(double(e) - double((mean(e(:)) * ones(ey,ex)))) / double(std(double(e(:))));
%normalize "e" subimage mean and standard deviation
[fy,fx] = size(f); %size of "f" image
f =(double( f) - double((mean(f(:)) * ones(fy,fx)))) / double(std(double(f(:))));
%normalize "f" subimage mean and standard deviation
[ly,lx] = size(l); %isze of "l" image
l =(double(l) - double((mean(l(:)) * ones(ly,lx)))) / double(std(double(l(:))));
%normalize "l" subimage mean and standard deviation
l = flipdim(l,2); %flips image "l" to correct anti aliasing
[eout,nume] = foundit(p,e,t); %finds occurances of "e" in image "p"
[lout,numl] = foundit(p,l,t); %finds occurances of "l" in image "p"
[fout,numf] = foundit(p,f,t); %finds occurances of "f" in image "p"
%imshow(lout);
howbige = size(nume); %size of "e" row/index matrix
howbigf = size(numf); %size of "f" row/index matrix
howbigl = size(numl); %size of "l" row/index matrix
%largeness = [howbige,howbigf,howbigl];
%overallsize = max(largeness);
%r = zeros(overallsize,3);
% if (howbige <overallsize)
%     for i = howbige + 1:overallsize
%         nume(i) = 0;
%     end
% end
% if (howbigl <overallsize)
%     for i = howbigl + 1:overallsize
%         numl(i) = 0;
%     end
% end
% if (howbigf <overallsize)
%     for i = howbigf + 1:overallsize
%         numf(i) = 0;
%     end
% end
%r = [101,nume;102,numf;108,numl];
weird = [nume(1,1:howbige(2)),numl(1,1:howbigl(2)),numf(1,1:howbigf(2))];
%matrix to hold all column indices

bigness = howbige+howbigf+howbigl; %total number of letters found (not 
% found contributes one to count)
r = zeros(bigness(2),3); %intialize output matrix
sumf = 0; %initialize sumf counter
suml = 0; %initialize suml counter
sume = 0; %initialize sume counter
for i = 1:bigness(2)
    if (i <= howbige(2)) %all "e" characters are first in weird matrix
        r(i,1) = 101; %ASCII "e" character
        sume = sume+1; %counter
        if (howbige > 1) %if letters are found
            r(i,2) = nume(2,sume); %set row index
        else
            r(i,2) = -1; %negative one for does not exist
        end
    end
    if ((i >= (bigness(2)))) %all "f" characters are last in weird matrix
        r(i,1) = 108; %set ASCII "f" character
        sumf = sumf+1; %counter
        if (howbigf > 1) %if letters are found
        r(i,2) = numf(2,sumf); %set row index
        else
            r(i,2) = -1; % negative one for does not exist
        end
        %disp(i);
    end
    if ((i > howbige(2)) && i < (bigness(2))) % all "l" characters are in
        %the middle of the weird matrix
    r(i,1) = 102; %ASCII "l" character
    suml = suml+1; %counter
        if(howbigl > 1) %if letters are found
         r(i,2) = numl(2,suml); %set row index
        else
            r(i,2) = -1; %negative one for does not exist
        end
   % disp(i);
    end 
    r(i,3) = weird(i); %displays column index
end
%disp(weird);
%disp(bigness);
%disp(howbige(2));
%disp(howbigf);
%disp(howbigl);
%disp(nume);
%disp(numl);
%disp(numf);
end

function [out,numletters] = foundit(p,f,t)
%much of this code was created by Dr. Reichenbach

[pM, pN] = size(p); %determine bounds of input image
[fM, fN] = size(f); %determine bounds of subimage
fsz = fM*fN; %determine amount of pixels in subimage
z = 0; %initialize output to zero
out = uint8(zeros(pM,pN)); %for debugging/remanant from Dr. R code
    sum = 0; %initizlize sum to zero
    for pn = 1:666 - fN +1 %search last sentence
        for pm = 125:159 - fM + 1 %search last sentence
        psum = 0; % For the mean of each region in p
        pssq = 0; % For the std of for each region in p
        prod = 0; % For the correlation with f of each region in p
            % Loop (fm,fn) to index each pixel of region at (pm,pn)
            for fn=1:fN %find autocorrelation in x dimension
             idx = pn+fn-1;%create index in x direction
                for fm=1:fM %find autocorrelation in y dimension
                 val = double(p(pm+fm-1,idx)); %retrieve value from input
                 %picture, index in y direction is created inside this
                 %keeps the search area inside the subimage bounds
                 psum = psum+val; %sum together all values in 
                 %input image within subimage dimensions
                 pssq = pssq+val*val; %Add the sqaure of the sums
                 prod = prod+val*f(fm,fn);
                end;
            end;
            % Test the cross correlation; set mask if match is good enough
            denom = fsz*pssq-psum*psum;
            if (denom > 0)
                crossCorrelation = prod / sqrt(denom);
                if (crossCorrelation >= t)
                    out(pm:pm+fM-1,pn:pn+fN-1) = p(pm:pm+fM-1,pn:pn+fN-1);
                    if(sum ==0)
                    sum = 1;
                    z(1,sum) = pn;  
                    z(2,sum) = pm;
                    %disp(z);
                    end
                    if(sum ~= 0 && z(1,sum) < pn - 10)
                    sum = 1 + sum;
                    z(1,sum) = pn;  
                    z(2,sum) = pm;
                   % disp(z); 
                    end
                end;
            end;
        end;
    end;
    numletters = z;
    %disp('numletters:');
    %disp(numletters);
    imshow(out);
   % figure;
end
%extracted image
%e - x: 43 - 55, y: 46 - 60
%f - x:203 - 210, y: 41 - 59
%l - x:395 - 397, y: 42 - 59
%line y:5 = 129:146  ,x: 1:666
