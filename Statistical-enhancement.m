%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Roman Nelson                                                            %
% 23 September 2013                                                       %
%                                                                         %
% Image Enhancement Through Statistical Analysis                          %
% Program takes as an input an input along with a desired mean and        %
% standard deviation and gives as an output three different structures p, %
% q, and r, along with their respective statistics.                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p,q,r] = hw2(a, desiredmean, desiredstd)
p = original(a); %calls input image statistic function
q = meanfixer(a, desiredmean, desiredstd, p); %calls Brightness and Contrasst Mapping function
r = stdfixer(a, desiredmean, desiredstd, p); %calls Histogram Specification function
end

function p =  original(image) %function to create input image statistical structure
p.pix = image; %assigns image to variable in the structure
[p.M, p.N] = size(p.pix); %determines image size
p.minimum = min(p.pix(:)); %image minimum graysccale level
p.maximum = max(p.pix(:)); %image maximum grayscale level
p.range = p.maximum - p.minimum; %grayscale range
p.mean = mean(p.pix(:)); %image grayscale mean
p.var = var(double(p.pix(:)),1); %image grayscale variance
p.std = sqrt(p.var); %image grayscale standard deviation
p.msqr = meansquare(double(p.pix(:))); %image grayscale meansquare
p.rms = sqrt(p.msqr); %image grayscale root mean square
p.skewness = skewness(double(p.pix(:))); %image grayscale skewness
p.kurtosis = kurtosis(double(p.pix(:))); %image grayscale kurtosis
p.fdf = histc(double(p.pix(:)), p.minimum:p.maximum); %image grayscale fdf
p.rdf = p.fdf / (p.M * p.N); %image grayscale rdf
p.crdf = cumlativerdf(p.fdf,p.rdf, p.minimum, p.maximum, p.range); %image grayscale cumlative rdf
p.entropy = entropy(p.pix); % image entropy
end

function q = meanfixer(image, desiredmean, desiredstd, p) % function to create brightness and contrast mapped image statistical structure
q.pix = image; %assigns image to variable in the structure
[q.M, q.N] = size(uint8(q.pix)); %determines image size
q.pix = stretchy(desiredstd, p.std, desiredmean, p.mean, q.M, q.N, q.pix, p.maximum, p.minimum); % Brightness and Contrast Mapping 
q.var = var(double(q.pix(:)),1); %image grayscale variance
q.minimum = min(q.pix(:)); %image minimum graysccale level
q.maximum = max(q.pix(:)); %image maximum grayscale level
q.range = q.maximum - q.minimum; %grayscale range
q.mean = mean(q.pix(:)); %image grayscale mean
q.std = sqrt(q.var); %image grayscale standard deviation
q.msqr = meansquare(double(q.pix(:))); %image grayscale meansquare
q.rms = sqrt(q.msqr); %image grayscale root mean square
q.skewness = skewness(double(q.pix(:))); %image grayscale skewness
q.kurtosis = kurtosis(double(q.pix(:))); %image grayscale kurtosis
q.fdf = histc(double(uint8(q.pix(:))), 0:255); %image grayscale fdf
q.rdf = q.fdf / (q.M * q.N); %image grayscale rdf
q.crdf = cumlativerdf(q.fdf,q.rdf,q.minimum, q.maximum, q.range); %image grayscale cumlative rdf
q.entropy = entropy(q.pix); % image entropy
end

function r = stdfixer(image, desiredmean, desiredstd, p) %function to create a Histogram Specified image statistical structure
r.pix = image; %assigns image to variable in the structure
r.dstd = desiredstd; %assigns desired standard deviation to variable in the structure
r.dmean = desiredmean;%assigns desired mean to variable in the structure
[r.M, r.N] = size(uint8(r.pix)); %determines image size
r.desiredCRDF = [normcdf(1:255,r.dmean,r.dstd),1];%CRDFdesire (r.dstd, r.dmean, r.M, r.N, p.minimum, p.maximum); %creates the desired CRDF
%r.pix = histspec(p.crdf, r.M, r.N, r.pix, r.desiredCRDF, p.range); %creats and assigns histogram specified image to a structure variable
r.pix = histeq(uint8(r.pix), r.desiredCRDF);
r.minimum = min(uint8(r.pix(:))); %image minimum graysccale level
r.maximum = max(uint8(r.pix(:))); %image maximum grayscale level
r.range = r.maximum - r.minimum; %grayscale range
r.mean = mean(r.pix(:)); %image grayscale mean
r.var = var(double(r.pix(:)),1); %image grayscale variance
r.std = sqrt(r.var); % image grayscale variance
r.msqr = meansquare(double(r.pix(:))); %image grayscale meansquare
r.rms = sqrt(r.msqr); % image grayscale root mean square
r.skewness = skewness(double(r.pix(:))); % image grayscale skewness
r.kurtosis = kurtosis(double(r.pix(:))); % image grayscale kurtosis
r.fdf = histc(double(r.pix(:)), 0:255); %image grayscale fdf
r.rdf = r.fdf / (r.M * r.N); % image grayscale relative distrobution function
r.crdf = cumlativerdf(r.fdf,r.rdf,r.minimum, r.maximum, r.range); %image grayscale cumlative rdf
r.entropy = entropy(r.pix); %image grayscale entropy
end

function msqr = meansquare(input)
[M,N] = size(input);
varsum = 0;
for m = 1:M
    for n = 1:N
        x = double(input(m,n));
        varsum = varsum + x*x;
    end
end
msqr = varsum / (M * N);
end % function to determine the mean square

function crdf = cumlativerdf(fdf,rdf, minimum, maximum, range) %function to determine cumlative rdf
for indx = 2:(range + 1)
    cfdf(1) = fdf(1);
    crdf(1) = rdf(1);
    cfdf(indx) = cfdf(indx - 1) + fdf(indx);
    crdf(indx) = cfdf(indx) / (minimum * maximum);
end
end

function bob = stretchy(desiredstd, inputstd, desiredmean, inputmean, M, N, image, max, min) % function to perform contrast and brightness mapping
scale = desiredstd / double(inputstd);
offset = desiredmean - double(inputmean) * scale;
newimage = zeros(M,N);
for m = 1:M
    for n = 1:N
        tmpval = scale * (double(image(m,n)) - inputmean) + desiredmean;
        if(tmpval < 0)
            tmpval = 0;
        elseif(tmpval > 255)
            tmpval = 255;
        end
        newimage(m,n) = tmpval;
    end
    bob = uint8(newimage);
    
end
end

function dc = CRDFdesire (sigma, mu, M, N, min, max) % function to create the desired CRDF distribution
x = 0:255;
for i = 1:256
rdf(i) = double(exp(-1*(x(i) - mu)^2 / (2 * sigma^2)) / (sigma * sqrt(2*pi))); % Gaussian distribution
end
fdf = rdf * M * N;
dc = uint8(cumlativerdf(fdf, rdf, min, max, 255));
end

function bro = histspec(crdf, M, N, image, desiredCRDF, max) % function to perform histogram specification
prev = 0.0; %initializes prev variable
pos = 1; %initializes pos variable
hi = zeros(M,N);
for indx = 1:(max+1)
    mid = (prev + crdf(indx)) / 2.0; % finds mean between sum of adjacent crdf values
    while (desiredCRDF(pos) < mid) % while loop to determine graylevel (value of pos) for gmap
        pos = pos + 1; %increments graylevel until mid is less than the desiredCRDF at pos
    end
    gmap(indx) = pos - 1; %assigns a grayscale level to gmap coresponding to the desiredCRDF level
    prev = crdf(indx); %assigns next crdf level
end

for m = 1:M % for loops for assigning grayscale values from gmap to new image
    for n = 1:N
        indx = uint8(image(m,n)) +1;
        hi(m,n) = gmap(indx);
        
bro = uint8(hi);
    end
end
end


