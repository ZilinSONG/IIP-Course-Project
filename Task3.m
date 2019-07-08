function Task3


function [Lambda1,Lambda2,Ix,Iy]=eig2image(Dxx,Dxy,Dyy) % get eigen values
% This function eig2image calculates the eigen values from the
% hessian matrix, sorted by abs value. And gives the direction
% of the ridge (eigenvector smallest eigenvalue) .
% 
% [Lambda1,Lambda2,Ix,Iy]=eig2image(Dxx,Dxy,Dyy)
%
%
% | Dxx  Dxy |
% |          |
% | Dxy  Dyy |
% Compute the eigenvectors of J, v1 and v2
tmp = sqrt((Dxx - Dyy).^2 + 4*Dxy.^2);
v2x = 2*Dxy; v2y = Dyy - Dxx + tmp;

% Normalize
mag = sqrt(v2x.^2 + v2y.^2); i = (mag ~= 0);
v2x(i) = v2x(i)./mag(i);
v2y(i) = v2y(i)./mag(i);

% The eigenvectors are orthogonal
v1x = -v2y; 
v1y = v2x;

% Compute the eigenvalues
mu1 = 0.5*(Dxx + Dyy + tmp);
mu2 = 0.5*(Dxx + Dyy - tmp);

% Sort eigen values by absolute value abs(Lambda1)<abs(Lambda2)
check=abs(mu1)>abs(mu2);

Lambda1=mu1; Lambda1(check)=mu2(check);
Lambda2=mu2; Lambda2(check)=mu1(check);

Ix=v1x; Ix(check)=v2x(check);
Iy=v1y; Iy(check)=v2y(check);
end

function [Dxx,Dxy,Dyy] = Hessian2D(I,Sigma)   % get Hessian matrix
%  This function Hessian2 Filters the image with 2nd derivatives of a 
%  Gaussian with parameter Sigma.
% 
% [Dxx,Dxy,Dyy] = Hessian2(I,Sigma);
% 
% inputs,
%   I : The image, class preferable double or single
%   Sigma : The sigma of the gaussian kernel used
%
% outputs,
%   Dxx, Dxy, Dyy: The 2nd derivatives
%
% example,
%   I = im2double(imread('moon.tif'));
%   [Dxx,Dxy,Dyy] = Hessian2(I,2);
%   figure, imshow(Dxx,[]);
%
% Function is written by D.Kroon University of Twente (June 2009)

if nargin < 2, Sigma = 1; end

% Make kernel coordinates
[X,Y]   = ndgrid(-round(3*Sigma):round(3*Sigma));

% Build the gaussian 2nd derivatives filters
DGaussxx = 1/(2*pi*Sigma^4) * (X.^2/Sigma^2 - 1) .* exp(-(X.^2 + Y.^2)/(2*Sigma^2));
DGaussxy = 1/(2*pi*Sigma^6) * (X .* Y)           .* exp(-(X.^2 + Y.^2)/(2*Sigma^2));
DGaussyy = DGaussxx';

Dxx = imfilter(I,DGaussxx,'conv');
Dxy = imfilter(I,DGaussxy,'conv');
Dyy = imfilter(I,DGaussyy,'conv');
% Dxx = imfilter(I,DGaussxx);
% Dxy = imfilter(I,DGaussxy);
% Dyy = imfilter(I,DGaussyy);
end

function [outIm,whatScale,Direction] = FrangiFilter2D(I, options)   % Frangi2D filter

% This function FRANGIFILTER2D uses the eigenvectors of the Hessian to  
% compute the likeliness of an image region to vessels, according  
% to the method described by Frangi:2001 (Chapter 2).  
% [J,Scale,Direction] = FrangiFilter2D(I, Options)  
% inputs,  
%   I : The input image (vessel image)  
%   Options : Struct with input options,  
%       .FrangiScaleRange : The range of sigmas used, default [1 8]  
%       .FrangiScaleRatio : Step size between sigmas, default 2  
%       .FrangiBetaOne : Frangi correction constant, default 0.5  
%       .FrangiBetaTwo : Frangi correction constant, default 15  
%       .BlackWhite : Detect black ridges (default) set to true, for  
%                       white ridges set to false.  
%       .verbose : Show debug information, default true  
% outputs,  
%   J : The vessel enhanced image (pixel is the maximum found in all scales)  
%   Scale : Matrix with the scales on which the maximum intensity   
%           of every pixel is found  
%   Direction : Matrix with directions (angles) of pixels (from minor eigenvector)     
%  
% Written by Marc Schrijver, 2/11/2001  
% Re-Written by D.Kroon University of Twente (May 2009)  
  
defaultoptions = struct('FrangiScaleRange', [1 4], 'FrangiScaleRatio', 0.11, 'FrangiBetaOne', 3.75, 'FrangiBetaTwo', 4.5, 'verbose',false,'BlackWhite',true);  
  
% Process inputs  
if(~exist('options','var'))
    options=defaultoptions;   
else  
    tags = fieldnames(defaultoptions);  
    for i=1:length(tags)  
         if(~isfield(options,tags{i})),  options.(tags{i})=defaultoptions.(tags{i}); end  
    end  
    if(length(tags)~=length(fieldnames(options)))   
        warning('FrangiFilter2D:unknownoption','unknown options found');  
    end  
end  
  
  
sigmas=options.FrangiScaleRange(1):options.FrangiScaleRatio:options.FrangiScaleRange(2);  
sigmas = sort(sigmas, 'ascend');  
  
beta  = 2*options.FrangiBetaOne^2;  
c     = 2*options.FrangiBetaTwo^2;  
  
% Make matrices to store all filterd images  
ALLfiltered=zeros([size(I) length(sigmas)]);  
ALLangles=zeros([size(I) length(sigmas)]);  
  
% Frangi filter for all sigmas  
for i = 1:length(sigmas) 
    % Show progress  
    if(options.verbose)  
        disp(['Current Frangi Filter Sigma: ' num2str(sigmas(i)) ]);  
    end  
      
    % Make 2D hessian  
    [Dxx,Dxy,Dyy] = Hessian2D(I,sigmas(i));  
      
    % Correct for scale  
    Dxx = (sigmas(i)^2)*Dxx;  
    Dxy = (sigmas(i)^2)*Dxy;  
    Dyy = (sigmas(i)^2)*Dyy;  
     
    % Calculate (abs sorted) eigenvalues and vectors  
    [Lambda2,Lambda1,Ix,Iy]=eig2image(Dxx,Dxy,Dyy);  
  
    % Compute the direction of the minor eigenvector  
    angles = atan2(Ix,Iy);  
  
    % Compute some similarity measures  
    Lambda1(Lambda1==0) = eps;  
    Rb = (Lambda2./Lambda1).^2;  
    S2 = Lambda1.^2 + Lambda2.^2;  
     
    % Compute the output image  
    Ifiltered = exp(-Rb/beta) .*(ones(size(I))-exp(-S2/c));  

    if(options.BlackWhite)  
        Ifiltered(Lambda1<0)=0;  
    else  
        Ifiltered(Lambda1>0)=0;  
    end  
    % store the results in 3D matrices  
    ALLfiltered(:,:,i) = Ifiltered;  
    ALLangles(:,:,i) = angles;  
end  
  
% Return for every pixel the value of the scale(sigma) with the maximum   
% output pixel value  
if length(sigmas) > 1
    [outIm,whatScale] = max(ALLfiltered,[],3);  
    outIm = reshape(outIm,size(I)); 
    if(nargout>1)  
        whatScale = reshape(whatScale,size(I));
    end  
    if(nargout>2)  
        Direction = reshape(ALLangles((1:numel(I))'+(whatScale(:)-1)*numel(I)),size(I));  
    end  
else  
    outIm = reshape(ALLfiltered,size(I));  
    if(nargout>1)  
            whatScale = ones(size(I)); 
    end  
    if(nargout>2)  
        Direction = reshape(ALLangles,size(I));  
    end
end 
end

function task3main()
    img = imread('1.tif'); % read image file
    mask = imread('mask_1.tif'); % read mask file
    img_gray = rgb2gray(img); % Convert to grayscale
    gausFilter = fspecial('gaussian',[5,5],0.4); % gaussian filter template
    img_gray =imfilter(img_gray,gausFilter,'replicate'); % gaussian filter
    gausFilter = fspecial('gaussian',[7,7],0.66); % gaussian filter template
    img_gray =imfilter(img_gray,gausFilter,'replicate'); % gaussian filter
    img=double(img_gray);
    pn_lab =0; % the number of blood vessel pixels in label file
    pn_image = 0; % the number of blood vessel pixels in tested file 
    pn_lab_back = 0; % the number of background pixels in label file
    pn_img_back = 0; % the number of background pixels in tested file
    p_img=FrangiFilter2D(img); % use frangi filter
    p_img(mask == 0) = 0; % remove the edge
    p_img = imbinarize(p_img,0.055);
    se=strel('disk',2);
    p_img=imclose(p_img,se);
    label_image = imread('label_1.tif'); % read label file to calculate P N T
    [row,col] = size(label_image);
    for r = 1:row
        for c = 1:col
            if(label_image(r,c) == 255) % blood vessel pixels check (label)
               pn_lab = pn_lab + 1;
            end
            if(mask(r,c) == 255 && label_image(r,c) == 0 )
                pn_lab_back = pn_lab_back + 1; % background pixels check (label)
            end
            if(p_img(r,c) > 0 && label_image(r,c) == 255 ) % check if is correct vessel pixels
                pn_image = pn_image +1;
                %p_img(r,c) = 1;
            end
            if(p_img(r,c) ==0  && mask(r,c) == 255 && label_image(r,c) == 0) %check if is correct background pixels
                pn_img_back = pn_img_back + 1;
            end
        end
    end
    P = pn_image/pn_lab %  percentage of blood vessel pixels that is being correctly classified as blood vessel
    N = pn_img_back/pn_lab_back % The percentage of background pixels that is being correctly classified as background
    T = (pn_image + pn_img_back)/(pn_lab+pn_lab_back) % The percentages ofpixels are beingcorrectly classified
    figure('name','Final result');
    imshow(p_img);
end

task3main();

end