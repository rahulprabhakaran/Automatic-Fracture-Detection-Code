function contShearletSystem = CSHRMgetContRidgeSystem(rows,cols,waveletEffSupp,gaussianEffSupp,scalesPerOctave,shearLevel,alpha,scales)
% CSHRMgetContRidgeSystem Construct a system of complex-valued shearlets to detect ridges in a 2D grayscale image.
% 
% Usage (optional parameters are enclosed in angle brackets):
% 
%  contShearletSystem = CSHRMgetContRidgeSystem(rows,cols,<waveletEffSupp>,<gaussianEffSupp>,<scalesPerOctave>,<shearLevel>,<alpha>,<scales>)
% 
% Example:
% 
%  img = double(imread('barbara.jpg'));
%  shearletSystem = CSHRMgetContRidgeSystem(size(img,1),size(img,2));
%  ridges = CSHRMgetRidges(img,shearletSystem);
%  imshow(CSHRMgetOverlay(img,ridges));
% 
% See also: CSHRMgetRidges

    if (nargin < 3), waveletEffSupp = min([rows,cols])/7; end
    if (nargin < 4), gaussianEffSupp = min([rows,cols])/20; end
    if (nargin < 5), scalesPerOctave = 2; end
    if (nargin < 6), shearLevel = 3; end
    if (nargin < 7), alpha = 0.5; end
    if (nargin < 8), scales = 1:(scalesPerOctave*3.5); end
   
    nOris = 2^shearLevel+2;
    nShearlets = length(scales)*nOris*2;
    shearlets = zeros(rows,cols,nShearlets);
    [~,coneh,~] = CSHRMgetConeOris(shearLevel);



    for j = 1:length(scales)
        scale = scales(j);
        for offOrigin = [1,0]
             for ori = 1:nOris        
                shearlet = CSHRMgetContShearlet(rows,cols,waveletEffSupp,gaussianEffSupp,scalesPerOctave,shearLevel,alpha,offOrigin,scale,ori);
                shearlet = real(fftshift(ifft2(ifftshift(shearlet))));        

                if ismember(ori,coneh)
                    shearlet = hilbert(shearlet);
                else
                    shearlet = hilbert(shearlet')';
                end
                if ori == 1
                    normalizationFactor = sum(sum(abs(real(shearlet))));

                    shearlet = shearlet/normalizationFactor;

                    centerRidge = real(shearlet(floor(end/2)+1,(floor(end/2)+1):end));
                    positiveWidths(2*(j-1)+1 + offOrigin) = (find(centerRidge < 0,1,'first')-2+offOrigin)*2+(1-offOrigin);
                    for k = offOrigin:((max(positiveWidths((1+offOrigin):2:(end-1+offOrigin)))-1+offOrigin)/2)
                        expectedCoeffs(2*k+(1-offOrigin),j) = sum(sum(real(shearlet(:,(floor(end/2)+1-k):(floor(end/2)+(1-offOrigin)+k)))));
                    end
                else
                    shearlet = shearlet/normalizationFactor;
                end   
                shearlets(:,:,2*nOris*(j-1)+2*(ori-1)+1+offOrigin) = fftshift(fft2(ifftshift(shearlet)));
            end
        end
    end
    contShearletSystem = struct('shearlets',shearlets,'size',[rows,cols],'waveletEffSupp',waveletEffSupp,'gaussianEffSupp',gaussianEffSupp,'scalesPerOctave',scalesPerOctave,'shearLevel',shearLevel,'scales',scales,'nShearlets',nShearlets,'alpha',alpha,'detectRidges',1,'positiveWidths',positiveWidths, 'expectedCoeffs',expectedCoeffs);
end

%  Copyright (c) 2016. Rafael Reisenhofer
%
%  Part of CoShREM Toolbox v1.1
%  Built Mon, 11/01/2016
%  This is Copyrighted Material