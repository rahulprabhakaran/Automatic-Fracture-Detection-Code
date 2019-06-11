function coeffs = CSHRMsheardec(X,contShearletSystem)
% CSHRMsheardec Compute a 2D complex-valued shearlet transform.
% 
% Usage:
% 
%  coeffs = CSHRMsheardec(X,contShearletSystem)
% 
% Example:
% 
%  img = double(imread('lena.jpg'));
%  shearletSystem = CSHRMgetContEdgeSystem(size(img,1),size(img,2));
%  coeffs CSHRMsheardec(img,shearletSystem);
% 
% 
% See also: CSHRMgetContEdgeSystem, CSHRMgetContRidgeSystem

    X = fftshift(fft2(ifftshift(X)));
    coeffs = zeros(contShearletSystem.size(1),contShearletSystem.size(2),contShearletSystem.nShearlets);
    %coeffs = gpuArray(zeros(contShearletSystem.size(1),contShearletSystem.size(2),contShearletSystem.nShearlets));

    for k = 1:contShearletSystem.nShearlets
        coeffs(:,:,k) = conj(fftshift(ifft2(ifftshift(X.*conj(contShearletSystem.shearlets(:,:,k))))));
    end
end

%  Copyright (c) 2016. Rafael Reisenhofer
%
%  Part of CoShREM Toolbox v1.1
%  Built Mon, 11/01/2016
%  This is Copyrighted Material