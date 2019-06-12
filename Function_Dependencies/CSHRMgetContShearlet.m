function [shearlet] = CSHRMgetContShearlet(rows,cols,waveletEffSupp,gaussianEffSupp,scalesPerOctave,shearLevel,alpha,sampleWaveletOffOrigin,scale,ori)
% CSHRMgetContShearlet Construct a 2D real-valued even-symmetric shearlet in the frequency domain.
% 
% Usage:
% 
%  shearlet = CSHRMgetContShearlet(rows,cols,waveletEffSupp,gaussianEffSupp,scalesPerOctave,shearLevel,alpha,sampleWaveletOffOrigin,scale,ori)
% 
% Example:
% 
% shearlet = CSHRMgetContShearlet(512,512,70,30,2,3,0.5,0,2,3);
% figure;
% subplot(1,2,1);
% imagesc(shearlet);axis tight;axis equal;axis off;
% title('Frequncy Domain');
% subplot(1,2,2);
% imagesc(fftshift(ifft2(ifftshift(shearlet))));axis tight;axis equal;axis off;
% title('Time Domain');
% 
% See also: CSHRMgetContEdgeSystem, CSHRMgetContRidgeSystem

    [~,coneh,ks] = CSHRMgetConeOris(shearLevel);
    isConeh = ismember(ori,coneh);
    if isConeh
        omegaWav = (63*double(waveletEffSupp)/512)*yapuls(rows);
        omegaGau = (74*double(gaussianEffSupp)/512)*yapuls(cols*(2^(shearLevel-2)));    
    else
        omegaWav = (63*double(waveletEffSupp)/512)*yapuls(cols);
        omegaGau = (74*double(gaussianEffSupp)/512)*yapuls(rows*(2^(shearLevel-2)));
    end

    %scaling omega
    omegaGau = omegaGau/(2^((scale-1)/scalesPerOctave))^(alpha);
    omegaWav = omegaWav/(2^((scale-1)/scalesPerOctave));

    %mexican hat wavelet in the frequncy domain
    wavFreq = 2*pi*(omegaWav.^2).*exp(-omegaWav.^2 / 2 );
    wavDimOdd = mod(length(wavFreq),2);
    if sampleWaveletOffOrigin == 1
        wavFreq = SLpadArray(fftshift(wavFreq),[1,length(wavFreq)*2]);
        wavTime = real(fftshift(ifft(ifftshift(wavFreq))));
        if wavDimOdd
            if isConeh
                wavFreq = fft(ifftshift(fliplr(wavTime(1:2:(end-1)))));
            else
                wavFreq = fft(ifftshift(wavTime(1:2:(end-1))));
            end
        else
            if isConeh
                wavFreq = fft(ifftshift(fliplr(wavTime(2:2:end))));
            else
                wavFreq = fft(ifftshift(wavTime(2:2:end)));
            end
        end
    end
    gauFreq = exp(-omegaGau.^2 / 2 );

    if isConeh
        shearlet = fftshift(wavFreq'*gauFreq);
        shearlet = SLdshear(shearlet,-ks(ori),2);
        shearlet = shearlet(:,1:(2^(shearLevel-2)):end);
    else
        shearlet = fftshift((gauFreq)'*wavFreq);
        shearlet = SLdshear(shearlet,-ks(ori),1);
        shearlet = shearlet(1:(2^(shearLevel-2)):end,:);
    end
end


%  Copyright (c) 2016. Rafael Reisenhofer
%
%  Part of CoShREM Toolbox v1.1
%  Built Mon, 11/01/2016
%  This is Copyrighted Material