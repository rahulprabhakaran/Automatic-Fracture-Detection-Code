function out = SLdshear(in,k,axis)

    if k == 0
        out = in;
        return;
    end
    rows = size(in,1);
    cols = size(in,2);

    out = zeros(size(in));
    if axis == 1
        for col = 1:cols
            out(:,col) = circshift(in(:,col),[k * (floor(cols/2)+1-col) 0]);
        end
    else        
        for row = 1:rows
            out(row,:) = circshift(in(row,:),[0 k * (floor(rows/2)+1-row)]);
        end
    end
end

%
%  Copyright (c) 2014. Rafael Reisenhofer
%
%  Part of ShearLab3D v1.1
%  Built Mon, 10/11/2014
%  This is Copyrighted Material
