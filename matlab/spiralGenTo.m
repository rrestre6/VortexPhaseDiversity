function [spiralMask] = spiralGenTo(amp, topCharge)
% spiralGen sumary
%   Generete spiral phase masks for integer
%   topological charge.
%   This code is used to replace another kind of 
%   spiral phase masks that have a vertex as 
%   center instead of a pixel. For topological 
%   charge ==1 we divide the domain into two 
%   regions, the left and right hand side of the cartesian 
%   plane and the. We assume that the center 
%   of the   matrix is the center of the cartesian 
%   plane, this means that:
%      * when m = m/2, y = 0
%      * when n = n/2, x = 0
%IN:
%   amp:            nxm Matrix with amplitude values. 
%   topCharge:   integer, represents a number of times
%                       where the phase changes from 0 to 2pi
%                   
%OUT: 
%   spiralMask
%
% Santiago Echeverri, Camilo Cano, René Restrepo, Daniel Sierra.

[m, n] = size(amp);

spiralMask = ones(m,n);

%if topCharge==0
%   spiralMask = spiralMask;
%end

    for Row = 1:m
        for Column = 1:(n/2)+1
            y = (Row-1 - m/2);
            x = (Column-1 - n/2);
            spiralMask(Row,Column) = atan((y)/(x)) + pi/2;
        end
    end
    %REGION 2
    for Row = 1:m
        for Column = n/2+1:n
            y = (Row-1-m/2);
            x = (Column-1-n/2);
            spiralMask(Row, Column) = atan((y)/(x))+pi+pi/2;
        end
    end
    spiralMask=spiralMask*topCharge;

    spiralMask(m/2+1,n/2+1) = (spiralMask(m/2,n/2) + spiralMask(m/2 + 1,n/2) + spiralMask(m/2 + 2,n/2) + spiralMask(m/2 ,n/2 +1) + spiralMask(m/2 + 2,n/2 + 1) + spiralMask(m/2 ,n/2 + 2) + spiralMask(m/2 + 1,n/2 + 2) + spiralMask(m/2 +2,n/2 +2))/8;%+0.5i;

    for Row = 1:m
        for Column = 1:n
            spiralMask(Row,Column) = exp(1i*spiralMask(Row,Column));

        end
    end

end


