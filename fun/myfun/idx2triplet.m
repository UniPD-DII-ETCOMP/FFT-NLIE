%% Index to Triplet
%************************************************************************
%%This constructs map i -> [l m n] 
%in a regular 3D grid with N points for each dimension

%Input: -i index
%       -LMN = [L M N] array of number of points for each direction
%Output: - lmn = [l m n] triplet
%************************************************************************

function [lmn] = idx2triplet(i,L,M,N)

lmn(3) = floor((i-1)/(L*M)); %
lmn(2) = floor((i-1-lmn(3)*L*M)/L); %
%lmn(1) = floor((i-1-lmn(3)*L*M)/L); %
ii = (i-1)-lmn(3)*L*M; %
lmn(1) = ii - floor(ii/L)*L;
lmn = lmn + 1;

end
