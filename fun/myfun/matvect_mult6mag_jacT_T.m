%% System MVM with vector

%**************************************************************************
function [JOut_full] = matvect_mult6mag_jacT_T(JIn0,CP,z_realx,z_realy,z_realz,idxF,...
                       idxV,d,Aee,L,M,N,AeeR,idxVR,H)

dx = d(1); dy = d(2); dz = d(3);

%%Physical 

%%Dimensions of Problem
num_node=size(Aee,1); %num potential nodes, 
num_nodeR=size(AeeR,1); %num potential nodes, 
num_curr=size(Aee,2); %num current faces

%Circulant fft dimensions (same for CL, CP)
[LfN, MfN, NfN, ~] = size(CP);
%Grid dimensions
% [L, M, N] = size(z_real);

%%Allocate space
JIn = zeros(L, M, N, 3); %3 because we have 3 basis functions
JOut = zeros(L, M, N, 3);    
% translate from local (idx) to global (L,M,N) coordinates
JIn(idxF) = JIn0(1:num_curr);

%Initialize output vector
JOut_full = zeros(num_curr+num_nodeR,1);

%%1. Do (R)*Jin for Jin relative to currents and add to current outputs
% %iwL*Jin via FFT fast procedure
% %x component of JIn, store contribution on 3 components of Jout
% fJ = fftn(JIn(:,:,:,1),[LfN, MfN, NfN]);
% Jout1 = CL(:,:,:,1) .* fJ; % Lxx*Jx
% 
% %y component of JIn, add contribution on 3 components of Jout
% fJ = fftn(JIn(:,:,:,2),[LfN, MfN, NfN]);
% Jout2 = CL(:,:,:,2) .* fJ; % Lyy*Jy
% 
% %z component of JIn, store contribution on 2 components of Jout
% fJ = fftn(JIn(:,:,:,3),[LfN, MfN, NfN]);
% Jout3 = CL(:,:,:,3) .* fJ; % Lzz*Jz

%%Apply ifft 
% Jout1 = ifftn(Jout1);
% Jout2 = ifftn(Jout2);
% Jout3 = ifftn(Jout3);

%%Add R matrix
JOut(:,:,:,1) = (dx/(dy*dz)) .* z_realx .* JIn(:,:,:,1);% + Jout1(1:L,1:M,1:N);
JOut(:,:,:,2) = (dy/(dx*dz)) .* z_realy .* JIn(:,:,:,2);% + Jout2(1:L,1:M,1:N);
JOut(:,:,:,3) = (dz/(dx*dy)) .* z_realz .* JIn(:,:,:,3);% + Jout3(1:L,1:M,1:N);

%%Return local coordinates related to material positions
JOut = JOut(idxF);
JOut_full(1:num_curr) = JOut;
JOut_full(1:num_curr) = JOut_full(1:num_curr) + H'*JIn0(1:num_curr);


%%3. Do P*A*Jin for Jin relative to currents and add to potential nodes
%Intermediate step: obtain q = A*Jin
noq = JIn0(num_curr+1:end);
%Now multiplication with circulant via FFT
QIn = zeros(L,M,N);  
QIn(idxVR) = noq;
fJ = fftn(QIn(:,:,:),[LfN, MfN, NfN]);
Jout = CP(:,:,:) .* fJ; % P*Q
%%Apply ifft 
JOut = ifftn(Jout);
JOut = JOut(1:L,1:M,1:N);
%%Return local coordinates related to material positions
JOut = JOut(idxVR);

%%4. Do -Jin relative to potential nodes and add to potential nodes
JOut_full(1:num_curr) = JOut_full(1:num_curr) + (AeeR.'*JOut) ;

JOut_full(num_curr+1:end) = ...
    (AeeR*JIn0(1:num_curr) - (JIn0(num_curr+1:end)));

JOut_full=JOut_full.';

end

