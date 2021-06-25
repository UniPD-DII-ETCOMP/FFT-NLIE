%% Post Processing H
%************************************************************************
%Function computing post processing on scattered magnetic field H
%************************************************************************
function [H,h_dof] = post_processing_H(Jout,xyz,d,Ae1x,Ae1y,Ae1z,CK)

%%Grid dimensions
[L,M,N,~] = size(xyz);
dx = d(1); dy = d(2); dz = d(3);

idxVe=ones(L,M,N);
idxVe=find(idxVe);
indPotv=[];
valPotv=[];

%%Computing Incidence matrix

%%Compute Circulant-Vector product via FFT
JIn = Jout; %3 because we have 3 basis functions

[LfN, MfN, NfN, ~] = size(CK);

% x component of JIn, store contribution on 3 components of Jout
fJ = fftn(JIn(:,:,:,1),[LfN, MfN, NfN]); 
% Jout1 = 0 .* fJ; % First component Vout: 0 * Vin_x
Jout2 = -CK(:,:,:,3) .* fJ; % Second component Vout: -fK_z * Vin_x
Jout3 = CK(:,:,:,2) .* fJ; % Third component Vout: +fK_y * Vin_x

% y component of JIn, add contribution on 3 components of Jout
fJ = fftn(JIn(:,:,:,2),[LfN, MfN, NfN]); 
Jout1 = CK(:,:,:,3) .* fJ;  % First component Vout: +fK_z * Vin_y
% Jout2 = Jout2 + 0 .* fJ; % Second component Vout: 0 * Vin_y
Jout3 = Jout3 - CK(:,:,:,1) .* fJ; % Third component Vout: -fK_x * Vin_y

% z component of JIn, add contribution on 3 components of Jout
fJ = fftn(JIn(:,:,:,3),[LfN, MfN, NfN]); 
Jout1 = Jout1 - CK(:,:,:,2) .* fJ;  % First component Vout: -fK_y * Vin_z
Jout2 = Jout2 + CK(:,:,:,1) .* fJ; % Second component Vout: +fK_x * Vin_z
% Jout3 = Jout3 + 0 .* fJ; % Third component Vout: 0 * Vin_z

%%Apply ifft 
Jout1 = ifftn(Jout1);
Jout2 = ifftn(Jout2);
Jout3 = ifftn(Jout3);

%%Extrapolate voxels
h_dof(:,:,:,1) = Jout1(1:L,1:M,1:N);
h_dof(:,:,:,2) = Jout2(1:L,1:M,1:N);
h_dof(:,:,:,3) = Jout3(1:L,1:M,1:N);

%%Return local coordinates related to material positions
MM(1:L*M*N,1)          =1/(dx/(dy*dz));
MM(L*M*N+1:2*L*M*N,1)  =1/(dy/(dx*dz)) ;
MM(2*L*M*N+1:3*L*M*N,1)=1/(dz/(dx*dy));
h=MM.*reshape(h_dof,3*L*M*N,1);

%%Post Processing
H = reshape(h,L, M, N, 3); %3 because we have 3 basis functions
[H,~] = fun_my_postRT2(H,L*M*N,Ae1x,Ae1y,Ae1z,xyz,L,M,N,d);

end