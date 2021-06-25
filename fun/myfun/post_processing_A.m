%% Post Processing A
%************************************************************************
%Function computing A via FFT product of Circulant tensor
%INPUT: - vsol: currents [A]
%       -xyz: 4D coordinates
%       -d: array of grid steps
%       -Integration: type of integral (NumNum or NumAn)
%       -circulant: circulant tensor (if circulant = [] reconstruct inside function)
%
%************************************************************************
%function [A,a_dof] = post_processing_A(vsol,CL,d,xyz)
function [A,a_dof] = post_processing_A(vsol,xyz,d,circulant,Ae1x,Ae1y,Ae1z)

dx=d(1); dy=d(2); dz=d(3);
[L, M, N, ~] = size(vsol);

CL = circulant;
%%Compute Circulant-Vector product via FFT
JIn = vsol;
[LfN, MfN, NfN, ~] = size(CL);


%
idxVe=ones(L,M,N);
idxVe=find(idxVe);
indPotv=[];
valPotv=[];

% %%Computing Incidence matrix
% disp(' Start computing incidence');
% mytic_incid = tic;
% [~,~,~,~,~,~,Ae1x,Ae1y,Ae1z] = incidence_matrix3(L*M*N,[L M N],idxVe);
% mytoc_incid = toc(mytic_incid);
% disp([' Time for computing incidence ::: ' ,num2str(mytoc_incid)]);

if dx == dy && dy == dz
    %%Apply fft and mv-op for each of the components of JIn
    %x component of JIn, store contribution on 3 components of Jout
    fJ = fftn(JIn(:,:,:,1),[LfN, MfN, NfN]);
    Jout1 = (CL(:,:,:,1)) .* fJ; % Gxx*Jx

    %y component of JIn, add contribution on 3 components of Jout
    fJ = fftn(JIn(:,:,:,2),[LfN, MfN, NfN]);
    Jout2 = (CL(:,:,:,1)) .* fJ; % Gyy*Jy

    %z component of JIn, store contribution on 2 components of Jout
    fJ = fftn(JIn(:,:,:,3),[LfN, MfN, NfN]);
    Jout3 = (CL(:,:,:,1)) .* fJ; % Gzz*Jz
else
    %%Apply fft and mv-op for each of the components of JIn
    %x component of JIn, store contribution on 3 components of Jout
    fJ = fftn(JIn(:,:,:,1),[LfN, MfN, NfN]);
    Jout1 = (CL(:,:,:,1)) .* fJ; % Gxx*Jx

    %y component of JIn, add contribution on 3 components of Jout
    fJ = fftn(JIn(:,:,:,2),[LfN, MfN, NfN]);
    Jout2 = (CL(:,:,:,2)) .* fJ; % Gyy*Jy

    %z component of JIn, store contribution on 2 components of Jout
    fJ = fftn(JIn(:,:,:,3),[LfN, MfN, NfN]);
    Jout3 = (CL(:,:,:,3)) .* fJ; % Gzz*Jz
    
end
%%Apply ifft 
Jout1 = ifftn(Jout1);
Jout2 = ifftn(Jout2);
Jout3 = ifftn(Jout3);

%%Extrapolate voxels
a_dof(:,:,:,1) = Jout1(1:L,1:M,1:N);
a_dof(:,:,:,2) = Jout2(1:L,1:M,1:N);
a_dof(:,:,:,3) = Jout3(1:L,1:M,1:N);

%%Return local coordinates related to material positions
MM(1:L*M*N,1)          =1/(dx/(dy*dz));
MM(L*M*N+1:2*L*M*N,1)  =1/(dy/(dx*dz)) ;
MM(2*L*M*N+1:3*L*M*N,1)=1/(dz/(dx*dy));
a=MM.*reshape(a_dof,3*L*M*N,1);

%%Post Processing
A = reshape(a,L, M, N, 3); %3 because we have 3 basis functions
[A,~] = fun_my_postRT2(A,L*M*N,Ae1x,Ae1y,Ae1z,xyz,L,M,N,d);


end