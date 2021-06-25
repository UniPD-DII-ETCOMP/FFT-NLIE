%%
% This function has been inspired by 
% lse_matvect_mult.m which is available at
% https://github.com/acyucel/VoxHenry
%%
function [JOut_full] = multiplyMATVECT(JIn0,zloc,Aee)
% dx = d(1); dy = d(2); dz = d(3);
num_node=size(Aee,1); %num potential nodes, 
num_curr=size(Aee,2); %num current faces
% [LfN, MfN, NfN, ~] = size(CL);
% JIn = zeros(L, M, N, 3); %3 because we have 3 basis functions
% JOut = zeros(L, M, N, 3);    
% JIn(idxF) = JIn0(1:num_curr);
JOut_full = zeros(num_curr+num_node,1);
% fJ = fftn(JIn(:,:,:,1),[LfN, MfN, NfN]);
% Jout1 = CL(:,:,:,1) .* fJ; % Lxx*Jx
% fJ = fftn(JIn(:,:,:,2),[LfN, MfN, NfN]);
% Jout2 = CL(:,:,:,2) .* fJ; % Lyy*Jy
% fJ = fftn(JIn(:,:,:,3),[LfN, MfN, NfN]);
% Jout3 = CL(:,:,:,3) .* fJ; % Lzz*Jz
% Jout1 = ifftn(Jout1);
% Jout2 = ifftn(Jout2);
% Jout3 = ifftn(Jout3);
% JOut(:,:,:,1) = (dx/(dy*dz)) .* z_realx .* JIn(:,:,:,1);
% JOut(:,:,:,2) = (dy/(dx*dz)) .* z_realy .* JIn(:,:,:,2);
% JOut(:,:,:,3) = (dz/(dx*dy)) .* z_realz .* JIn(:,:,:,3);
JOut = zloc.*JIn0(1:num_curr);
JOut_full(1:num_curr) = JOut;
JOut_full(1:num_curr) = JOut_full(1:num_curr) + (Aee.'*JIn0(num_curr+1:num_curr+num_node)) ;
%
JOut_full(num_curr+1:end) =  Aee*JIn0(1:num_curr);
end

