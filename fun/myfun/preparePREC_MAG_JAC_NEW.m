%% Sparse Preconditioner Preparation: Used For Capacitive Effects
%*************************************************************************
%*************************************************************************
function [Y_inv,P_diag,D_diag,LL,UU,PP,QQ,RR] = preparePREC_MAG_JAC_NEW(d,z_realF,idxFx,idxFy,...
                    idxFz,st_sparse_preconP,AeeR,Aee,Kt,H,prec_mod,fl_precon_type)


dx = d(1); dy = d(2); dz = d(3);




%%Number of unknowns
num_node = size(Aee,1); % tutti, anche quelli a potenziale vincolato
num_nodeR = size(AeeR,1); % senza quelli a potenziale vincolato
num_curr = size(Aee,2);
num_curr_one3rd = num_curr/3; %since we have 3 basis functions for current
num_currx=length(idxFx);
num_curry=length(idxFy);
num_currz=length(idxFz);

%%Compute Inverse of Z matrix, Only Diagonal
%Actually this is Y matrix in VoxHenry Article [1]
tic
% get the actual values of 'OneoverSigma_e'. This allows different conductivites for each voxel.
diag_pulse=zeros(num_curr,1);
diag_pulse(1:num_currx,1)                                        =1./(z_realF(idxFx)*dx/(dy*dz) );
diag_pulse(num_currx+1:num_currx+num_curry,1)                    =1./(z_realF(Kt+idxFy)*dy/(dz*dx) );
diag_pulse(num_currx+num_curry+1:num_currx+num_curry+num_currz,1)=1./(z_realF(Kt+Kt+idxFz)*dz/(dx*dy) );

%%Sparse 'Y_inv' formation
inds=zeros(num_curr,3);
% rows for sparse 'Y_inv' formation
inds(1:num_curr,1)=[1:1:num_curr];
% columns for sparse 'Y_inv' formation
inds(1:num_curr,2)=inds(1:num_curr,1);
% values for sparse 'Y_inv' formation, split in three sets:
% first set is from 1 to 3/5 of num_curr, i.e. Ix, Iy, Iz
if strcmp(prec_mod,'diag')
    inds(:,3)=1./(1./(diag_pulse)+diag(H));
    % inds(1:num_curr_one3rd,3)=abs(diag_pulse(1:num_curr_one3rd));
    % inds(num_curr_one3rd+1:2*num_curr_one3rd,3)=abs(diag_pulse);
    % inds(2*num_curr_one3rd+1:3*num_curr_one3rd,3)=abs(diag_pulse);
    %%Create Sparse 'Y_inv'
    Y_inv=sparse(inds(:,1),inds(:,2),inds(:,3));
elseif strcmp(prec_mod,'allH')
    error('to do')
    inds(:,3)=(diag_pulse);
    % inds(1:num_curr_one3rd,3)=abs(diag_pulse(1:num_curr_one3rd));
    % inds(num_curr_one3rd+1:2*num_curr_one3rd,3)=abs(diag_pulse);
    % inds(2*num_curr_one3rd+1:3*num_curr_one3rd,3)=abs(diag_pulse);
    %%Create Sparse 'Y_inv'
%     Y_inv=sparse(inds(:,1),inds(:,2),1./(1./inds(:,3)+(H))); % this sucks, may use a sparse inv    
end
%Memory requirements

%%Compute Diagonal of P matrix
%%Schur complement in [2]Â eq. 2.2 is S = D - Cinv(A)B
%inv(A) is computed above, using only diagonal. C = P*A_r in our problem. 
%QUESTION: can we use only diagonal of P in this product?
tic 
%%Entries of diagonal 'P_diag'
diag_pulse = zeros(num_nodeR,1); 
diag_pulse(1:num_nodeR,1) = st_sparse_preconP(1); %diagonal of P

%%Sparse 'P_diag' formation
inds = zeros(num_nodeR,3);
inds(1:num_nodeR,1) = [1:1:num_nodeR];
inds(1:num_nodeR,2)=inds(1:num_nodeR,1);
inds(:,3)=(diag_pulse);

%%Create Sparse 'P_diag' 
P_diag =sparse(inds(:,1),inds(:,2),inds(:,3));
% global Ptot
% P_diag=Ptot;
%Memory requirements

%%Sparse Diagonal of D matrix: D = -iwI
tic 
%%Entries of diagonal 'D_diag'
diag_pulse = zeros(num_nodeR,1); 
diag_pulse(1:num_nodeR,1) = -1; %diagonal of P

%%Sparse 'D_diag' formation
inds = zeros(num_nodeR,3);
inds(1:num_nodeR,1) = [1:1:num_nodeR];
inds(1:num_nodeR,2)=inds(1:num_nodeR,1);
inds(:,3)= diag_pulse;

%%Create Sparse 'D_diag' 
D_diag =sparse(inds(:,1),inds(:,2),inds(:,3));
%Memory requirements

%%Obtain Schur complement
tic
    
% as the rows of Ae corresponding to ground or excitation nodes have been removed,
% obtaining what in the VoxHenry paper is called 'Ar', we can proceed with the standard
% Schur complement from the system (but in the code 'Ar' is actually 'Ae' without the zeroed rows)
% [A  B] = [Z       Ae']
% [C  D]   [P*Ae  -iwI ]

% symm. voltage source and current source will use this
Sch_comp= D_diag - P_diag*(AeeR*Y_inv*AeeR.');


            
%%Computing Schur Inversion
[LL,UU,PP,QQ,RR] = lu(Sch_comp);


end