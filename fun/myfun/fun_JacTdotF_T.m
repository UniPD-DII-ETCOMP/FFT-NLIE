function [res] = fun_JacTdotF_T(x,d,idxV,idxF,idxFx,idxFy,...
                            idxFz,Kt,Ae,Ae1x,Ae1y,Ae1z,L,M,N,...
                            fun_J_to_coeff,...
                            opCirculantP_all,Aee,AeeR,idxVR,F,...
                            fun_J2_to_dcoeff,K)
%%
num_curr=length(idxF);
%% 
Jout = zeros(L,M,N,3);
Jout(idxF) = x(1:num_curr) ; % return to global variables
[~,normJ2] = fun_my_postRT2_bis(Jout,Kt,Ae1x,Ae1y,Ae1z,d);
%% 
rho=fun_J_to_coeff([normJ2.^0.5;0]);
%% 
rho_eF=0.5*(abs(Ae(:,:)).'*rho(1:end-1));
z_realF=rho_eF;
indFneq=setdiff([1:3*Kt].',idxF);
z_realF(indFneq,:)=0;
z_realx=zeros(L,M,N);
z_realx(idxFx)=z_realF(idxFx);
z_realy=zeros(L,M,N);
z_realy(idxFy)=z_realF(Kt+idxFy);
z_realz=zeros(L,M,N);
z_realz(idxFz)=z_realF(2*Kt+idxFz);
%% 
[H] = fun_constructH_for_Jacobian_bis(x,d,idxV,idxF,idxFx,idxFy,...
         idxFz,K,Kt,Ae,L,M,N,fun_J2_to_dcoeff,num_curr,normJ2);
%% 
% JAC=[R+H       AeeR.';...
%      P*AeeR  -eye(size(AeeR,1))];
res =  matvect_mult6mag_jacT_T(F,opCirculantP_all,z_realx,z_realy,z_realz,...
                        idxF,idxV,d,Aee,L,M,N,AeeR,idxVR,...
                        H);
end

