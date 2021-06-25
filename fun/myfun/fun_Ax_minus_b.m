function [FUN] = fun_Ax_minus_b(x,d,idxV,idxF,idxFx,idxFy,...
                            idxFz,Kt,Ae,Ae1x,Ae1y,Ae1z,L,M,N,...
                            fun_J_to_coeff,...
                            opCirculantP_all,Aee,AeeR,idxVR,rhs_vectR)
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
FUN=multiplyMATVECT_MAG(x,opCirculantP_all,z_realx,z_realy,z_realz,idxF,idxV,d,Aee,L,M,N,AeeR,idxVR);
FUN=FUN-rhs_vectR;
end

