function [vsol] = fun_resolve_Jac(x,d,idxV,idxF,idxFx,idxFy,...
                            idxFz,Kt,Ae,Ae1x,Ae1y,Ae1z,L,M,N,...
                            fun_J_to_coeff,...
                            opCirculantP_all,Aee,AeeR,idxVR,FF,...
                            fun_J2_to_dcoeff,K,st_sparse_preconP,fl_precon_type,solver_flag)
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
fMVM = @(J) matvect_mult6mag_jac(J,opCirculantP_all,z_realx,z_realy,z_realz,...
                        idxF,idxV,d,Aee,L,M,N,AeeR,idxVR,...
                        H);
%% 
prec_mod = 'diag'; % 'diag' 'allH'; only 'diag' works now
[Y_inv,P_diag,D_diag,LL,UU,PP,QQ,RR] = preparePREC_MAG_JAC_NEW(d,z_realF,idxFx,idxFy,idxFz,st_sparse_preconP,AeeR,Aee,Kt,...
                          H,prec_mod,fl_precon_type);
%fPMV = @(JOut_full_in)sparse_precon_multiply2(JOut_full_in, Aee, [], [], prectol,AeeR);
fPMV = @(JOut_full_in)multiplyPREC_CAP_NEW(JOut_full_in,AeeR,fl_precon_type,Y_inv,P_diag,LL,UU,PP,QQ,RR);
%% solving 
global NLcounter
if NLcounter==1
tol = 1e-4;
inner_it = 100;
outer_it = 1;
else
tol = 1e-4;
inner_it = 100;
outer_it = 1;    
end


if strcmp(solver_flag,'gmres')
    [vsol, flag, relres, iter, resvec] =  pgmres(@(J)fMVM(J),FF, inner_it, tol, outer_it, @(JOut_full_in)fPMV(JOut_full_in) );
elseif strcmp(solver_flag,'tfqmr') 
    [vsol, flag, relres, iter, resvec] =  RTtfqmr(@(J)fMVM(J),FF, tol, inner_it, @(JOut_full_in)fPMV(JOut_full_in) );        
elseif strcmp(solver_flag,'bicgstab') 
    [vsol, flag, relres, iter, resvec] =  RTbicgstab(@(J)fMVM(J),FF, tol, inner_it, @(JOut_full_in)fPMV(JOut_full_in)); % per chimarlo con un x0 ,@(JOut_full_in) JOut_full_in,vsol);        
end




NLcounter = NLcounter+1;
end

