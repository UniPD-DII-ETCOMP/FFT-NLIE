%% FFT-Non_Linear_Inductance_Extractor: 
close all
clear global
clear
clc
restoredefaultpath
warning on
format short
warning off
delete log_Green_K.txt
warning on
%%
%% BEGIN USER SETTINGS
%%
%% Directory
name_dir='test';
%% Selections
plot_vectorsJ_flag = 1; %quiver plot of real and imag of J
plot_potential_flag = 1; %color plot of phi real and imag
paraview_export_flag = 1; % export to paraviw
refine.flag = 0; refine.x=1; refine.y=1; refine.z=1; % refine
Integration_flag = 'NumAn'; %'NumAn'; 'NumNum' (Integration: NumericalNumerical or AnalyticalNumerical)
%% Solver parameters
tol = 1e-6;
inner_it = 40;
outer_it = 5;
NR_tol = 1e-3;
NR_iter=50;
%%
%% END USER SETTINGS
%%
%% Add Path
dad = pwd;
cd('fun'); addpath(genpath(pwd)); cd(dad)
cd('fortran'); addpath(pwd); cd(dad)
cd('data'); cd(name_dir); load('data.mat'); 
fileList = dir('*.stl');
figure
hold on
xmin=[];xmax=[];ymin=[];ymax=[];zmin=[];zmax=[];ccolor=distinguishable_colors(size(fileList,1));
for ii = 1:size(fileList,1)
    [stlcoords] = READ_stl(fileList(ii).name);
    xco = squeeze( stlcoords(:,1,:) )';
    yco = squeeze( stlcoords(:,2,:) )';
    zco = squeeze( stlcoords(:,3,:) )';
    [hpat] = patch(xco,yco,zco,ccolor(ii,:));
    axis equal
    xlabel('x');
    ylabel('y');
    zlabel('z');
    view(3)
    title('stl (original, not scaled)')
    drawnow
end
cd(dad)
modelname = name_dir;
%% refine
if refine.flag
    warning('refine on')
    mymod=1;
    for ii = 1:refine.x
        [Ind,L,M,N,xyz,smeshx,smeshy,smeshz,Nmat,nVoxel] = fun_refine(Ind,xyz,smeshx,smeshy,smeshz,Nmat,L,M,N,1,mymod);
    end
    for ii = 1:refine.y
        [Ind,L,M,N,xyz,smeshx,smeshy,smeshz,Nmat,nVoxel] = fun_refine(Ind,xyz,smeshx,smeshy,smeshz,Nmat,L,M,N,2,mymod);
    end   
    for ii = 1:refine.z
        [Ind,L,M,N,xyz,smeshx,smeshy,smeshz,Nmat,nVoxel] = fun_refine(Ind,xyz,smeshx,smeshy,smeshz,Nmat,L,M,N,3,mymod);
    end      
end
%% EM constants
mu = 4*pi*1e-7;
co = 299792458;
eo = 1/co^2/mu;
%% extract data information 
rhoVoxel=zeros(nVoxel,1);
idxVe=[]; rhomin=Inf; ind_c=[]; val_c=[]; k=1;
idxVm=[];
for ii = 1:Nmat
    Ind(ii).ind=reshape(Ind(ii).ind,length(Ind(ii).ind),1);
    if strcmp(Ind(ii).tag,'air') ||  strcmp(Ind(ii).tag,'diel')
        % nothing to do here (?)
    elseif strcmp(Ind(ii).tag,'cond')
        idxVe=[idxVe;Ind(ii).ind];  
        rhoVoxel(Ind(ii).ind,1)=Ind(ii).rho;  
        rhomin=min([rhomin,Ind(ii).rho]);
    elseif strcmp(Ind(ii).tag,'mag')
        idxVm=[idxVm;Ind(ii).ind];  
        BH_tab=Ind(ii).BH_tab;
    elseif strcmp(Ind(ii).tag,'port')
        idxVe=[idxVe;Ind(ii).ind]; 
        rhoVoxel(Ind(ii).ind,1)=Ind(ii).rho;
        rhomin=min([rhomin,Ind(ii).rho]);
        ind_c(k)=Ind(ii).ind(1);
        val_c(k)=Ind(ii).cur;
        k=k+1;
    end
end
idxVe=unique(idxVe);
idxVm=unique(idxVm);
%% Grid Definition
disp('----DOMAIN--------------------------------')
%%Grid resolution
disp([' Number of voxels in x direction: ', num2str(L)])
disp([' Number of voxels in y direction: ', num2str(M)])
disp([' Number of voxels in z direction: ', num2str(N)])
disp(' Resolution:')
dx = smeshx; dy = smeshy; dz = smeshz;
disp([' dx = ',num2str(dx),' m']); disp([' dy = ',num2str(dy),' m']); disp([' dz = ',num2str(dz),' m'])
d = [dx dy dz]; 
Kt = nVoxel; %total number of voxels
K = length(idxVe); %number of non-empty voxels
Km = length(idxVm); %number of non-empty voxels
%% Set Material Properties
rho_eV=reshape(rhoVoxel,L,M,N); % 
clear rhoVoxel
%%
disp([' Total number of voxels: ', num2str(Kt)])
disp([' Number of non-empty voxels: ', num2str(K)])
disp(' ')
%% Incidence Matix A
disp('----COMPUTING INCIDENCE--------------------------------')
mytic=tic;
[Ae,Aee,idxF,idxFx,idxFy,idxFz,Ae1x,Ae1y,Ae1z] = ...
    incidence_matrix3(Kt,[L M N],idxVe);
[~,Amm,idxFm,idxFxm,idxFym,idxFzm,~,~,~] = ...
    incidence_matrix3(Kt,[L M N],idxVm);
disp([' Number of DoFs (cond): ', num2str(size(Aee,1)+size(Aee,2))])
disp([' Number of DoFs (mag): ', num2str(size(Amm,1)+size(Amm,2))])
disp([' Time for computing incidence ::: ' ,num2str(toc(mytic))]);
disp(' ')
%% Forcing Term
Vx=zeros(L*M*N,1);
Vy=zeros(L*M*N,1);
Vz=zeros(L*M*N,1);
%% Matrices Z_real and Z_imag
rho_eF=0.5*(abs(Ae(:,:)).'*rho_eV(:)); clear rho_eV
z_realF=rho_eF;
indFneq=setdiff([1:3*Kt].',idxF);
z_realF(indFneq,:)=0;
z_realx=zeros(L,M,N);
z_realx(idxFx)=z_realF(idxFx);
z_realy=zeros(L,M,N);
z_realy(idxFy)=z_realF(Kt+idxFy);
z_realz=zeros(L,M,N);
z_realz(idxFz)=z_realF(2*Kt+idxFz);
%
%% Compute Green Tensor 
disp('----COMPUTING GREEN TENSOR--------------------------------')
mytic_G=tic;
[Gmn] = computeGREEN(d,L,M,N,Integration_flag);
[Kmn] = computeGREEN_GK(d,L,M,N,Integration_flag,'mex');
[opCirculantL_all,st_sparse_preconL] = computeCIRCULANT(Gmn,d,'L');
opCirculantL_all = (mu)*opCirculantL_all;
disp([' Time for getting Green tensor ::: ' ,num2str(toc(mytic_G))]);
disp(' ')
%% Compute Circulant Tensors
disp('----COMPUTING CIRCULANT TENSOR--------------------------------')
disp(' Circulant Tensors related to P,L matrices')
mytic_cir=tic;
[opCirculantP_all,st_sparse_preconP] = computeCIRCULANT(Gmn,d,'P');
opCirculantP_all=opCirculantP_all/mu;
st_sparse_preconP=st_sparse_preconP/mu;
[opCirculantK_all] = computeCIRCULANT_K(Kmn);
%%Add constants to Circulants
disp([' Time for getting circulant tensors ::: ' ,num2str(toc(mytic_cir))])
clear Gmn %Green tensor is not used anymore
disp(' ')
%% Generating RHS vector
num_node = size(Aee,1); %all potential nodes in non-empty voxels 
num_curr = size(Aee,2); %all currens in non-empty voxels 
num_nodem = size(Amm,1); %all potential nodes in non-empty voxels 
num_currm = size(Amm,2); %all currens in non-empty voxels 
%%Define RHS: (injected currents)
iinj=zeros(L*M*N,1);
iinj(ind_c)=val_c;
iinj=iinj(idxVe);
rhs_vect = [Vx(idxFx);Vy(idxFy);Vz(idxFz);-iinj]; 
clear Vx Vy Vz
%% Computing Preconditioner
disp('----COMPUTING PRECONDITIONER--------------------------------')
mytic_prec=tic;
[Y_inv,LL,UU,PP,QQ,RR] = preparePREC_NEW(d,z_realF,idxFx,...
    idxFy,idxFz,Aee,Kt);
fPMV = @(JOut_full_in)multiplyPREC_NEW(JOut_full_in,Aee,Y_inv,LL,UU,PP,QQ,RR);
disp([' Time for computing preconditioner ::: ' ,num2str(toc(mytic_prec))]);
disp(' ')
%% Solution of Linear System
disp('----SOLVING LINEAR SYSTEM-------------------------------')
zloc=[(dx/(dy*dz))*z_realx(idxFx);...
      (dy/(dx*dz))*z_realy(idxFy);...
      (dz/(dx*dy))*z_realz(idxFz)];
fMVM = @(J) multiplyMATVECT(J,zloc,Aee);
mytic_solver=tic;
[vsol] = pgmres_mod(@(J)fMVM(J),rhs_vect, inner_it, tol, outer_it, @(JOut_full_in)fPMV(JOut_full_in) );  
disp([' Time for solving system with gmres ::: ' ,num2str(toc(mytic_solver))]);
disp(' ')
%% extract solution
Jout = zeros(L,M,N,3);
Jout(idxF) = vsol(1:num_curr) ; % return to global variables
%%
[A,a_dof] = post_processing_A(Jout,xyz,d,opCirculantL_all,Ae1x,Ae1y,Ae1z);
Wm=0.5*0.5*sum((vsol(1:num_curr).*conj(a_dof(idxF))));
L1=2*2*Wm/(abs(val_c(1))^2);
disp('-------------------------------------------------------------------')
disp(' ')
disp(['  Coil Inductance ' ,num2str(L1),' H [without magnetic core]']);
disp(' ')
disp('-------------------------------------------------------------------')
%%
[H] = post_processing_H(Jout,xyz,d,Ae1x,Ae1y,Ae1z,opCirculantK_all);   
Hx = reshape(H(:,1),L,M,N);
Hy = reshape(H(:,2),L,M,N);
Hz = reshape(H(:,3),L,M,N);
Gram = dx*dy*dz; %volume of cubic element
Vmx = (Gram.*Hx)./(dy*dz);
Vmy = (Gram.*Hy)./(dx*dz);
Vmz = (Gram.*Hz)./(dx*dy);
clear Hx Hy Hz 
rhs_vectm = [Vmx(idxFxm);Vmy(idxFym);Vmz(idxFzm);zeros(num_nodem,1)]; 
%%
[fun_J_to_coeff,fun_J2_to_dcoeff] = fun_from_BH_to_funs3(BH_tab);
fun2.F=@(x) fun_Ax_minus_b ... 
            (x,d,idxVm,idxFm,idxFxm,idxFym,...
                            idxFzm,Kt,Ae,Ae1x,Ae1y,Ae1z,L,M,N,...
                            fun_J_to_coeff,...
                            opCirculantP_all,Amm,Amm,idxVm,rhs_vectm);
fun2.solvJ=@(x,FF) fun_resolve_Jac(x,d,idxVm,idxFm,idxFxm,idxFym,...
                            idxFzm,Kt,Ae,Ae1x,Ae1y,Ae1z,L,M,N,...
                            fun_J_to_coeff,...
                            opCirculantP_all,Amm,Amm,idxVm,FF,...
                            fun_J2_to_dcoeff,Km,st_sparse_preconP,'schur_invert','gmres'); % funzione che risolve il problema con jacobiano Jac(x)\b
fun2.JTFF_T = @(x,FF) fun_JacTdotF_T(x,d,idxVm,idxFm,idxFxm,idxFym,...
                            idxFzm,Kt,Ae,Ae1x,Ae1y,Ae1z,L,M,N,...
                            fun_J_to_coeff,...
                            opCirculantP_all,Amm,Amm,idxVm,FF,...
                            fun_J2_to_dcoeff,Km); 
global NLcounter
NLcounter=1;
options = optimset('TolFun',NR_tol,'MaxIter',NR_iter); % set TolX
mytic_sol_NR=tic;
[vsolm, resnorm, f, exitflag, ...
    output, jacob,resnormstore,xstore] ...
    = newtonraphson_mod2(fun2, zeros(num_currm+num_nodem,1), options);
disp([' Total time for solving NON LINEAR SYSTEM ::: ' ,num2str(toc(mytic_sol_NR))])                        
%% extract solution
Joutm = zeros(L,M,N,3);
Joutm(idxFm) = vsolm(1:num_currm) ; % return to global variables
%%
[~,a_dof2] = post_processing_H(Joutm,xyz,d,Ae1x,Ae1y,Ae1z,opCirculantK_all);
Wm=0.5*0.5*sum(((vsol(1:num_curr)).*conj(a_dof(idxF)+a_dof2(idxF))));
L2=2*2*Wm/(abs(val_c(1))^2);
disp('-------------------------------------------------------------------')
disp(' ')
disp(['  Coil Inductance ' ,num2str(L2),' H [with magnetic core]']);
disp(' ')
disp('-------------------------------------------------------------------')
%%
%% POST PROCESSING
%%
%% Post Processing J
disp('----POST PROCESSING J------------------------------')
mytic_prec=tic;
[J,XYZ] = fun_my_postRT2(Jout,Kt,Ae1x,Ae1y,Ae1z,xyz,L,M,N,d);
[MM,XYZ] = fun_my_postRT2(Joutm,Kt,Ae1x,Ae1y,Ae1z,xyz,L,M,N,d);
potval=zeros(Kt,1);
potval(idxVe)=vsol(num_curr+1:end);
disp([' Total time for post processing J and mu0M ::: ' ,num2str(toc(mytic_prec))]);
disp(' ')
%% Plot Vectors  J
if plot_vectorsJ_flag
    jjR = (J);
    figure
    normJR=sqrt(jjR(:,1).^2+jjR(:,2).^2+jjR(:,3).^2);
    quiver3_c_scal(XYZ(:,1),XYZ(:,2),XYZ(:,3),jjR(:,1),jjR(:,2),jjR(:,3),...
              normJR,4);
    axis equal
    c1=colorbar;
    caxis([min(normJR) max(normJR)]);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('Current Density Vector ')
    c1.Location = 'southoutside';
    xlim([min(XYZ(:,1))-dx max(XYZ(:,1))+dx])
    ylim([min(XYZ(:,2))-dy max(XYZ(:,2))+dy])
    zlim([min(XYZ(:,3))-dz max(XYZ(:,3))+dz])
    %% Plot Vectors H
    hhR = H;%(M);
    figure
    normHR=sqrt(hhR(:,1).^2+hhR(:,2).^2+hhR(:,3).^2);
    quiver3_c_scal(XYZ(:,1),XYZ(:,2),XYZ(:,3),hhR(:,1),hhR(:,2),hhR(:,3),...
              normHR,4);
    axis equal
    c1=colorbar;
    caxis([min(normHR) max(normHR)]);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('H produced by coil only')
    c1.Location = 'southoutside';
    xlim([min(XYZ(:,1))-dx max(XYZ(:,1))+dx])
    ylim([min(XYZ(:,2))-dy max(XYZ(:,2))+dy])
    zlim([min(XYZ(:,3))-dz max(XYZ(:,3))+dz])
    %
    %% Plot Vectors H
    mmR = MM;
    figure
    normMR=sqrt(mmR(:,1).^2+mmR(:,2).^2+mmR(:,3).^2);
    quiver3_c_scal(XYZ(:,1),XYZ(:,2),XYZ(:,3),mmR(:,1),mmR(:,2),mmR(:,3),...
              normMR,4);
    axis equal
    c1=colorbar;
    caxis([min(normMR) max(normMR)]);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('\mu_0M Vector')
    c1.Location = 'southoutside';
    xlim([min(XYZ(:,1))-dx max(XYZ(:,1))+dx])
    ylim([min(XYZ(:,2))-dy max(XYZ(:,2))+dy])
    zlim([min(XYZ(:,3))-dz max(XYZ(:,3))+dz])
    %
end
warning off
delete log_Green_K.txt
warning on
