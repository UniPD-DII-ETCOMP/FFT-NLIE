function [H] = fun_constructH_for_Jacobian_bis(x,d,idxV,idxF,idxFx,idxFy,idxFz,K,Kt,Ae,L,M,N,fun_mat,...
                                      num_curr,normJ2)
%  
% H = Hxx Hxy Hxz
%     Hxy Hyy Hyz
%     Hxz Hyz Hzz
%% 
drho_dnormJ2=fun_mat([normJ2;0]);
% 
drho_dnormJ2_medio=0.5*(abs(Ae(:,idxF)).'*drho_dnormJ2(1:end-1));
%%
dx = d(1); 
dy = d(2); 
dz = d(3);
Ax=(dy*dz);
Ay=(dx*dz);
Az=(dy*dx);
Vol=dx*dy*dz;

Nfx=length(idxFx);
Nfy=length(idxFy);
Nfz=length(idxFz);

coeffxx=(1/Ax^4)*2*Vol;
coeffyy=(1/Ay^4)*2*Vol;
coeffzz=(1/Az^4)*2*Vol;

%%
Hxx = sparse(1:Nfx,1:Nfx,coeffxx.*(x(1:Nfx).^2).*drho_dnormJ2_medio(1:Nfx),Nfx,Nfx);
Hyy = sparse(1:Nfy,1:Nfy,coeffyy.*(x(Nfx+1:Nfx+Nfy).^2).*drho_dnormJ2_medio(Nfx+1:Nfx+Nfy),Nfy,Nfy);
Hzz = sparse(1:Nfz,1:Nfz,coeffzz.*(x(Nfx+Nfy+1:Nfx+Nfy+Nfz).^2).*drho_dnormJ2_medio(Nfx+Nfy+1:Nfx+Nfy+Nfz),Nfz,Nfz);
%% 
bHxy=(Kt+1)*ones(4,Kt);
vHxy=(Kt+1)*ones(4,Kt);
bHxz=(Kt+1)*ones(4,Kt);
vHxz=(Kt+1)*ones(4,Kt);
bHyz=(Kt+1)*ones(4,Kt);
vHyz=(Kt+1)*ones(4,Kt);
for iii = 1:K
    ii = idxV(iii);
    
    lmn = idx2triplet(ii,L,M,N); %map i -> [l m n] triplet
    l=lmn(1);
    m=lmn(2);
    n=lmn(3);
    % i = l+(m-1)*L+(n-1)*L*M
    
    % Hxy
    bHxy(1,ii) =  (l+(m-1)*L+(n-1)*L*M);
    vHxy(1,ii) = (l+(m-1)*L+(n-1)*L*M);
    if l+1 <= L
        j = (l+1)+(m-1)*L+(n-1)*M*L;
        bHxy(2,ii) =  j;
        vHxy(2,ii) = (l+1+(m-1)*L+(n-1)*L*M);
    end
    if m-1 > 0
        j = (l)+(m-1-1)*L+(n-1)*M*L;
        bHxy(3,ii) =  j;
        vHxy(3,ii) = (l+(m-1)*L+(n-1)*L*M);
    end
    if (l+1 <= L) && (m-1 > 0)
        j = (l+1)+(m-1-1)*L+(n-1)*M*L;
        bHxy(4,ii) =  j;
        vHxy(4,ii) = (l+1+(m-1)*L+(n-1)*L*M);
    end    
    
    % Hxz
    bHxz(1,ii) =  (l+(m-1)*L+(n-1)*L*M);
    vHxz(1,ii) = (l+(m-1)*L+(n-1)*L*M);
    if l+1 <= L
        j = (l+1)+(m-1)*L+(n-1)*M*L;
        bHxz(2,ii) = j;
        vHxz(2,ii) = (l+1+(m-1)*L+(n-1)*L*M);
    end
    if n-1 > 0
        j = (l)+(m-1)*L+(n-1-1)*M*L;
        bHxz(3,ii) =  j;
        vHxz(3,ii) = (l+(m-1)*L+(n-1)*L*M);
    end
    if (l+1 <= L) && (n-1 > 0)
        j = (l+1)+(m-1)*L+(n-1-1)*M*L;
        bHxz(4,ii) =  j;
        vHxz(4,ii) = (l+1+(m-1)*L+(n-1)*L*M);
    end       
    
    % Hyz
    bHyz(1,ii) = (l+(m-1)*L+(n-1)*L*M);
    vHyz(1,ii) = (l+(m-1)*L+(n-1)*L*M);
    if m+1 <= M
        j = (l)+(m+1-1)*L+(n-1)*M*L;
        bHyz(2,ii) = j;
        vHyz(2,ii) = (l+(m+1-1)*L+(n-1)*L*M);
    end
    if n-1 > 0
        j = (l)+(m-1)*L+(n-1-1)*M*L;
        bHyz(3,ii) =  j;
        vHyz(3,ii) = (l+(m-1)*L+(n-1)*L*M);
    end
    if (m+1 <= M) && (n-1 > 0)
        j = (l)+(m+1-1)*L+(n-1-1)*M*L;
        bHyz(4,ii) =  j;
        vHyz(4,ii) = (l+(m+1-1)*L+(n-1)*L*M);
    end 
end
%% 

x_x=x(1:Nfx,1);
x_y=x(Nfx+1:Nfx+Nfy,1);
% x_z=x(Nfx+Nfy+1:Nfx+Nfy+Nfz,1);
x_all=zeros(3*Kt+1,1);
x_all(idxF,1)=x(1:num_curr);

Hxy=sparse(Nfx,Kt+1);
for ii = 1:4
coeff=(1/(Ax^2))*(1/(Ay^2))*2*(x_x.*x_all(Kt+bHxy(ii,idxFx)))*Vol.*0.25.*drho_dnormJ2(vHxy(ii,idxFx));
Hxy=Hxy+sparse(1:Nfx,bHxy(ii,idxFx),coeff,Nfx,Kt+1);    
end
Hxy=Hxy(:,idxFy);

Hxz=sparse(Nfx,Kt+1);
for ii = 1:4
coeff=(1/(Ax^2))*(1/(Az^2))*2*(x_x.*x_all(2*Kt+bHxz(ii,idxFx)))*Vol.*0.25.*drho_dnormJ2(vHxz(ii,idxFx));
Hxz=Hxz+sparse(1:Nfx,bHxz(ii,idxFx),coeff,Nfx,Kt+1);    
end
Hxz=Hxz(:,idxFz);


Hyz=sparse(Nfy,Kt+1);
for ii = 1:4
coeff=(1/(Ay^2))*(1/(Az^2))*2*(x_y.*x_all(2*Kt+bHyz(ii,idxFy)))*Vol.*0.25.*drho_dnormJ2(vHyz(ii,idxFy));
Hyz=Hyz+sparse(1:Nfy,bHyz(ii,idxFy),coeff,Nfy,Kt+1);    
end
Hyz=Hyz(:,idxFz);

H=[Hxx   Hxy   Hxz;...
   Hxy.' Hyy   Hyz;...
   Hxz.' Hyz.' Hzz];

end

