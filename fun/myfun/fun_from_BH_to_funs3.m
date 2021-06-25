function [fun_coeffmmu0M,fun_dcoeffm_dmu0M2] = fun_from_BH_to_funs3(BH_tab)
mytype='linearinterp';%'cubicinterp'; % 'linearinterp'
         
%         H             B 
% BH_tab=[0	            0
        % 663.146	    1
        % 1067.5	    1.1
        % 1705.23	    1.2
        % 2463.11	    1.3
        % 3841.67	    1.4
        % 5425.74	    1.5
        % 7957.75	    1.6
        % 12298.3	    1.7
        % 20462.8	    1.8
        % 32169.6	    1.9
        % 61213.4	    2
        % 111408	    2.1
        % 188487.757	2.2
        % 267930.364	2.3
        % 347507.836	2.4];

mu_0=4*pi*1e-7;
%%
figure
plot(BH_tab(:,1),BH_tab(:,2),'x-')
hold on 
title('BH')
xlabel('H')
ylabel('B')
%% from fun_BH to fun_murH -->mur=fun_murH(H)
murH_tab=[BH_tab(:,1),BH_tab(:,2)./BH_tab(:,1)/mu_0];
murH_tab(1,2)=murH_tab(2,2);
% figure
% plot(murH_tab(:,1),murH_tab(:,2),'x-')
% hold on 
% title('fun \mu_r H')
% ylabel('\mu_r')
% xlabel('H')
%% from fun_murH to fun_chimH -->mur-1=fun_chimH(H)
chimH_table=murH_tab; 
chimH_table(:,2)=chimH_table(:,2)-1;
% figure
% plot(chimH_table(:,1),chimH_table(:,2),'x-')
% hold on 
% title('fun \chi_m H')
% ylabel('\chi_m')
% xlabel('H')
%% from fun_chimH fun_mu0MH  mu0M=mu0*(mur-1)*H --> mu0M=fun_mu0MH(H)
mu0MH_table = chop([chimH_table(:,1),chimH_table(:,1).*chimH_table(:,2)*mu_0],6);
% okk=ones(size(mu0MH_table,1),1);
% for ii = 2:size(mu0MH_table,1)
%     if (mu0MH_table(ii,2)-mu0MH_table(ii-1,2))<0
%        okk(ii)=0; 
%     end
% end
% figure
% plot(mu0MH_table(:,1),mu0MH_table(:,2),'x-')
% hold on 
% title('fun \mu_0 M H')
% ylabel('\mu_0 M')
% xlabel('H')
%% from fun_mu0MH to fun_Hmu0M --> H= fun_Hmu0M(mu0M)
% mu0Mrange=mu_0*Hrange.*fun_chimM(Hrange);
% mu0Mrange=fun_mu0MH(Hrange);
% [mu0Mrange,id]=unique(mu0Mrange);
% Hrange=Hrange(id);
% fun_Hmu0M = fit(mu0Mrange,Hrange,mytype);
Hmu0M_table=mu0MH_table(:,[2,1]);
% figure
% plot(mu0Mrange,fun_Hmu0M(mu0Mrange),'.')
% plot(Hmu0M_table(:,1),Hmu0M_table(:,2),'x-')
% hold on 
% title('fun H \mu_0 M')
% ylabel('H')
% xlabel('\mu_0 M')
%% from fun_Hmu0M to fun_coeffm_mu0M   H=coeffm*mu0Mrange --> 
% coeffm = fun_coeffmmu0M(mu0M)     coeffm=1/(mu0*(mur-1)) 
% cresce_mu0Mrange diminuisce_(mur-1) cresce_coeffm
% Hrange=linspace(0,max(BH_tab(:,1))*10,10000).';
% mu0Mrange=fun_mu0MH(Hrange);
% [mu0Mrange,id]=unique(mu0Mrange);
% Hrange=Hrange(id);
% fun_coeffmmu0M = fit(mu0Mrange,Hrange./mu0Mrange,mytype);
% Mrange2=linspace(min(mu0Mrange),max(mu0Mrange),50000);
coeffmmu0M_table=[Hmu0M_table(:,1),Hmu0M_table(:,2)./Hmu0M_table(:,1)];
coeffmmu0M_table(1,2)=coeffmmu0M_table(2,2);
% figure
% plot(Mrange2,fun_coeffmmu0M(Mrange2),'bx')
% plot(coeffmmu0M_table(:,1),coeffmmu0M_table(:,2),'x-')
% hold on % plot(mu0Mrange,fun_Hmu0M(mu0Mrange)./mu0Mrange,'ro')
% title('fun coeffmM')
% ylabel('coeffM')
% xlabel('\mu_0 M')
%% from fun_coeffmmu0M to fun_dcoeffm_dmu0M -->
% dcoeffm_dmu0M = fun_dcoeffm_dmu0M(mu0M)
% Mrange2=linspace(min(mu0Mrange),max(mu0Mrange),10000);
% delta=(Mrange2(2)-Mrange2(1))/100;
% fun_dcoeffm_dmu0M=@(yy)( (fun_coeffmmu0M(yy) - fun_coeffmmu0M(yy-delta))./(delta));
dcoeffm_dmu0M=(coeffmmu0M_table(2:end,2)-coeffmmu0M_table(1:end-1,2))./...
              (coeffmmu0M_table(2:end,1)-coeffmmu0M_table(1:end-1,1));
dcoeffm_dmu0M_table=[0.5*(coeffmmu0M_table(2:end,1)+coeffmmu0M_table(1:end-1,1)),...
                     dcoeffm_dmu0M];
% figure
% plot(dcoeffm_dmu0M_table(:,1),dcoeffm_dmu0M_table(:,2),'x-')
% hold on 
% title('fun dcoeffm dmu0M')
% ylabel('fun dcoeffm dmu0M')
% xlabel('\mu_0 M')
%% from  fun_coeffm_mu0M  fun_coeffm_mu0M^2 -->
% coeffm = fun_coeffmmu0M(mu0M.^2)
% fun_coeffmM2 = fit((Mrange2.^2).',fun_coeffmmu0M(Mrange2),mytype);
coeffmmu0M2_table=[coeffmmu0M_table(:,1).^2,coeffmmu0M_table(:,2)];
% figure
% plot(coeffmmu0M2_table(:,1),coeffmmu0M2_table(:,2),'x-')
% hold on
% title('fun coeffmM^2')
% ylabel('coeffM')
% xlabel('||\mu_0 M||^2')
%% from fun_coeffmmu0M2 to fun_dcoeffm_dmu0M2 -->
% dcoeffm_dmu0M = fun_dcoeffm_dmu0M2(mu0M.^2)
% Mrange3=linspace(min(mu0Mrange.^2),max(mu0Mrange.^2),10000);
% delta=(Mrange3(2)-Mrange3(1))/100;
% fun_dcoeffm_dmu0M2=@(yy)( (fun_coeffmM2(yy) - fun_coeffmM2(yy-delta))./delta);
dcoeffm_dmu0M2=(coeffmmu0M2_table(2:end,2)-coeffmmu0M2_table(1:end-1,2))./...
              (coeffmmu0M2_table(2:end,1)-coeffmmu0M2_table(1:end-1,1));
dcoeffm_dmu0M2_table=[0.5*(coeffmmu0M2_table(2:end,1)+coeffmmu0M2_table(1:end-1,1)),...
                     dcoeffm_dmu0M2];

% figure
% plot(dcoeffm_dmu0M2_table(:,1),dcoeffm_dmu0M2_table(:,2),'x-')
% hold on 
% title('fun dcoeffm dmu0M2')
% ylabel('fun dcoeffm dmu0M2')
% xlabel('||\mu_0 M||^2')
%%
[~,id]=unique(coeffmmu0M_table(:,1));
% Hrange=Hrange(id);
fun_coeffmmu0M= fit(coeffmmu0M_table(id,1),coeffmmu0M_table(id,2),mytype);
% figure
% plot(Mrange2,fun_coeffmmu0M(Mrange2),'bx')
% plot(coeffmmu0M_table(:,1),coeffmmu0M_table(:,2),'x-')
% hold on
% plot([coeffmmu0M_table(:,1);max(coeffmmu0M_table(:,1))*1.001],...
%     fun_coeffmmu0M([coeffmmu0M_table(:,1);max(coeffmmu0M_table(:,1))*1.001]),'or')
% title('fun coeffmM')
% ylabel('coeffM')
% xlabel('\mu_0 M')
%
%%
% id=~isinf(dcoeffm_dmu0M2_table(:,2));
% Hrange=Hrange(id);

dcoeffm_dmu0M2=(coeffmmu0M2_table(2:end,2)-coeffmmu0M2_table(1:end-1,2))./...
              (coeffmmu0M2_table(2:end,1)-coeffmmu0M2_table(1:end-1,1));
dcoeffm_dmu0M2_table=[0.5*(coeffmmu0M2_table(2:end,1)+coeffmmu0M2_table(1:end-1,1)),...
                     dcoeffm_dmu0M2];

 deltas=(coeffmmu0M2_table(2:end,1)-coeffmmu0M2_table(1:end-1,1));
 delta=Inf;
 for ii = 1:length(deltas)
     if deltas(ii)<delta && deltas(ii)>0
        delta= deltas(ii);
     end
 end
 delta=delta/100;
 points=zeros(2*size(coeffmmu0M2_table,1)-2,1);
points(1:2:end,1)=coeffmmu0M2_table(1:end-1,1)+delta;                 
points(2:2:end,1)=coeffmmu0M2_table(2:end,1)-delta;  
points(1,1)=coeffmmu0M2_table(1,1);
points(end,1)=coeffmmu0M2_table(end,1);
val=zeros(2*size(coeffmmu0M2_table,1)-2,1);
val(1:2:end,1)=dcoeffm_dmu0M2_table(1:end,2);                 
val(2:2:end,1)=dcoeffm_dmu0M2_table(1:end,2);

id=~isinf(val);

points2=points(id);%([id;id]);
val2=val(id);%([id;id]);

fun_dcoeffm_dmu0M2= fit(points2,val2,mytype);
% figure
% plot(Mrange2,fun_coeffmmu0M(Mrange2),'bx')
% plot(points2,val2,'x-')
% hold on
% plot(points2,...
%     fun_dcoeffm_dmu0M2(points2),'or')
% title('fun dcoeffm dmu0M2')
% ylabel('fun dcoeffm dmu0M2')
% xlabel('||\mu_0 M||^2')
%

end

