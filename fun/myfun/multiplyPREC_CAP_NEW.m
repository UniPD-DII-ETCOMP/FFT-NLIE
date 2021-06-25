%% Sparse Preconditioner Multiplication: Used For Capacitive Effects
%*************************************************************************

function [JOut_full_out] = multiplyPREC_CAP_NEW(JOut_full_in,AeR,fl_precon_type,Y_inv,P_diag,LL,UU,PP,QQ,RR)

fl_profile = 0;

tic

if (strcmp(fl_precon_type, 'no_prec') == 1)
    JOut_full_out = JOut_full_in;
    return
end

num_nodeR=size(AeR,1); % 
num_curr=size(AeR,2);
JOut_full_out=zeros(num_nodeR+num_curr,1);
   
%     if (fl_volt_source == 1) % voltage source
%         disp('Change preconditioner type !!! ')
%         error('These decompositions can not work with fl_volt_source=1 !!!')
%     end
           
%%Calculate output vectors:
%d = inv(S)*[b - P*A*inv(Z)*a]
%c = inv(Z)*[a - A'*d]
               
%%Calculate d vector
warning off
JOut_full_out(num_curr+1:num_curr+num_nodeR) = ...
    QQ * (UU \ (LL \ (PP * (RR \ (JOut_full_in(num_curr+1:num_curr+num_nodeR) - ( P_diag * (AeR * (Y_inv * JOut_full_in(1:num_curr))) ) )))));
warning on
%%Calculate c vector
JOut_full_out(1:num_curr) = Y_inv*(JOut_full_in(1:num_curr) - ((AeR.')*JOut_full_out(num_curr+1:num_curr+num_nodeR)));
                       
if(fl_profile == 1); disp(['Time for matvect - sparse preconditioner part::: ',num2str(toc)]); end

end