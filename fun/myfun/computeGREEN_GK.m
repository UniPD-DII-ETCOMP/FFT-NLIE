%% Compute K Tensor
%**************************************************************************
%INPUT: - d = [dx dy dz] array of grid steps
%       - L,M,N: grid dimensions (num of voxels)
%       - Integration: method can be analytical/numerical or numerical/numerical
%**************************************************************************
function [K] = computeGREEN_GK(d,L,M,N,Integration,mex_flag)
%%Integration Properties
NGx = 3; %number of integration points x-dir
NGy = 3; %number of integration points y-dir
NGz = 3; %number of integration points z-dir
K = zeros(L,M,N,3); 
dx = d(1); dy = d(2); dz = d(3);
infofN = whos('K'); memestimated = 2*infofN.bytes/(1024*1024);
disp([' Memory for temporarily storing K (MB) ::: ' , num2str(memestimated)]);
disp(' Start computing K tensor')
mytic_k = tic;
if strcmp(Integration,'NumAn') 
    [face,face_diff,r_f,r_e,u_e,triang_face,idxf_dot,idxf_cross] = prepareINTEGRATION(d);
    [x,wx]=lgwt(NGx,-dx/2,dx/2);
    [y,wy]=lgwt(NGy,-dy/2,dy/2);
    [z,wz]=lgwt(NGz,-dz/2,dz/2);
    try 
        if strcmp(mex_flag,'mex')
            [K1,K2,K3]=computeGREEN_K_f90_mexed(L,M,N,x,NGx,face,r_f,r_e,u_e,triang_face,idxf_dot,...
                                    idxf_cross,wx,d.',22,NGy,NGz,wy,wz,y,z); 
            K(:,:,:,1)=reshape(K1,L,M,N);
            K(:,:,:,2)=reshape(K2,L,M,N);
            K(:,:,:,3)=reshape(K3,L,M,N);
        else
            error('mex disabled')
        end
    catch
        warning('mex function not supported, try to re-mex it')
        for mx = 1:L
                for my = 1:M
                    for mz = 1:N
                        m = [mx,my,mz];
                        %Centre of the observation voxel
                        r_m = (m-1) .* d;
                        % Kxy
                        xyz_trg=zeros(NGy,3);
                        xyz_trg(:,1) = r_m(1); %Gauss points 
                        xyz_trg(:,2) = y + r_m(2); %Gauss points 
                        xyz_trg(:,3) = r_m(3)+dz/2; %Gauss points 
                        res1 =  Integrate_NumAn3(face,face_diff,r_f,r_e,u_e,triang_face,idxf_dot,idxf_cross,xyz_trg,wy);
                        xyz_trg(:,3) = r_m(3)-dz/2; %Gauss points 
                        res2 =  Integrate_NumAn3(face,face_diff,r_f,r_e,u_e,triang_face,idxf_dot,idxf_cross,xyz_trg,wy);
%                         Kxy(mx,my,mz,:) = -res1+res2;
                        Kxy = (-res1+res2)...
                                          *(1/(dx*dz))...% shape function in the source volume   
                                          *(dx/(dy*dz)); % projection;
                        K(mx,my,mz,3) = Kxy;                        
                        % Kxz
                        xyz_trg=zeros(NGz,3);
                        xyz_trg(:,1) = r_m(1); %Gauss points 
                        xyz_trg(:,2) = r_m(2)+dy/2; %Gauss points 
                        xyz_trg(:,3) = z + r_m(3); %Gauss points 
                        res1 =  Integrate_NumAn3(face,face_diff,r_f,r_e,u_e,triang_face,idxf_dot,idxf_cross,xyz_trg,wz);
                        xyz_trg(:,2) = r_m(2)-dy/2; %Gauss points  
                        res2 =  Integrate_NumAn3(face,face_diff,r_f,r_e,u_e,triang_face,idxf_dot,idxf_cross,xyz_trg,wz);
%                         Kxz(mx,my,mz,:) = res1-res2;  
                        Kxz = (res1-res2)...
                                          *(1/(dx*dy))...% shape function in the source volume   
                                          *(dx/(dy*dz)); % projection;
                        K(mx,my,mz,2) = -Kxz;                         
                        % Kyz
                        xyz_trg=zeros(NGz,3);
                        xyz_trg(:,1) = r_m(1)+dx/2; %Gauss points 
                        xyz_trg(:,2) = r_m(2); %Gauss points 
                        xyz_trg(:,3) = z + r_m(3); %Gauss points 
                        res1 =  Integrate_NumAn3(face,face_diff,r_f,r_e,u_e,triang_face,idxf_dot,idxf_cross,xyz_trg,wz);
                        xyz_trg(:,1) = r_m(1)-dx/2; %Gauss points  
                        res2 =  Integrate_NumAn3(face,face_diff,r_f,r_e,u_e,triang_face,idxf_dot,idxf_cross,xyz_trg,wz);
                        Kyz = (-res1+res2)...
                                          *(1/(dx*dy))...% shape function in the source volume   
                                          *(dy/(dx*dz)); % projection;
                        K(mx,my,mz,1) = Kyz;
                    end
                end
        end
    end
    %Remove unnecessary memory consuming elements
    clear face face_diff r_f r_e u_e triang_face idxf_dot idxf_cross
    clear xyz_src xyz_trg w
elseif strcmp(Integration,'NumNum')
vol = prod(d);
        for mx = 1:L
            for my = 1:M
                for mz = 1:N
                    m = [mx,my,mz];
                    %Centre of the observation voxel
                    r_m = ((m-1) .* d)';
                    % Kxy
                    xyz_trg=zeros(1,3);
                    xyz_trg(:,1) = r_m(1); %Gauss points 
                    xyz_trg(:,2) = r_m(2); %Gauss points 
                    xyz_trg(:,3) = r_m(3)+dz/2; %Gauss points 
                    res1 = vol/(4*pi)/(norm([0 0 0]-xyz_trg));%Integrate_NumNum_K(xyz_trg,d,[0 0 0].');
                    xyz_trg(:,3) = r_m(3)-dz/2; %Gauss points 
                    res2 = vol/(4*pi)/(norm([0 0 0]-xyz_trg));
                    Kxy = (-res1+res2)...
                                  *dy... line integral
                                  *(1/(dx*dz))...% shape function in the source volume   
                                  *(dx/(dy*dz)); % projection;
                    K(mx,my,mz,3) = Kxy;                          
                    % Kxz
                    xyz_trg=zeros(1,3);
                    xyz_trg(:,1) = r_m(1); %Gauss points 
                    xyz_trg(:,2) = r_m(2)+dy/2; %Gauss points 
                    xyz_trg(:,3) = r_m(3); %Gauss points 
                    res1 = vol/(4*pi)/(norm([0 0 0]-xyz_trg));
                    xyz_trg(:,2) = r_m(2)-dy/2; %Gauss points 
                    res2 = vol/(4*pi)/(norm([0 0 0]-xyz_trg));
                    Kxz = (res1-res2)...
                                     *dz...
                                     *(1/(dx*dy))...% shape function in the source volume   
                                     *(dx/(dy*dz)); % projection;
                    K(mx,my,mz,2) = -Kxz;                          
                    % Kyz
                    xyz_trg=zeros(1,3);
                    xyz_trg(:,1) = r_m(1)+dx/2; %Gauss points 
                    xyz_trg(:,2) = r_m(2); %Gauss points 
                    xyz_trg(:,3) = r_m(3); %Gauss points 
                    res1 = vol/(4*pi)/(norm([0 0 0]-xyz_trg));
                    xyz_trg(:,1) = r_m(1)-dx/2; %Gauss points 
                    res2 = vol/(4*pi)/(norm([0 0 0]-xyz_trg));
                    Kyz = (-res1+res2)...
                                     *dz...
                                     *(1/(dx*dy))...% shape function in the source volume   
                                     *(dy/(dx*dz)); % projection
                    K(mx,my,mz,1) = Kyz;                         
                end
            end
        end
    %Remove unnecessary memory consuming elements
    clear xyz1 ww1 xyz2 ww2
    clear xyz_trg
end   
mytoc_k = toc(mytic_k);
disp([' Time for computing K Green tensor ::: ',num2str(mytoc_k)])
end

