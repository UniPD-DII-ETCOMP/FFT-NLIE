-------------------------------------------------------------------------------------------------------
How to create a new user-defined test-case:

1. duplicate the directory test1 and rename it (e.g. "user_test")
2. copy the user's files "stl1.stl", "stl2.stl",  etc inside "user_test"
3. customize the simulation data between %% BEGIN USER SETTINGS and %% END USER SETTINGS
    Description/Example:
    %% BEGIN USER SETTINGS
    paraview_export_flag = 1; % if 1 export paraview files in "res_para" directory
    x_ray_flag = 1;           % if 1 make x_ray plot
    model_name='pcb';         % name of the model
    %
    stl_files(1).name = 'pcb_coil.stl'; % name of the first stl file to load
    stl_files(1).tag = 'cond';          % tag for the material (write 'cond' for condutive media, "port" if you want to impose the current, "mag" for magnetic media )
    stl_files(1).cur=[];                % injected current value, only active if stl_files(1).tag='port';
    stl_files(1).rho=1/57e6;            % resistivity of the medium
    stl_files(1).BHcurve =[];           % BH curve, only active if stl_files(1).tag="mag"; 
    %
    stl_files(2).name = 'pcb_port1.stl'; % name of the second stl file to load
    stl_files(2).tag = 'port';
    stl_files(2).cur=1; % current value in ampere
    stl_files(2).rho=1/57e6;
    stl_files(1).BHcurve =[];
    %
    stl_files(3).name = 'pcb_port2.stl';
    stl_files(3).tag = 'port';
    stl_files(3).cur=-1;
    stl_files(3).rho=1/57e6;
    stl_files(1).BHcurve =[];
    %
    stl_files(4).name = 'test_core.stl'; 
    stl_files(4).tag = 'mag';
    stl_files(4).cur=[];
    stl_files(4).rho=[];%1/57e6;
    stl_files(4).BHcurve=  [0	        0
                        663.146	    1
                        1067.5	    1.1
                        1705.23	    1.2
                        2463.11	    1.3
                        3841.67	    1.4
                        5425.74	    1.5
                        7957.75	    1.6
                        12298.3	    1.7
                        20462.8	    1.8
                        32169.6	    1.9
                        61213.4	    2
                        111408	    2.1
                        188487.757	2.2
                        267930.364	2.3
                        347507.836	2.4];
    % to scale a stl file from any unit to meters
    scal_geomery.x=1/1000; scal_geomery.y=1/1000; scal_geomery.z=1/1000; % if you want to scale the dimension of the stl data
    % Box 
    % number of voxels in the x y z directions
    Nx=100; number of voxels in the x direction
    Ny=100; number of voxels in the y direction
    Nz=8;   number of voxels in the z direction
    % corners
    flag_auto=1; % if 1, user_data below are ignored and the box is created automatically (suggested)
    % user_data   corners of the Box
    meshXmin = -6e-3;      % (m)        
    meshXmax = 4.5e-3;     % (m)
    meshYmin = -4e-3;      % (m)
    meshYmax = 3.5e-3;     % (m)
    meshZmin = 0;          % (m)
    meshZmax = 0.4e-3;     % (m)
    %% END USER SETTINGS
3. run "mainVOXELISE.m"

File "data.mat" is now created and you can select "user_dir" in "FFT_PEEC_COND.m" 
to select the user test case
-------------------------------------------------------------------------------------------------------
