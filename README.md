<img src="docs/ind.png" width="220" height="200"><img src="docs/JM.png" width="200" height="200">  

# FFT-Non_Linear_Inductance_Extractor 

This directory contains a FFT-Non_Linear_Inductance_Extractor code for extracting inducatnce of non-linear inductors 

This code has been inspired by VoxHenry https://github.com/acyucel/VoxHenry

-------------------------------------------------------------------

# Description
 
FFT_NLIE.m is the main file you must run to start the code. 

All user-settable quantities are contained in the block identified by the 
BEGIN USER SETTINGS / END USER SETTINGS comments.

Available test cases
--------------------
Three simple test cases are contained in separate directories under "data". 
Set the "name_dir" variable in "FFT_NLIE.m"  to the appropriate directory.

User-defined test cases
-----------------------
Follow the instuctions given in "README.txt" inside the "data" directory.

Results
--------------------
Th inductance value is shown in the Command Window.  

Credits
--------------------
If you use FFT_NLIE, please consider citing:

 [1] [R. Torchio, "A Volume PEEC Formulation Based on the Cell Method for Electromagnetic Problems From Low to High Frequency," in IEEE Transactions on Antennas and Propagation, vol. 67, no. 12, pp. 7452-7465, Dec. 2019, doi: 10.1109/TAP.2019.2927789](https://ieeexplore.ieee.org/document/8764572)

 [2] [Bettini et al., "Fast Fourier transform-volume integral: a smart approach for the electromagnetic design of complex systems in large fusion devices", Plasma Physics and Controlled Fusion, Volume 63, Number 2, doi: 10.1088/1361-6587/abce8f](https://iopscience.iop.org/article/10.1088/1361-6587/abce8f)
 
 [3] R. Torchio et al., "FFT-PEEC: A Fast Tool From CAD to Power Electronics Simulations," in IEEE Transactions on Power Electronics, doi: 10.1109/TPEL.2021.3092431
 
and FFT-FFT_NLIE itself

 [4] R. Torchio, "FFT-FFT_NLIE toolbox", https://github.com/UniPD-DII-ETCOMP/FFT-FFT_NLIE
 
Contacts & Authors
-----------------------
For any questions or need help to solve your problem, please contact us

Riccardo Torchio (riccardo.torchio@unipd.it)
