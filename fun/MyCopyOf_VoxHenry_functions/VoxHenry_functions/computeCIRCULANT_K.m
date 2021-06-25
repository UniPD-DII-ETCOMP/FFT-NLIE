function [Circ] = computeCIRCULANT_K(K)
type = 'FFTofC'; %default FFT of Circulant operator
mytic_circ_K = tic; 
[Kp_mn] = circulant_kop(K);
mytoc_circ_K = toc(mytic_circ_K);
disp([' Time for computing circulant tensor ::: ',num2str(mytoc_circ_K)])
if strcmp(type ,'FFTofC')
    disp(' Computing FFT of circulant tensor');
    mytic_fft_K = tic;
    Circ = fft_operator(Kp_mn);
    mytoc_fft_K = toc(mytic_fft_K);
    disp([' Time for FFT of circulant tensor ::: ',num2str(mytoc_fft_K)])
else
    Circ = Kp_mn;
    infofN = whos('Circ'); memestimated = 2*infofN.bytes/(1024*1024);
    disp('  Two circulant tensors will be held in memory for fast sweep!');
    disp(['  Memory for storing two circulant tensors (MB) ::: ' , num2str(memestimated)]);
end
end