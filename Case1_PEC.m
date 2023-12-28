%% Define frequencies range and solving parameters
fbeg=0.5;
fstp=0.05;
fend=7;
freq=fbeg:fstp:fend;% unit: GHz
a=16.6667e-3;% radius, unit: mm
k=2*pi*freq*1e9/3e8;% wave number
deg=17;% degree of SVWFs
%% Calculate translation of SVWFs
ro1=50e-3;% the center of sphere 1
ro2=100e-3;% the center of sphere 2. ro2 > ro1;
ParaList=[ro1 ro2 deg fbeg fend fstp];
save('TransInfo.txt','ParaList','-ascii');%Interface Files
% system('Trans_Solver.exe');
%% sweep frequency
for i=1:length(freq)
    f=freq(i)
    % calculate the T matrix for single sphere
    T1=Functions.Tsph_PEC(deg,f,a);
    T2=T1;% the two sphere are identical
    % load translation matrices
    R1=Functions.ReadBinMat(['TranlationMatrix\' num2str(f*1e9),'.000000R1.dat']);
    R2=Functions.ReadBinMat(['TranlationMatrix\' num2str(f*1e9),'.000000R2.dat']); 
    Yr=Functions.ReadBinMat(['TranlationMatrix\' num2str(f*1e9),'.000000Yr.dat']);
    Yi=Functions.ReadBinMat(['TranlationMatrix\' num2str(f*1e9),'.000000Yi.dat']);
    Y=Yr+1j*Yi; 
    I=eye(size(R1));
    
    T=R1*inv(I-T1*Y*T2*Y.')*(T1*R1'+T1*Y*T2*R2')...
        +R2*inv(I-T2*Y.'*T1*Y)*(T2*R2'+T2*Y.'*T1*R1');% total T matrix of the two spheres
    
    S=I+2*T;
    S0=I+2*R2*T2*R2.';% set sphere 2 as background
    
    [F,Val]=eigs((S*S0'-I)/2,40);% extract the first 40 eigenvectors
    Val=diag(Val);
    save(['Case1_Data\' num2str(f*1e9),'.000000result.mat'],'F','Val');
end
%% Postprocess
Functions.Tracking('Case1_Data',fbeg,fend,fstp,a);