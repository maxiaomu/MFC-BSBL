% if the imaging result is poor,please check for the delta_v (+ -)
clear;clc;

% Load data 
% Y = J * Phi * X + V
% Y;M * L    
% J;M * N    
% Phi;N * GH 
% X: GH*L    
% V;M * L    

% load J_sen;     %
% load deltav_danfang1;
% load deltav11
% delta_v = deltav';  % move_ele

% load 'simulation_data/delta_vhvi_sigma-20_1fang_J.mat'; use_true =1;%simulation
load delta_vhvi_sigma-80_1fang_J.mat
use_true = 0; 
load J_sen.mat
M=8; 
N=6272;
h = 4; % overlap
g = N/h; % number of blocks

Phi_ = J;
% Phi = J;
% Phi_ = zscore(Phi);

groupStartLoc = 1:h:N;

%% run....
tic
%% ==== Compress signal ======
% Y = Phi_ * X_ /4.0;
Frate = 91;
% Y_ = delta_v(:,Frate:Frate+M-1);
% Y = zscore(Y_);
Y = delta_v(:,Frate:Frate+M-1);

%%
% Y1 = delta_v(:,1:11)';%:2
% % X= sqrt(2)*sin(0:pi/1000000:6*pi);                % generate sin
% X1=Y1;
% Y2 = awgn(X1,55,'measured');                          % generate noise 
% sigPower = sum(abs(X1).^2)/length(X1);            % sigpower
% noisePower = sum(abs(Y2-X1).^2)/length(Y2-X1);   % noisepower
% SNR = 10*log10(sigPower/noisePower);          % add noise db
% Y = Y2'; %  
  
Result = MFC_BSBL(Phi_,Y,groupStartLoc,2,'LEARNTYPE', 1,'prune_gamma',1e-3,'max_iters',8,'EPSILON',1e-5);
signal_hat = Result.x; 
if use_true == 1
    signal_hat(signal_hat<0) = 0;
    x = signal_hat;
else
    x = signal_hat;
end
toc
save recover_x x; 
       
 