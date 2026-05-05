%==========================================================================
% Ying Wang, wangyingstu@163.com
% (c) 2025 by Ying Wang.   
%==========================================================================
clear;clc;close all;
%% Parameters
NrRepetitions       = 100;               % Number of Monte Carlo repetitions
M_SNR_dB            = 0:5:30;           % Signal-to-noise ratio

NrVerticalAnetnna    = 16;     % number of vertical anetnnas
NrHorizentalAntenna  = 8;      % number of horizental antennas
NrBSAntenna          = NrVerticalAnetnna*NrHorizentalAntenna;
NrUser               = 8;     % number of single antenna users
NrRF                 = 16;
Np                   = 3;     % Number of resolvable paths
Nrpilot              = 16;    % number of pilot blocks
ArrayGeometry        = 'UPA';
fc                   = 60e9; % carrier frequency
tic;
for i_rep = 1:NrRepetitions     
%% Lens Array Array
% Saleh-Valenzuela channel model
U = Channel.LensArray(NrVerticalAnetnna,NrHorizentalAntenna,0.5,ArrayGeometry);
H = zeros(NrBSAntenna,NrUser); % the beamspace channel 
for q=1:NrUser
  [H(:,q)] = U.'*Channel.Saleh_Valenzuelachannel(NrVerticalAnetnna,NrHorizentalAntenna,Np,fc,ArrayGeometry);
end
for i_SNR = 1:length(M_SNR_dB)
    SNR_dB = M_SNR_dB(i_SNR); 
    Pn     = 10.^(-SNR_dB/10);
    alpha  = 2*Pn^(log10(exp(1)/2));
    W      = (2*randi([0,1],Nrpilot*NrRF,NrBSAntenna) - 1)/sqrt(Nrpilot*NrRF);
    noise  = sqrt(Pn/2)*(randn(NrBSAntenna,NrUser )+1j*randn(NrBSAntenna,NrUser ));
    y      = W*(H + noise);
    [H_hat_MiAMP] = Estimator.MiAMP(y,W,alpha);
    H_hat_OMP=[];
    H_hat_SSD=[];
    H_hat_AMP=[];
    for i=1:NrUser
        yi=y(:,i);
        h_hat_OMP = Estimator.OMP(yi, W, 40);
        h_hat_SSD = Estimator.SSD(yi, W, Np, NrBSAntenna, 1, 1e9, fc, 40);
        h_hat_AMP = Estimator.AMP(yi, W, 1.2);
        H_hat_OMP =[H_hat_OMP,h_hat_OMP];
        H_hat_SSD =[H_hat_SSD,h_hat_SSD];
        H_hat_AMP =[H_hat_AMP,h_hat_AMP]; 
    end
    NMSE_MiAMP(i_SNR,i_rep) = 10*log10((norm(H_hat_MiAMP-H,"fro"))^2/(norm(H,"fro"))^2);
    NMSE_OMP(i_SNR,i_rep)   = 10*log10((norm(H_hat_OMP-H, "fro"))^2/(norm(H,"fro"))^2);
    NMSE_SSD(i_SNR,i_rep)   = 10*log10((norm(H_hat_SSD-H, "fro"))^2/(norm(H,"fro"))^2);
    NMSE_AMP(i_SNR,i_rep)   = 10*log10((norm(H_hat_AMP-H, "fro"))^2/(norm(H,"fro"))^2);
    NMSE_IAM(i_SNR,i_rep)   = (norm(H_hat_OMP-H, "fro"))^2/(norm(H,"fro"))^2;
end
TimePassed = toc;
if mod(i_rep,10)==0
disp(['Realization ' int2str(i_rep) ' of ' int2str(NrRepetitions) '. Time left: ' int2str(TimePassed/i_rep*(NrRepetitions-i_rep)/60) 'minutes']);
end
end
%% Plot Results
MarkerSize = 10;
Linewidth  = 1.5;
figure(1);
plot(M_SNR_dB,mean(NMSE_SSD,2),      '-*','Color',[0 1 0]*0.75,'MarkerSize',MarkerSize,'Linewidth',Linewidth);
hold on;grid on;
plot(M_SNR_dB,mean(NMSE_OMP,2),      '-x','Color',[1 0 1]*0.75,'MarkerSize',MarkerSize,'Linewidth',Linewidth);
plot(M_SNR_dB,mean(NMSE_AMP,2),      '-^','Color',[0 0 1]*0.75,'MarkerSize',MarkerSize,'Linewidth',Linewidth);
plot(M_SNR_dB,mean(NMSE_MiAMP,2),    '-o','Color',[1 0 0]*0.75,'MarkerSize',MarkerSize,'Linewidth',Linewidth);
xlabel('SNR, $P_{s}/P_{n}$ (dB)','Interpreter','latex');
ylabel('NMSE (dB)');
legend('SSD','OMP','AMP','Multi-d-AMP','Location','southwest');
set(gca,'FontName','Times New Roman','FontSize',16,'LooseInset', [0,0,0,0]);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);

