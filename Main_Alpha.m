%==========================================================================
% Ying Wang, wangyingstu@163.com
% (c) 2025 by Ying Wang.   
%==========================================================================
clear;clc;close all;
%% Parameters
NrRepetitions       = 100;               % Number of Monte Carlo repetitions
M_SNR_dB            = [0,10,20];         % Signal-to-noise ratio

NrVerticalAnetnna    = 16;     % number of vertical anetnnas
NrHorizentalAntenna  = 8;     % number of horizental antennas
NrBSAntenna          = NrVerticalAnetnna*NrHorizentalAntenna;
NrUser               = 8;     % number of single antenna users
NrRF                 = 16;
Np                   = 3;     % Number of resolvable paths
Nrpilot              = 16;    % number of pilot blocks
ArrayGeometry        = 'UPA';
fc                   = 300e9; % carrier frequency
tic;
for i_rep = 1:NrRepetitions     
%% Lens Array Array
% Saleh-Valenzuela channel model
U = Channel.LensArray(NrVerticalAnetnna,NrHorizentalAntenna,0.5,ArrayGeometry);
H=zeros(NrBSAntenna,NrUser); % the beamspace channel 
for q=1:NrUser
  [H(:,q)] = U.'*Channel.Saleh_Valenzuelachannel(NrVerticalAnetnna,NrHorizentalAntenna,Np,fc,ArrayGeometry);
end
alpha = 0.5:0.025:3;
for i_alpha = 1:length(alpha)
    for i_SNR = 1:length(M_SNR_dB)
        SNR_dB = M_SNR_dB(i_SNR); 
        Pn     = 10.^(-SNR_dB/10);
        W      = (2*randi([0,1],Nrpilot*NrRF,NrBSAntenna) - 1)/sqrt(Nrpilot*NrRF);
        noise  = sqrt(Pn/2)*(randn(NrBSAntenna,NrUser )+1j*randn(NrBSAntenna,NrUser ));
        y      = W*(H + noise);
        [H_hat_MiAMP] = Estimator.MiAMP(y,W,alpha(i_alpha));
        NMSE_MiAMP(i_alpha,i_rep,i_SNR) = 10*log10((norm(H_hat_MiAMP-H,"fro"))^2/(norm(H,"fro"))^2);
    end
end
TimePassed = toc;
if mod(i_rep,10)==0
disp(['Realization ' int2str(i_rep) ' of ' int2str(NrRepetitions) '. Time left: ' int2str(TimePassed/i_rep*(NrRepetitions-i_rep)/60) 'minutes']);
end
end
%% Plot Results
MarkerSize = 10;
Linewidth  = 1.5;
[Nmes_0db,index_0db]  =min(mean(NMSE_MiAMP(:,:,1),2));
[Nmes_10db,index_10db]=min(mean(NMSE_MiAMP(:,:,2),2));
[Nmes_20db,index_20db]=min(mean(NMSE_MiAMP(:,:,3),2));
figure(1);
b1 = plot(alpha,mean(NMSE_MiAMP(:,:,1),2),      '-','Color',[1 0 0]*0.75,'MarkerSize',MarkerSize,'Linewidth',Linewidth);hold on;grid on;
plot([alpha(index_0db),alpha(index_0db)],  [Nmes_0db-2.5,Nmes_0db+2.5],'Color',[1 1 1]*0.25,'Linewidth',1.5);
b2 = plot(alpha,mean(NMSE_MiAMP(:,:,2),2),      '-','Color',[0 1 0]*0.75,'MarkerSize',MarkerSize,'Linewidth',Linewidth);
plot([alpha(index_10db),alpha(index_10db)],[Nmes_10db-2.5,Nmes_10db+2.5],'Color',[1 1 1]*0.25,'Linewidth',1.5);
b3 = plot(alpha,mean(NMSE_MiAMP(:,:,3),2),      '-','Color',[0 0 1]*0.75,'MarkerSize',MarkerSize,'Linewidth',Linewidth);
plot([alpha(index_20db),alpha(index_20db)],[Nmes_20db-2.5,Nmes_20db+2.5],'Color',[1 1 1]*0.25,'Linewidth',1.5);
xlabel('$\varsigma $','Interpreter','latex');
ylabel('NMSE (dB)');
legend([b1,b2,b3],{'SNR =  0 dB','SNR = 10dB','SNR = 20dB'},'Location','southeast');
set(gca,'FontName','Times New Roman','FontSize',16,'LooseInset', [0,0,0,0]);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);

