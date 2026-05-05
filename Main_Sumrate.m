%==========================================================================
% Ying Wang, wangyingstu@163.com
% (c) 2025 by Ying Wang.   
%==========================================================================
clc; clear; close all;
%% Parameters
M_SNR_dB      = 0:10:40;  
NrRepetitions = 100; 

NrVerticalAnetnna    = 16;     % number of vertical anetnnas
NrHorizentalAntenna  = 8;      % number of horizental antennas
NrBSAntenna          = NrVerticalAnetnna*NrHorizentalAntenna;
NrUser               = 8;      % number of single antenna users
NrRF                 = 16;
Np                   = 3;      % Number of resolvable paths
Nrpilot              = 16;     % number of pilot blocks
ArrayGeometry        = 'UPA';
fc                   = 60e9;  % carrier frequency
for i_SNR=1:length(M_SNR_dB) 
    temp = 0; temp_MiAMP = 0; temp_OMP = 0; temp_AMP = 0; temp_SSD = 0;
    for iter = 1:NrRepetitions
        %% Lens Array Array
        % Saleh-Valenzuela channel model
        U = Channel.LensArray(NrVerticalAnetnna,NrHorizentalAntenna,0.5,ArrayGeometry);
        H=zeros(NrBSAntenna,NrUser); % the beamspace channel 
        for q=1:NrUser
            [H(:,q)] = U.'*Channel.Saleh_Valenzuelachannel(NrVerticalAnetnna,NrHorizentalAntenna,Np,fc,ArrayGeometry);
        end
        H_beam = H; 
        SNR_dB = M_SNR_dB(i_SNR); 
        Pn     = 10.^(-SNR_dB/10);
        alpha  = 2*Pn^(log10(exp(1)/2));
        W      = (2*randi([0,1],Nrpilot*NrRF,NrBSAntenna) - 1)/sqrt(Nrpilot*NrRF);
        noise  = sqrt(Pn/2)*(randn(NrBSAntenna,NrUser)+1j*randn(NrBSAntenna,NrUser));
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
        %%%%%%% perfect CSI, IA-BS
        [Hr,Hr_e] = Estimator.IA_BS(H_beam,H_beam); 
        F = Hr_e/(Hr_e'*Hr_e);
        beta = sqrt(NrUser/trace(F*F'));
        H_eq=Hr'*F;
        for k=1:NrUser
            sum_inf=sum(abs(H_eq(:,k)).^2)-abs(H_eq(k,k))^2;
            temp=temp+log2(1+abs(H_eq(k,k))^2/(sum_inf+NrUser/((1/Pn)*beta^2)));
        end

        %%%%%% OMP, IA-BS
        [Hr_OMP,Hr_OMP_e] = Estimator.IA_BS(H_hat_OMP,H_beam); 
        F = Hr_OMP_e/(Hr_OMP_e'*Hr_OMP_e);
        beta_OMP = sqrt(NrUser/trace(F*F'));
        H_eq=Hr_OMP'*F;
        for k=1:NrUser
            sum_inf=sum(abs(H_eq(:,k)).^2)-abs(H_eq(k,k))^2;
            temp_OMP=temp_OMP+log2(1+abs(H_eq(k,k))^2/(sum_inf+NrUser/((1/Pn)*beta_OMP^2)));
        end

        %%%%%% SSD, IA-BS
        [Hr_SSD,Hr_SSD_e] = Estimator.IA_BS(H_hat_SSD,H_beam); 
        F = Hr_SSD_e/(Hr_SSD_e'*Hr_SSD_e);
        beta_SSD = sqrt(NrUser/trace(F*F'));
        H_eq=Hr_SSD'*F;
        for k=1:NrUser
            sum_inf=sum(abs(H_eq(:,k)).^2)-abs(H_eq(k,k))^2;
            temp_SSD=temp_SSD+log2(1+abs(H_eq(k,k))^2/(sum_inf+NrUser/((1/Pn)*beta_SSD^2)));
        end

        %%%%%% MiAMP, IA-BS
        [Hr_MiAMP,Hr_MiAMP_e] = Estimator.IA_BS(H_hat_MiAMP,H_beam); 
        F = Hr_MiAMP_e/(Hr_MiAMP_e'*Hr_MiAMP_e);
        beta_MiAMP = sqrt(NrUser/trace(F*F'));
        H_eq1=Hr_MiAMP'*F;
        for k=1:NrUser
            sum_inf=sum(abs(H_eq1(:,k)).^2)-abs(H_eq1(k,k))^2;
            temp_MiAMP=temp_MiAMP+log2(1+abs(H_eq1(k,k))^2/(sum_inf+NrUser/((1/Pn)*beta_MiAMP^2)));
        end

        %%%%%% AMP1, IA-BS
        [Hr_AMP,Hr_AMP_e] = Estimator.IA_BS(H_hat_AMP,H_beam); 
        F = Hr_AMP_e/(Hr_AMP_e'*Hr_AMP_e);
        beta_AMP = sqrt(NrUser/trace(F*F'));
        H_eq1=Hr_AMP'*F;
        for k=1:NrUser
            sum_inf=sum(abs(H_eq1(:,k)).^2)-abs(H_eq1(k,k))^2;
            temp_AMP=temp_AMP+log2(1+abs(H_eq1(k,k))^2/(sum_inf+NrUser/((1/Pn)*beta_AMP^2)));
        end
    end
    C_PerCSI(i_SNR) = temp/NrRepetitions;
    C_SSD(i_SNR)    = temp_SSD/NrRepetitions;
    C_OMP(i_SNR)    = temp_OMP/NrRepetitions;
    C_AMP(i_SNR)    = temp_AMP/NrRepetitions;
    C_MiAMP(i_SNR)   = temp_MiAMP/NrRepetitions;

    fprintf('SNR = %.2f dB complete with avg. C_PerCSI = %.2f bits/s/Hz, C_Multi-d-AMP = %.2f bits/s/Hz, C_AMP = %.2f bits/s/Hz, C_0MP = %.2f bits/s/Hz, C_SSD = %.2f bits/s/Hz \n', ...
    M_SNR_dB(i_SNR), C_PerCSI(i_SNR),C_MiAMP(i_SNR),C_AMP(i_SNR), C_OMP(i_SNR), C_SSD(i_SNR));
end
%% Plot Results
figure(1)
plot(M_SNR_dB,C_PerCSI,'-o','Color',[0 0 0]*0.90,'Linewidth',1.5,'MarkerSize',10);grid on; hold on;
plot(M_SNR_dB,C_MiAMP, '-+','Color',[1 0 0]*0.90,'Linewidth',1.5,'MarkerSize',10);
plot(M_SNR_dB,C_AMP,   '-^','Color',[1 0 1]*0.90,'Linewidth',1.5,'MarkerSize',10);
plot(M_SNR_dB,C_OMP,   '-s','Color',[0 1 0]*0.80,'Linewidth',1.5,'MarkerSize',10);
plot(M_SNR_dB,C_SSD,   '-*','Color',[0 0 1]*0.90,'Linewidth',1.5,'MarkerSize',10);

xlabel('SNR, $P_{s}/P_{n}$ (dB)','Interpreter','latex');
ylabel('Achievable sum-rate (bits/s/Hz)');
L1=legend('Perfect CSI','Multi-d-AMP-estimated CSI','AMP-estimated CSI','OMP-estimated CSI','SSD-estimated CSI','Location','northwest');
set(gca,'FontSize',16, 'FontName','Times New Roman','LooseInset', [0,0,0,0])
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
