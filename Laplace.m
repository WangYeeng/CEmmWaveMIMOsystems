%==========================================================================
% Ying Wang, wangyingstu@163.com
% (c) 2025 by Ying Wang.   
%==========================================================================
clear;clc;close all;rng(1.5)
%% Parameters
NrVerticalAnetnna    = 16;     % number of vertical anetnnas
NrHorizentalAntenna  = 8;      % number of horizental antennas
NrBSAntenna          = NrVerticalAnetnna*NrHorizentalAntenna;
NrUser               = 8;     % number of single antenna users
Np                   = 3;     % Number of resolvable paths
ArrayGeometry        = 'UPA';
fc                   = 300e9; % carrier frequency
% Parameters of the Laplace distribution
mu = 66;          % Position parameters
b  = 1;      % Dimensional parameters
%% Lens Array Array
% Saleh-Valenzuela channel model
U = Channel.LensArray(NrVerticalAnetnna,NrHorizentalAntenna,0.5,ArrayGeometry);
H = zeros(NrBSAntenna,NrUser); % the beamspace channel 
for q=1:NrUser
  [H(:,q)] = U.'*Channel.Saleh_Valenzuelachannel(NrVerticalAnetnna,NrHorizentalAntenna,Np,fc,ArrayGeometry);
end

x_vals = linspace(0, 128, 1000);
pdf_vals1 = (2*b/(2)) * exp(-abs(x_vals - mu)/b);
%% Plot
figure()
plot(x_vals, pdf_vals1,"LineWidth",1.2);hold on;grid on;
legend({"$\sigma = 2,\ {{m}} = 66$"},'Interpreter','latex');
xlabel('${{m}}$','Interpreter','latex'); ylabel('$p\left( h_m \right)$','Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',16,'LooseInset', [0,0,0,0]);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
figure()
plot(abs(H(:,1)),"LineWidth",1.2);hold on;grid on;
plot(abs(H(:,3)),"LineWidth",1.2);
plot(abs(H(:,6)),"LineWidth",1.2);
xlabel('${{m}}$','Interpreter','latex'); ylabel('$|{\bf h}_q|$','Interpreter','latex');
legend({"UE#1","UE#3","UE#6"});
set(gca,'FontName','Times New Roman','FontSize',16,'LooseInset', [0,0,0,0]);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);



