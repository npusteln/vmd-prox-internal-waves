% demo_VMDprox.m
%
% E. Horne, J. Schmitt, N. Pustelnik, S. Joubaud. P. Odier
% Variational Mode Decomposition for estimating critical reflected internal wave in stratified fluid,
% submitted, 2020. 
%
% J. Schmitt. Version: 25-Mars-2021.

clear all
close all
addpath(genpath('.'));

%% Load example of experimental image and experimental parameters
load('experiment/U.mat')
load('experiment/parameters_experiment.mat')
% N0:       buoyancy frequency s^{_1}
% lambda:   vertical wavelength
% gamma:    angle of the slope respect to horizontal 
% gammacm:  coeficient of mapping pixel to cm 
% piv_grid: resolution in pixels of the PIV method 

%% Creates synthetic image based on experiment and D-Y model
[Usy,Usy_inc,Usy_refl]=synthetic_critical_reflected_image(N0,lambda,gamma,gammacm,U_inc,size(U,1),piv_grid,size(U,2),size(U,1)); % in cm/s and cm

%% Add noise to the synthetic field
noise_coef = 0.1; 
Usy  =  Usy + noise_coef*std(Usy(:))*randn(size(Usy));

%% Adding zero padding
tmp = zeros(2*size(Usy));
tmp1 = zeros(2*size(Usy));
tmp2 = zeros(2*size(Usy));
tmp(round(size(Usy,1)-size(Usy,1)/2)+1:round(size(Usy,1)+size(Usy,1)/2),round(size(Usy,2)-size(Usy,2)/2)+1:round(size(Usy,2)+size(Usy,2)/2)) = Usy;
tmp1(round(size(Usy_inc,1)-size(Usy_inc,1)/2)+1:round(size(Usy_inc,1)+size(Usy_inc,1)/2),round(size(Usy_inc,2)-size(Usy_inc,2)/2)+1:round(size(Usy_inc,2)+size(Usy_inc,2)/2)) = Usy_inc;
tmp2(round(size(Usy_refl,1)-size(Usy_refl,1)/2)+1:round(size(Usy_refl,1)+size(Usy_refl,1)/2),round(size(Usy_refl,2)-size(Usy_refl,2)/2)+1:round(size(Usy_refl,2)+size(Usy_refl,2)/2)) = Usy_refl;

Usyn = tmp;
Usyn_inc = tmp1;
Usyn_refl = tmp2;

%% VMD decomposition to synthetic field using and optimizing VMD parameters
u = Usyn;

rho = 10;
eta = 10;
beta = 10;
K = 2;
N = 1000;
alpha0 = 100;
alpha_z2 = 10;
Alpha = [alpha0 alpha0; alpha0 alpha_z2];
tol = 2e-06;

omega   = initialize_omega(0,N,K);           % Initialize the spatial frequencies around the unit circle. Chose init=1 to initialize frequencies randomly
[Uvmd, Uvmd_hat, omega_vmd, crit_vmd] = VMD_2D_prox_proj_zero(u, Alpha, rho, eta, beta, K, omega, tol, N);

%% Display results
figure,
subplot(3,2,1),imagesc(Usyn),title('U_{syn}= U^{inc}_{syn}+U^{refl}_{syn}')
subplot(3,2,3),imagesc(Usyn_inc),title('U^{inc}_{syn}')
subplot(3,2,5),imagesc(Usyn_refl),title('U^{refl}_{syn}')
subplot(3,2,2),imagesc(Uvmd(:,:,1)+Uvmd(:,:,2)),title('U_{VMD}= mode 1 + mode 2')
subplot(3,2,4),imagesc(Uvmd(:,:,1)),title('mode 1')
subplot(3,2,6),imagesc(Uvmd(:,:,2)),title('mode 2')

subplots = get(gcf,'Children'); % Get each subplot in the figure
for i=1:length(subplots) % for each subplot
    caxis(subplots(i),[-0.1,0.1]); % set the clim
end


figure, hold on
plot(mean(abs(Usyn),2),'--k','linewidth',2)
plot(mean(abs(Usyn_inc),2),'--b','linewidth',2)
plot(mean(abs(Usyn_refl),2),'--r','linewidth',2)
plot(mean(abs(Uvmd(:,:,1)+Uvmd(:,:,2)),2),'k')
plot(mean(abs(Uvmd(:,:,1)),2),'b')
plot(mean(abs(Uvmd(:,:,2)),2),'r')
set(gca, 'XDir','reverse')
legend('U^{inc}_{syn}+U^{refl}_{syn}','U^{inc}_{syn}','U^{refl}_{syn}','mode 1 + mode 2','mode 1','mode 2')
xlabel('z'),ylabel('<|u|>_{x}')

 %% Computes the SNR 
snr_inc = 20*log10(norm(Usyn_inc(1:size(Uvmd,1),:),'fro')/norm(Uvmd(:,:,1)-Usyn_inc(1:size(Uvmd,1),:),'fro'));
snr_refl = 20*log10(norm(Usyn_refl(1:size(Uvmd,1),:),'fro')/norm(Uvmd(:,:,2)-Usyn_refl(1:size(Uvmd,1),:),'fro'));




