%-- Code written by Aditya Varma Muppala for the paper titled: Fast-Fourier Time-Domain SAR Reconstruction for
%-- Millimeter-Wave FMCW 3-D Imaging. Last edited on 03/05/2024. 

clc
clear
clear path

%-- Radar System Operational Parameters
fBW = 8e9;                 % bandwidth
fc = 79.6e9;                 % carrier frequency
wc = 2*pi*fc;
c  = 3e8;                  % RF propagation speed
theta_b = 60;              % Antenna Beamwidth in degrees

wBW = 2*pi*fBW;
RR = c/(2*fBW);

%-- Fast-Time domain parameters and arrays
fs = 50*1e6;                % Sampling rate;
T = 48*1e-6;                % in microseconds. Last 2 microseconds are dead zone

N = fs * T;                 % number of fast-time samples
gamma = wBW/T;              % Chirp rate

dt = 1/fs;
Tp  = (N-1)*dt;             % Chirp pulse duration
t = linspace(-Tp/2,Tp/2,N); % time array for data acquisition

range = 0:RR:RR*(N-1);

%-- Slow-Time domain parameters and arrays
spacing = 1.6e-3;
xp = -200e-3:spacing:200e-3;
yp = -200e-3:spacing:200e-3;
zp = 0;

Mx = length(xp);           
My = length(yp);           
M = Mx*My;                 % Number of slow time measurements

% Raw data acquisition (Try different data sets available in Imaging_new_data folder, but change the background file as well.)
fileID = fopen('Imaging_raw_data/SiemensStencil_Measurement_1.ch0','rb');
ch1 = fread(fileID,'double');
fclose(fileID);

fileID2 = fopen('Imaging_raw_data/Bkg_SiemensStencil_Measurement_1.ch0','rb');
ch1_bkg = fread(fileID2,'double');
fclose(fileID2);

% Re-shaping to fast time X slow time array
valid_samples_per_trigger = fs*T;
ch1b = squeeze(reshape(ch1,valid_samples_per_trigger,[],M));

sdata_time = ch1b-ch1_bkg; % Background subtraction

% sdata_time = sdata_time.*hanning(valid_samples_per_trigger); % Windowing in time (not required)

%% Calibration step 2
% Non-linear phase calibration is applied during data collection and is included in the raw data.
% Doing it here would make things too slow so it is better to apply it during data acquisition.

% Phase center calibration
del_lt = 19.1e-2; % Obtained by sweeping
exp_shift = transpose(exp(-1i*(wc*2*del_lt/c+2*gamma*del_lt*t/c)));
sdata_time2 = sdata_time.*exp_shift;

% %% Reformatting to complex exponential (Not required)
% % s1f = fft(sdata_time,[],1);
% % s1f2 = ifft(s1f(1:length(s1f)/2,:),[],1);
% % sdata_time = resample(s1f2,2,1);

%% Rearranging raster scanned data to XY - grid

disp('Reformatting raster scan to XY grid');
F1 = zeros(Mx,My,N);
main_ind = 1;

for Myind = My:-1:1
    if rem(Myind,2) == 0
        for Mxind = 1:Mx
            F1(Mxind,Myind,:) = sdata_time2(:,main_ind);
            main_ind = main_ind+1;
        end
    else
        for Mxind = Mx:-1:1
            F1(Mxind,Myind,:) = sdata_time2(:,main_ind);
            main_ind = main_ind+1;
        end
    end
end

%% Reconstruction

disp('Starting algorithm. Est. time: 40 seconds');
disp('First 2-D FFT');

Kr = -2*gamma*t/c+2*wc/c; % Negative sign since downchirp is used
k_limit = (2*wc/c)*sind(theta_b/2); % Light cone boundary

dkx = 2*pi/(2*max(xp));
dky = 2*pi/(2*max(yp));
kx = dkx*(-(Mx-1)/2:(Mx-1)/2);
ky = dky*(-(My-1)/2:(My-1)/2);

% Truncating data that lies inside antenna beam light cone
kx_trunc_ind = abs(kx)<=k_limit;
ky_trunc_ind = abs(kx)<=k_limit;
kx_trunc = kx(kx_trunc_ind);
ky_trunc = ky(ky_trunc_ind);

Mxn = length(kx_trunc); Myn = length(ky_trunc);

kx3D = repmat(kx_trunc.',[1,Myn,N]);
ky3D = repmat(ky_trunc,[Mxn,1,N]);
Qt1D(1,1,:) = Kr;
Qt3D = repmat(Qt1D,[Mxn,Myn,1]);
kz3D = sqrt(Qt3D.^2-kx3D.^2-ky3D.^2);

kz1D_uni = kz3D((Mxn+1)/2,(Myn+1)/2,:); % Uniform sampling of Kz taken along Kz axis

indX_range = 1:Mx; indY_range = 1:Mx;
truncX_range = indX_range(kx_trunc_ind);
truncY_range = indY_range(ky_trunc_ind);

tic
F2 = fty(ftx(F1)); % First 2D FFT

%% Non-uniform to Uniform sampling and Range FFT
% This section can be sped up by a 1-D NUFFT in C using a MEX file. See 1-D NUFFT by Flatiron institute.

disp('NUFFT: Resampling to uniform grid and Range FFT');

fhat = zeros(Mxn,Myn,N);

for indX = 1:Mxn
%     indX
    for indY = 1:Myn
        fhat(indX,indY,:) = interp1(squeeze(kz3D(indX,indY,:)),squeeze(F2(truncX_range(indX),truncY_range(indY),:)),kz1D_uni);
    end
end

fhat(isnan(fhat))=0; % Truncating extrapolated values lying outside the grid to zero.
fhat_F = fft(fhat,[],3); % 1-D IFFT in Z. FFT is used due to sign convention in t.

%% Image gen
% Zero padding improves image quality when using shading interp. It is not required though.
% Can remove it to speed up the reconstruction.

disp('Last 2-D IFFT with zero-padding'); 

Mxn2 = 2^(nextpow2(Mxn));
pad_size = floor((Mxn2-Mxn)/2);
fhat2 = padarray(fhat_F,[pad_size pad_size],0,'both');
szNew = size(fhat2);
Mxn3 = szNew(1); Myn3 = szNew(2);

fhat3 = ifty(iftx(fhat2));
time = toc;

disp(['Total computation time: ',num2str(time), 'seconds']);

dx = 2*pi/(max(2*kx_trunc));
dy = 2*pi/(max(2*ky_trunc));
dz = 2*pi/(max(kz3D((Mxn+1)/2,(Myn+1)/2,:))-min(kz3D((Mxn+1)/2,(Myn+1)/2,:)));

xIm = dx*(-(Mxn3-1)/2:(Mxn3-1)/2);
yIm = dy*(-(Myn3-1)/2:(Myn3-1)/2);
distZ = dz*(1:N);

[X,Y] = ndgrid(xIm,yIm);

%% Plotting Final image

f = sum(abs(fhat3(:,:,20)),3);
% f = circshift(f,-12,1);
% f = circshift(f,-10,2);
figure(1); clf
surf(X*100,Y*100,abs(f/max(f(:))));
% surf(X*100,Y*100,20*log10(abs(f/max(f(:)))));
shading interp;
xlabel('X (cm)');
ylabel('Y (cm)');
% title('f(x,y) - Target Scene');
hold on;
axis equal; axis tight
colormap gray;
view([0 90])
colorbar
xlim([-15 15])
ylim([-15 15])
% set(gca, 'clim', [-50 0]);
set(gca, 'clim', [0 0.9]);
set(gcf,'position',[550,50,800,800])
set(gca,'FontSize',20)
set(gca,'GridAlpha',1)
set(gca,'GridLineStyle','--')
set(gca,'fontname','ariel')
set(gcf,'color','w');

