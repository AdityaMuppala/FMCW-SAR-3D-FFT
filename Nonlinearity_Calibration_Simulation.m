%% Phase Calibration simulation example

% Initialize Workspace
clc; clear;% close all;

%-- Radar System Operational Parameters
c  = 3e8;                  % RF propagation speed
fBW = 8e9;                 % bandwidth
fc = 79.6e9;                 %%%%%%%%%% carrier frequency

%-- Fast-Time domain parameters and arrays
fs = 50e6;               % Sampling rate;
N = 2^11;                % number of fast-time samples

dt = 1/fs;                  % Time domain sampling (guard band factor 2)
Ts = dt;                 % start time of sampling
Tf = Ts + N*dt;            % end time of sampling
Tp  = Tf - Ts;                  % Chirp pulse
t = linspace(Ts,Tf,N);          % time array for data acquisition
alpha = fBW/Tp;                 % Chirp rate

RR = c/(2*fBW);
range = (0:N-1)*RR;

%% Target Parameters
rng("default")
curvedegree = 5;  % Polynomial component of error
poly_coeff = 1e3; % Polynomial component of error
cos_coeff = 1e1;  % Sinusoidal component of error
cos_freq = 3e5;   % Sinusoidal component of error

coeffs = polyfit(t, poly_coeff*rand(size(t)), curvedegree);
coeffs(6) = 0;

xg = 0; yg = 0; zg = 0; % Radar location
xn = 0; yn = 0; zn = 20; fn = 1; % Target location


lt = sqrt((xn-xg)^2 + (yn-yg)^2 + (zn-zg)^2); % Time reference starts from zero
Tau  = 2*lt/c;
td = t - Tau;


epsi_t = polyval(coeffs, t)+cos_coeff*cos(t*cos_freq);
epsi_td = polyval(coeffs, td)+cos_coeff*cos(td*cos_freq);

s1 = exp(-1i*2*pi*(fc*Tau+alpha*t*Tau)); % Ideal linear no RVP
s2 = exp(-1i*2*pi*(fc*Tau+alpha*t*Tau-0.5*alpha*Tau^2)); % Ideal linear with RVP
s3 = exp(-1i*2*pi*(fc*Tau+alpha*t*Tau-0.5*alpha*Tau^2+epsi_t-epsi_td)); % Non linear

%% Phase non-linearity estimation (Note that absolute target location is not used)

xn1 = 0; yn1 = 0; zn1 = 1+1e-3; fn1 = 1; % Calib target 1 location
xn2 = 0; yn2 = 0; zn2 = 1; fn2 = 1; % Calib target 2 location

lt1 = sqrt((xn1-xg)^2 + (yn1-yg)^2 + (zn1-zg)^2); % Time reference starts from zero
Tau1  = 2*lt1/c;
td1 = t - Tau1;

lt2 = sqrt((xn2-xg)^2 + (yn2-yg)^2 + (zn2-zg)^2); % Time reference starts from zero
Tau2  = 2*lt2/c;
td2 = t - Tau2;

epsi_td1 = (polyval(coeffs, td1))+(cos_coeff*cos(td1*cos_freq));
epsi_td2 = (polyval(coeffs, td2))+(cos_coeff*cos(td2*cos_freq));


sIF_1 = exp(-1i*2*pi*(fc*Tau1+alpha*t*Tau1+epsi_t-epsi_td1));
sIF_2 = exp(-1i*2*pi*(fc*Tau2+alpha*t*Tau2+epsi_t-epsi_td2));

delTau = 2*(zn1-zn2)/c;

sdel = exp(-1i*2*pi*((fc+alpha*t)*delTau));

delsIF = sIF_1.*conj(sIF_2).*conj(sdel);

epsi_prime = -1*filloutliers(angle(delsIF),"linear")/delTau;
% epsi_prime = epsi_prime - mean(epsi_prime);

epsi_est = cumtrapz(t,epsi_prime)/2/pi;
s_epsi = exp(-1i*2*pi*epsi_est);

%% Upsampling and Calibration

del_t_new = 1/fBW; % Upsampling to RF bandwidth
t_new = del_t_new:del_t_new:Tp;
Nt_new = length(t_new); %length of sampling vector
fs_new = 1/del_t_new;

freq_t = linspace(-fs_new/2+1,fs_new/2,Nt_new);
timeDisp = exp(-1i*pi*((freq_t).^2)./alpha);

s3_new = resample(s3, Nt_new/N,1); % create high sampled data from low sampled on.
s_epsi_new = resample(s_epsi, Nt_new/N,1);
range_new = (0:Nt_new-1)*RR;

sb2 = s3_new.*conj(s_epsi_new); % Removing epsilon term
sb3 = ifft(fftshift(fft(sb2)).*(timeDisp)); % Range deskew and dispersive fileter
s_epsi_T = ifft(fftshift(fft(s_epsi_new)).*conj(timeDisp)); % Epsilon_RVP term
s3_clean = sb3.*s_epsi_T;

%% Plotting

sF1 = 20*log10(abs(ifft(s1)));
sF3 = 20*log10(abs(ifft(s3)));
sF3_clean = 20*log10(abs(ifft(s3_clean)));

figure(1)
clf
plot(range,sF3-max(sF1),'LineWidth',2,'LineStyle','--','Color','k'); hold on
plot(range_new,sF3_clean-max(sF3_clean),'LineWidth',2,'Color','k');
legend('Nonlinear','Corrected','Location','northwest')
xlabel('Range (m)');
ylabel('IF signal (dB)');
grid on
hold on;
axis tight
view([0 90])
xlim([zn-1.2 zn+1.2])
ylim([-50 0])
% set(gcf,'position',[50,50,800,800])
set(gca,'FontSize',20)
set(gca,'GridAlpha',1)
set(gca,'GridLineStyle','--')
set(gca,'fontname','ariel')
set(gcf,'color','w');

figure(2); clf
plot(t,(unwrap(epsi_t)),'LineWidth',1.2,'LineStyle','--','Color','k'); hold on
plot(t,(epsi_est),'LineWidth',1.2,'Color','k');
legend('actual phase error','estimated phase error','Location','northwest')
xlabel('time (sec)');
ylabel('Error (radians)');
grid on
hold on;
axis tight
view([0 90])
% set(gcf,'position',[950,50,300,300])
% set(gca,'FontSize',8)
set(gca,'GridAlpha',1)
set(gca,'GridLineStyle','--')
set(gca,'fontname','ariel')
set(gcf,'color','w');

