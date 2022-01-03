function spikepopulation_ECoG_FC(fnameIFR, dnameECoG)
% fnameIFR = 'C:\Users\ck1be\OneDrive - Ulster University\Documents\all data\Desktop\current_works\oxford works\Jose contrb\FINAL COUNTDOWN\IFR\iFr0.mat';
% dnameECoG = 'C:\Users\ck1be\OneDrive - Ulster University\Documents\all data\Desktop\current_works\oxford works\Jose contrb\FINAL COUNTDOWN\ECoG2\m0';
% spikepopulation_ECoG_FC(fnameIFR, dnameECoG)
%
% fnameIFR = 'C:\Users\ck1be\OneDrive - Ulster University\Documents\all data\Desktop\current_works\oxford works\Jose contrb\FINAL COUNTDOWN\IFR\ifr5ms\iFr0.mat';
% dnameECoG = 'C:\Users\ck1be\OneDrive - Ulster University\Documents\all data\Desktop\current_works\oxford works\Jose contrb\FINAL COUNTDOWN\ECoG2\m0';
% spikepopulation_ECoG_FC(fnameIFR, dnameECoG)

rng('default');Fs = 200; % Hz (ECoG sampling frequency)

flowcut = 0.05; % Hz
fhighcut = 40; % Hz
Nepochs = 1800; % define it

% Read IFR data
D = load(fnameIFR);
tmp = fieldnames(D);
if (length(tmp) > 1)
    error('Unexpected: only one variable should be saved in this mat file: %s', st.name);
end
X_IFRall = cell2mat(D.(tmp{1}));
%X_IFRall = mean(X_IFRall(1:3600*Fs,:),2);
  X_IFRall = sum(X_IFRall(1:1.8e5,:),2);
clear D;
Nsamples = length(X_IFRall);

% epoching IFR
Nwinsize = Nsamples/Nepochs;
X_IFR = reshape(X_IFRall,[Nwinsize Nepochs]);
X_IFR = bsxfun(@minus,X_IFR,mean(X_IFR)); % remove baseline

% filtering IFR
flowpass = 45; % FC was measured per second (window size) with an slide of 30 seconds
[B,A] = butter(5, flowpass/(Fs/2));
X_IFR = filtfilt(B, A, X_IFR);

% Read ECoG data and labels; epoch and filter data
Fs_ECoG = 3e4;
ns = Nwinsize*Nepochs*Fs_ECoG/Fs;
if (mod(Fs_ECoG,Fs)~=0)
    error('Fs of ECoG must be a multiple of IFR''s Fs and 200 Hz');
end
d = dir(fullfile(dnameECoG, '*.mat'));
Nchann_ECoG = length(d);
X_ECoG = zeros(Nwinsize, Nepochs, Nchann_ECoG);
[B,A] = butter(5, flowpass/(Fs_ECoG/2));
data = [];
ecog_label = cell(1,Nchann_ECoG);
for it = 1:Nchann_ECoG
    load(fullfile(d(it).folder, d(it).name)); %#ok<LOAD>
    % epoching
    data = reshape(data(1:ns),[], Nepochs);
    % filtering
    data = filtfilt(B, A, data);
    % downsample to Fs
    X_ECoG(:,:,it) = downsample(data, Fs_ECoG/Fs);
    clear data;
    ecog_label{it} = d(it).name(5:end-4);
end

% Tapering IFR and ECoG data
NFFT = Nwinsize;
winhann = hanning(NFFT);%Removes the bordering effect
Y_IFR = bsxfun(@times,X_IFR,winhann);
Y_ECoG = bsxfun(@times,X_ECoG,winhann);

% PSD analysis
freq = Fs/2*linspace(0,1,NFFT/2+1);
flagfreq = (freq > flowcut & freq < fhighcut);
freq = freq(flagfreq);
Nfreq = length(freq);
Y_IFR = fft(Y_IFR,NFFT);
Y_IFR = Y_IFR(flagfreq,:,:);
PSD_IFR = sqrt(squeeze(mean(abs(Y_IFR).^2,2)));
figure; semilogy(freq, PSD_IFR.^2);
annotation('textbox', [0.36 0.95 0.28 0.03], 'String', 'Power Spectral Density of all channels ', 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
xlabel('Frequency (Hz)');
ylabel('Signal power');
Y_ECoG = fft(Y_ECoG,NFFT);
Y_ECoG = Y_ECoG(flagfreq,:,:);
PSD_ECoG = sqrt(squeeze(mean(abs(Y_ECoG).^2,2)));
figure; semilogy(freq, PSD_ECoG.^2);
annotation('textbox', [0.36 0.95 0.28 0.03], 'String', 'Power Spectral Density of all channels ', 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
xlabel('Frequency (Hz)');
ylabel('Signal power');

% FC analysis
COH = bsxfun(@times, Y_IFR, conj(Y_ECoG));
den = bsxfun(@times, PSD_IFR, PSD_ECoG);
COH = squeeze(mean(COH,2))./den;                                            % Coherence
COH = abs(COH).^2;
figure; plot(freq, COH(:,:,1));

% maximum statistics
Nr = 1000; % Number of Monte-Carlo simulation or resampling.
COH_maxstat = zeros(1,Nr);
for iter = 1:Nr
    if mod(iter,100)==0, disp(iter); end
    ind = randperm(Nepochs);
    COHs = bsxfun(@times, Y_IFR, conj(Y_ECoG(:,ind)));
    den = bsxfun(@times, PSD_IFR, PSD_ECoG);
    COHs = squeeze(mean(COHs,2))./den;
    COHs = abs(COHs).^2;
    COH_maxstat(:,iter) = max(COHs(:));
end
th = prctile(COH_maxstat,95);

figure;
flag = all(COH < th);
if ~all(flag)
    plot(freq, COH(:,~flag));
end
hold on;
plot(freq, COH(:,flag), '-', 'Color', [0.5 0.5 0.5]);
plot(freq([1 end]), [th th], '--k', 'LineWidth', 2);
hold off;
title('Using P<0.05 coherence with ECoG channels');
legend(ecog_label);
xlabel('Frequency (Hz)');
ylabel('Coherence');
