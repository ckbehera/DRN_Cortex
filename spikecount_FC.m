function spikecount_FC(X, labels, Fs, Nepochs, flowpass, flowcut, fhighcut, fhighcut2, file_Type)
% This function conducts coherence analysis between IFR channels.
%
%
% % Examples:
% % 1.- For file_Type = 1 format
% fname = 'C:\Users\sb00722029\OneDrive - Ulster University\Documents\all data\Desktop\current_works\oxford works\Jose contrb\FINAL COUNTDOWN\IFR\ruairi original\Mouse number 0.csv';
% NUM = xlsread(fname);
% fname_labels = 'C:\Users\sb00722029\OneDrive - Ulster University\Documents\all data\Desktop\current_works\oxford works\Jose contrb\FINAL COUNTDOWN\IFR\labels0.csv';
% fid = fopen(fname_labels);
% labels = textscan(fid, '%s', 'delimiter', '\n');
% fclose(fid);
% labels = labels{1};
% Fs0 = 3e4;
% timebin = 10; % ms
% type = 'num_spike_sort'; % sort by the number of spikes (ascending)
% [ifr,labels] = spikecount(NUM, labels, Fs0, timebin, type);
% Fs = 1000/timebin;
% Nepochs = 600;
% flowpass = 45; % Hz
% flowcut = 0; % Hz
% fhighcut = 45; % Hz
% fhighcut2 = 5; % Hz
% spikecount_FC(ifr', labels, Fs, Nepochs, flowpass, flowcut, fhighcut, fhighcut2);
%
% % 2.- For file_Type = 2 format:
% fname = ['C:\Users\sb00722029\OneDrive - Ulster University\Documents\all data\Desktop\current_works\' ...
%     'oxford works\Jose contrb\FINAL COUNTDOWN\IFR\ruairi original\spiketimes.csv'];
% NUM = xlsread(fname);
% labels = [];
% file_Type = 2;
% [ifr,labels] = spikecount(NUM(:,2:3), labels, Fs0, timebin, type, file_Type);
% Fs = 1000/timebin;
% Nepochs = 600;
% flowpass = 45; % Hz
% flowcut = 0; % Hz
% fhighcut = 45; % Hz
% fhighcut2 = 5; % Hz
% spikecount_FC(ifr', labels, Fs, Nepochs, flowpass, flowcut, fhighcut, fhighcut2, file_Type);
%
% See also spikecount
rng('default');
Nchann = size(X,2);

if ~exist('file_Type','var') || isempty(file_Type)
    file_Type = 1;
end

% epoching IFR
Nwinsize = floor(size(X,1)/Nepochs);
X = reshape(X(1:Nwinsize*Nepochs,:),[Nwinsize Nepochs Nchann]);
X = bsxfun(@minus,X,mean(X)); % remove baseline

% colors for different types of spike neurons
if (file_Type == 1)
    cols = [1 0 0; ... % slow irregular
        0 1 0; ... % slow regular
        0 0 1; ... % fast irregular
        0 1 1];    % fast regular
    colnames = {'slow irregular'; 'slow regular'; 'fast irregular'; 'fast regular'};
elseif (file_Type == 2)
    cols = [         0    0.4470    0.7410
        0.8500    0.3250    0.0980
        0.9290    0.6940    0.1250
        0.4940    0.1840    0.5560
        0.4660    0.6740    0.1880
        0.3010    0.7450    0.9330
        0.6350    0.0780    0.1840];
    colnames = unique(labels);
    cols = cols(mod(0:length(colnames)-1,7)+1,:);
end

[~,iloc] = ismember(labels, colnames);
cols = cols(iloc,:);

% Tapering
NFFT = Nwinsize;
winhann = hanning(NFFT);%Removes the bordering effect
Y = bsxfun(@times,X,winhann);

% filtering IFR
[B,A] = butter(5, flowpass/(Fs/2)); % default: lowpass digital Butterworth filter
Y = filtfilt(B, A, Y);

% PSD analysis
freq = Fs/2*linspace(0,1,NFFT/2+1);
flagfreq = (freq >= flowcut & freq < fhighcut);
freq = freq(flagfreq);
Y = fft(Y,NFFT);
Y = Y(flagfreq,:,:);
PSD = sqrt(squeeze(mean(abs(Y).^2,2)));
figure; h = semilogy(freq, PSD.^2);
for it = 1:length(h)
    h(it).Color = cols(it,:);
end
annotation('textbox', [0.36 0.95 0.28 0.03], 'String', 'Power Spectral Density of all channels ', 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
xlabel('Frequency (Hz)');
ylabel('Signal power');
legend(colnames);

drawnow;
print('-dpng','-r300', 'PSD');
close all;

% apply high cut #2
flagfreq = (freq < fhighcut2);
freq = freq(flagfreq);
Nfreq = length(freq);
Y = Y(flagfreq,:,:);
PSD = sqrt(squeeze(mean(abs(Y).^2,2)));

% FC analysis
denPSD = bsxfun(@times, PSD, reshape(PSD, [Nfreq 1 Nchann]));
COH = bsxfun(@times, Y, reshape(conj(Y),[Nfreq Nepochs 1 Nchann]));
den = reshape(imag(COH), Nfreq, []);
den = reshape(abs(hilbert(den)), [Nfreq Nepochs Nchann Nchann]);
EIC = squeeze(mean(imag(COH),2))./squeeze(sqrt(mean(den.^2,2)));
EIC = reshape(EIC, Nfreq, []);
EIC = reshape(abs(hilbert(EIC)).^2, [Nfreq Nchann Nchann]);
COH = squeeze(mean(COH,2))./denPSD;                                            % Coherence
iCOH = abs(imag(COH));                                                         % Imaginary coherence
COH = abs(COH).^2;
EIC(isnan(EIC)) = 0;
iCOH(isnan(iCOH)) = 0;
COH(isnan(COH)) = 0;

% maximum statistics
Nr = 100; % Number of Monte-Carlo simulation or resampling.
COH_maxstat = zeros(1,Nr);
EIC_maxstat = zeros(1,Nr);
iCOH_maxstat = zeros(1,Nr);
for iter = 1:Nr
    if mod(iter,10)==0, disp(iter); end
    ind = randperm(Nepochs);
    COHs = bsxfun(@times, Y, reshape(conj(Y(:,ind,:)),[Nfreq Nepochs 1 Nchann]));
    den = reshape(imag(COHs), Nfreq, []);
    den = reshape(abs(hilbert(den)), [Nfreq Nepochs Nchann Nchann]);
    EICs = squeeze(mean(imag(COHs),2))./squeeze(sqrt(mean(den.^2,2)));
    EICs = reshape(EICs, Nfreq, []);
    EICs = reshape(abs(hilbert(EICs)).^2, [Nfreq Nchann Nchann]);
    COHs = squeeze(mean(COHs,2))./denPSD;
    iCOHs = abs(imag(COHs)); 
    COHs = abs(COHs).^2;
    EICs(isnan(EICs)) = 0;
    iCOHs(isnan(iCOHs)) = 0;
    COHs(isnan(COHs)) = 0;
    COH_maxstat(iter) = max(COHs(:));
    EIC_maxstat(iter) = max(EICs(:));
    iCOH_maxstat(iter) = max(iCOHs(:));
end
th = [prctile(COH_maxstat,95), prctile(EIC_maxstat,95), prctile(iCOH_maxstat,95)];

Mask = cell(1,3);
stat = {COH EIC iCOH};
stat_names = {'COH' 'EIC' 'iCOH'};
for it = 1:3
    Mask{it} = stat{it}(freq>0.5 & freq<15,:,:).*(stat{it}(freq>0.5 & freq<15,:,:) > th(it));
%     Mask{it} = stat{it}(freq>0.5 & freq<4,:,:).*(stat{it}(freq>0.5 & freq<4,:,:) > th(it));
end

% labels = [num2str((1:Nchann)') repmat('-',[Nchann 1]) char(labels)];
labels = strtrim(cellstr(labels));

labels = [num2str((1:Nchann)') repmat('-',[Nchann 1]) char(labels)];
labels = strtrim(cellstr(labels));

for it = 1:3
    fc = squeeze(max(Mask{it}));
    
    figure; imagesc(fc - diag(diag(fc))); axis image; colorbar
    set(gca, 'YTick', 1:Nchann);
    set(gca, 'YTicklabel', strrep(labels,'_','\_'));
    set(gca, 'XTick', 1:Nchann);
    set(gca, 'XTicklabel', strrep(labels,'_','\_'));
    set(gca, 'XTickLabelRotation', 90);
    title('Significant FC between 0.5 and 4.0 Hz', 'FontSize', 12)
    %title('Significant FC between 0.5 and 15 Hz', 'FontSize', 12)
    
    drawnow;
    print('-dpng','-r300', sprintf('FC_matrix_%s',stat_names{it}));
    close all;
end

for it = 1:3
    for ich = 1:Nchann
        figure;
        indchann = 1:Nchann;
        indchann(ich) = [];
        label = labels;
        channlab = label{ich};
        label(ich) = [];
        fc = stat{it}(:,:,ich);
        fc(:,ich) = [];
        colstmp = cols;
        colstmp(ich,:) = [];
        flag = all(fc < th(it));
        label = label(~flag);
        if ~all(flag)
            h = plot(freq, fc(:,~flag), 'LineWidth', 2);
            colstmp = colstmp(~flag,:);
            for k = 1:length(h)
                h(k).Color = colstmp(k,:);
            end
        end
        hold on;
        if any(flag)
            plot(freq, fc(:,flag), '-', 'Color', [0.5 0.5 0.5]);
        end
        plot(freq([1 end]), th(it)*[1 1], '--k', 'LineWidth', 2);
        hold off;
        % title(sprintf('Using P<0.05 coherence with Channel %s',channlab));
        title(strrep(channlab,'_','\_'));
        if ~all(flag)
            %legend(label);
            legend(strrep(strtrim(cellstr(label)),'_','\_'));
        end
        xlabel('Frequency (Hz)');
        ylabel('Coherence');
        
        %pause;
        
        drawnow;
        print('-dpng','-r300',sprintf('%s_Neuron_id_%s',stat_names{it},channlab));
        close all;
    end
end

fid = fopen('summary.txt','wt');
for it = 1:3
    fprintf(fid, '%s:\n\n', stat_names{it});

    fc = bsxfun(@minus, stat{it}, reshape(eye(Nchann),[1 Nchann Nchann]));
    maxfc = max(fc(:));
    minfc = min(fc(:));
    
    ind = find(fc == maxfc);
    [ii,jj,kk] = ind2sub(size(fc), ind(1));
    
    fprintf(fid, 'Maximum connectivity value: %.2f\n', maxfc);
    
    fprintf(fid, 'at frequency: %.2f\n', freq(ii));
    fprintf(fid, 'and between channels %d and %d.\n', jj, kk);
    
    ind1 = find(fc == minfc);
    [ii,jj,kk] = ind2sub(size(fc), ind1(1));
    fprintf(fid, 'Minimum connectivity value: %.2f\n', minfc);
    
    fprintf(fid, 'at frequency: %.2f\n', freq(ii));
    fprintf(fid, 'and between channels %d and %d.\n', jj, kk);
end
fclose(fid);