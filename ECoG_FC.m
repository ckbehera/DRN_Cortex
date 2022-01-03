function ECoG_FC(X, Fs, ecog_label, Nepochs, flowcut, fhighcut)
% labels = [];
% Nepochs = 180;
% flowcut = 0.05; % Hz
% fhighcut = 40; % Hz
% ECoG_FC(X, Fs, labels, Nepochs, flowcut, fhighcut)

Nwinsize = floor(size(X,1)/Nepochs);
Nchann = size(X,2);

ind = 1:2*Fs;
figure;
plot(ind/Fs, X(1e4+ind,:)); axis tight;
xlabel('Seconds');
ylabel('Signal');
title('Raw signal activity');

drawnow;
print('-dpng','-r300', 'signal_raw_activity');
close all;

X = reshape(X(1:Nwinsize*Nepochs,:), [Nwinsize Nepochs Nchann]);
X = bsxfun(@minus, X, mean(X)); % remove baseline

% Tapering ECoG data
NFFT = Nwinsize;
winhann = hanning(NFFT); % Removes the bordering effect
Y = bsxfun(@times,X,winhann);

% PSD analysis
freq = Fs/2*linspace(0,1,NFFT/2+1);
flagfreq = (freq > flowcut & freq < fhighcut);
freq = freq(flagfreq);
Nfreq = length(freq);
Y = fft(Y,NFFT);
Y = Y(flagfreq,:,:);
PSD = sqrt(squeeze(mean(abs(Y).^2,2)));
figure; semilogy(freq, PSD.^2);
annotation('textbox', [0.36 0.95 0.28 0.03], 'String', 'Power Spectral Density of all channels ', 'FontSize', 12, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
xlabel('Frequency (Hz)');
ylabel('Signal power');
legend(strrep(ecog_label,'_','\_'));

drawnow;
print('-dpng','-r300', 'PSD_signal');
savefig('PSD_signal.fig');
close all;

% FC analysis
denPSD = bsxfun(@times, PSD, reshape(PSD, [Nfreq 1 Nchann]));
COH = bsxfun(@times, Y, reshape(conj(Y),[Nfreq Nepochs 1 Nchann]));
den = reshape(imag(COH), Nfreq, []);
den = reshape(abs(hilbert(den)), [Nfreq Nepochs Nchann Nchann]);
EIC = squeeze(mean(imag(COH),2))./squeeze(sqrt(mean(den.^2,2)));
EIC = reshape(EIC, Nfreq, []);
EIC = reshape(abs(hilbert(EIC)).^2, [Nfreq Nchann Nchann]);
COH = squeeze(mean(COH,2))./denPSD;                                            % Coherence
iCOH = abs(imag(COH));                                                      % Imaginary coherence
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
end

% Mask = COH(freq>0.5 & freq<4,:,:).*(COH(freq>0.5 & freq<4,:,:) > th);

for it = 1:3
    fc = squeeze(max(Mask{it}));
    
    figure; imagesc(fc - diag(diag(fc))); axis image; colorbar
    set(gca, 'YTick', 1:Nchann);
    set(gca, 'YTicklabel', strrep(ecog_label,'_','\_'));
    set(gca, 'XTick', 1:Nchann);
    set(gca, 'XTicklabel', strrep(ecog_label,'_','\_'));
    set(gca, 'XTickLabelRotation', 90);
    title('Significant FC between 0.5 and 4.0 Hz', 'FontSize', 12)
    
    drawnow;
    print('-dpng','-r300', sprintf('FC_matrix_%s',stat_names{it}));
%     fnm=sprintf('FC_matrix_%s.fig',stat_names{it});
%     savefig(fnm);
    close all;
end

for it = 1:3
    for ich = 1:Nchann
        figure;
        indchann = 1:Nchann;
        indchann(ich) = [];
        label = ecog_label;
        channlab = [num2str(ich) ' ' label{ich}];
        label(ich) = [];
        fc = stat{it}(:,:,ich);
        fc(:,ich) = [];
        flag = all(fc < th(it));
        label = [num2str(indchann(~flag)') repmat(' ', [nnz(~flag) 1]) char(label(~flag))];
        if ~all(flag)
            plot(freq, fc(:,~flag));
        end
        hold on;
        if any(flag)
            plot(freq, fc(:,flag), '-', 'Color', [0.5 0.5 0.5]);
        end
        plot(freq([1 end]), th(it)*[1 1], '--k', 'LineWidth', 2);
        hold off;
        title(sprintf('Using P<0.05 coherence with Channel %s',channlab));
        title(strrep(channlab,'_','\_'));
        legend(strrep(strtrim(cellstr(label)),'_','\_'));
        xlabel('Frequency (Hz)');
        ylabel('Coherence');
        %pause;
        
        drawnow;
        print('-dpng','-r300',sprintf('%s_Channel_%s',stat_names{it},channlab));
        fnm1=sprintf('%s_Channel_%s.fig',stat_names{it},channlab);
        savefig(fnm1);
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