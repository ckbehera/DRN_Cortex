function spikecount_ECoG_FC(X_IFR, labels_IFR, Fs_IFR, X_ECoG, labels_ECoG, Fs_ECoG, Nepochs, flowpass, fhighcut, file_Type)
% This function conducts coherence analysis between IFR and ECoG.

rng('default');
Nchann_IFR = size(X_IFR,2);
Nchann_ECoG = size(X_ECoG,2);

if ~exist('file_Type','var') || isempty(file_Type)
    file_Type = 1;
end

if (Fs_IFR > Fs_ECoG)
    error('IFR sampling frequency must be lower than ECoG sampling frequency.');
end

% epoching IFR
Nwinsize = floor(size(X_IFR,1)/Nepochs);
X_IFR = reshape(X_IFR(1:Nwinsize*Nepochs,:),[Nwinsize Nepochs Nchann_IFR]);
X_IFR = bsxfun(@minus,X_IFR,mean(X_IFR)); % remove baseline

% putting ECoG in same frequencies as IFR data
[B,A] = butter(5, flowpass/(Fs_ECoG/2));
X_ECoG = filtfilt(B, A, X_ECoG);
X_ECoG = resample(X_ECoG, Fs_IFR, Fs_ECoG);

% epoching ECoG
X_ECoG = reshape(X_ECoG(1:Nwinsize*Nepochs,:), [Nwinsize Nepochs Nchann_ECoG]);
X_ECoG = bsxfun(@minus, X_ECoG, mean(X_ECoG)); % remove baseline

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
    colnames = unique(labels_IFR);
    cols = cols(mod(0:length(colnames)-1,7)+1,:);
end

[~,iloc] = ismember(labels_IFR, colnames);
cols = cols(iloc,:);

% Tapering IFR and ECoG data
NFFT = Nwinsize;
winhann = hanning(NFFT);%Removes the bordering effect
Y_IFR = bsxfun(@times,X_IFR,winhann);
Y_ECoG = bsxfun(@times,X_ECoG,winhann);

% PSD analysis
freq = Fs_IFR/2*linspace(0,1,NFFT/2+1);
Y_IFR = fft(Y_IFR,NFFT);
Y_ECoG = fft(Y_ECoG,NFFT);

% apply high cut
flagfreq = (freq < fhighcut);
freq = freq(flagfreq);
Nfreq = length(freq);
Y_IFR = Y_IFR(flagfreq,:,:);
PSD_IFR = sqrt(squeeze(mean(abs(Y_IFR).^2,2)));
Y_ECoG = Y_ECoG(flagfreq,:,:);
PSD_ECoG = sqrt(squeeze(mean(abs(Y_ECoG).^2,2)));

% FC analysis

% COH = zeros(Nfreq, Nepochs, Nchann_IFR, Nchann_ECoG);
% for it1 = 1:Nchann_IFR
%     for it2 = 1:Nchann_ECoG
%         COH(:,:,it1,it2) = Y_IFR(:,:,it1).*conj(Y_ECoG(:,:,it2));
%     end
% end

% COH = zeros(Nfreq, Nepochs, Nchann_IFR, Nchann_ECoG);
% for it = 1:Nchann_ECoG
%     COH(:,:,:,it) = Y_IFR.*repmat(conj(Y_ECoG(:,:,it)), [1 1 Nchann_IFR]);
% end

COH = bsxfun(@times, Y_IFR, reshape(conj(Y_ECoG),[Nfreq Nepochs 1 Nchann_ECoG]));
den = bsxfun(@times, PSD_IFR, reshape(PSD_ECoG, [Nfreq 1 Nchann_ECoG]));
COH = squeeze(mean(COH,2))./den;                                            % Coherence
COH = abs(COH).^2;
% figure; plot(freq, COH(:,:,1));

% maximum statistics
Nr = 100; % Number of Monte-Carlo simulation or resampling.
COH_maxstat = zeros(1,Nr);
for iter = 1:Nr
    if mod(iter,10)==0, disp(iter); end
    ind = randperm(Nepochs);
    COHs = bsxfun(@times, Y_IFR, reshape(conj(Y_ECoG(:,ind,:)),[Nfreq Nepochs 1 Nchann_ECoG]));
    den = bsxfun(@times, PSD_IFR, reshape(PSD_ECoG, [Nfreq 1 Nchann_ECoG]));
    COHs = squeeze(mean(COHs,2))./den;
    COHs = abs(COHs).^2;
    COH_maxstat(:,iter) = max(COHs(:));
end
th = prctile(COH_maxstat,95);
% th2 = prctile(COH_maxstat,99); %#ok<NASGU>

labels_IFR = [num2str((1:Nchann_IFR)') repmat('-',[Nchann_IFR 1]) char(labels_IFR)];
labels_IFR = strtrim(cellstr(labels_IFR));
    
% Mask = COH(freq>0.5 & freq<4,:,:).*(COH(freq>0.5 & freq<4,:,:) > th);
Mask = COH(freq>0.5 & freq<15,:,:).*(COH(freq>0.5 & freq<15,:,:) > th);

% fc = squeeze(mean(Mask));
fc = squeeze(max(Mask));

figure; imagesc(fc); axis image; colorbar
set(gca, 'YTick', 1:Nchann_IFR);
set(gca, 'YTicklabel', strrep(labels_IFR,'_','\_'));
set(gca, 'XTick', 1:Nchann_ECoG);
set(gca, 'XTicklabel', strrep(labels_ECoG,'_','\_'));
set(gca, 'XTickLabelRotation', 90);
title('Significant FC between 0.5 and 15 Hz', 'FontSize', 12)
% title('Significant FC between 0.5 and 4.0 Hz', 'FontSize', 12)

drawnow;
print('-dpng','-r300', 'XCoh_matrix');
close all;

figure;
for ich = 1:Nchann_ECoG
    channlab = labels_ECoG{ich};
    fc = COH(:,:,ich);
    flag = all(fc < th);
    label = labels_IFR(~flag);
    cla;
    if ~all(flag)
        h = plot(freq, fc(:,~flag), 'LineWidth', 2);
        colstmp = cols(~flag,:);
        for it = 1:length(h)
            h(it).Color = colstmp(it,:);
        end
    end
    hold on;
    if any(flag)
        plot(freq, fc(:,flag), '-', 'Color', [0.5 0.5 0.5]);
    end
    plot(freq([1 end]), [th th], '--k', 'LineWidth', 2);
    hold off;
    % title(sprintf('Using P<0.05 coherence with Channel %s',channlab));
    title(strrep(channlab,'_','\_'));
     legend(strrep(colnames,'_','\_'));
    if ~all(flag)
        % legend(label);
        legend(strrep(strtrim(cellstr(label)),'_','\_'));
    end
    xlabel('Frequency (Hz)');
    ylabel('Coherence');
    
    %pause;
    
    drawnow;
    print('-dpng','-r300',sprintf('coherence with Channel %s',channlab));
    savefig(sprintf('coherence with Channel %s',channlab));
    fnm=sprintf('coherence with Channel %s.fig',channlab);
    savefig(fnm);
    close all;
end

fid = fopen('summary.txt','wt');

fprintf(fid, '%s:\n\n', channlab);

fc = COH;

%     fc = bsxfun(@minus, COH, reshape(eye(Nchann_ECoG,Nchann_IFR),[1 Nchann_ECoG Nchann_IFR]));
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

fclose(fid);
