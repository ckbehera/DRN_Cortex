function [ifr,labels,tsec] = spikecount(NUM, labels, Fs0, timebin, type, file_Type, fplot, winsize)
% This function estimates the IFR of spike trains for user defined bin size.
%
% fplot = true then visualize ifr using winsize for the length of the
% window.
% Examples:
%
% % 1.- For file_Type = 1 format:
% fname = 'C:\Users\sb00722029\OneDrive - Ulster University\Documents\all data\Desktop\current_works\oxford works\Jose contrb\FINAL COUNTDOWN\IFR\ruairi original\Mouse number 0.csv';
% NUM = xlsread(fname);
% fname_labels = 'C:\Users\sb00722029\OneDrive - Ulster University\Documents\all data\Desktop\current_works\oxford works\Jose contrb\FINAL COUNTDOWN\IFR\labels0.csv';
% fid = fopen(fname_labels);
% labels = textscan(fid, '%s', 'delimiter', '\n');
% fclose(fid);
% labels = labels{1};
% Fs0 = 3e4;
% timebin = 10; % bin size in ms
% type = 'num_spike_sort'; % sort by the number of spikes (ascending)
% fplot = true;
% winsize = 0.2; % in seconds
% [ifr,labels,tsec] = spikecount(NUM, labels, Fs0, timebin, type, [], fplot, winsize); % by default file_Type = 1
%
% % 2.- For file_Type = 2 format:
% fname = ['C:\Users\sb00722029\OneDrive - Ulster University\Documents\all data\Desktop\current_works\' ...
%     'oxford works\Jose contrb\FINAL COUNTDOWN\IFR\ruairi original\spiketimes.csv'];
% NUM = xlsread(fname);
% labels = [];
% file_Type = 2;
% [ifr,labels] = spikecount(NUM(:,2:3), labels, Fs0, timebin, type, file_Type, fplot, winsize);
%

if ~exist('fplot','var')
    fplot = false;
    winsize = [];
end
if (fplot == false)
    winsize = [];
end
if ~exist('type','var')
    type = 'default';
end
if ~exist('file_Type','var') || isempty(file_Type)
    file_Type = 1; % This is two distinguish the 2 different formats of spiketrain data. 1: ticks in seconds, 2: ticks in samples.
end
Fs = 1000/timebin;
N = Fs0/Fs;
if (N - round(N) > 1e-6)
    error('Fs0 must be a multiple of 1000/timebin');
end
N = round(N);
if (file_Type == 1)
    tick = round(NUM(:,2)*Fs0);
elseif (file_Type == 2)
    tick = NUM(:,2);
end
label = unique(NUM(:,1));
if isempty(labels)
    labels = strtrim(cellstr(num2str(label)));
end
[~,iloc] = ismember(NUM(:,1),label);
ts = sparse(iloc,tick,1,length(label),max(tick));
if (nnz(ts>1) > 0)
    warning('There are repeated spikes: %d', nnz(ts>1));
end
ts = logical(ts);
ts = ts(:, 1:N*floor(size(ts,2)/N));
switch type
    case 'default'
        
    case 'num_spike_sort'
        % sorting channels by number of spikes
        numspikes = full(sum(ts,2));
        [~,iloc] = sort(numspikes);
        ts = ts(iloc,:);
        labels = labels(iloc);
end
% calculate IFR
[nchann, ns] = size(ts);
ifr = zeros(nchann, ns/N);
for it = 1:nchann
    disp(it);
    tmp = reshape(full(ts(it,:)), [N ns/N]);
    ifr(it,:) = sum(tmp);
end
ifr = 1000*ifr/timebin;
% plot if required
if fplot
    labeltmp = [num2str((1:nchann)') repmat('-',[nchann 1]) char(labels)];
    labeltmp = strtrim(cellstr(labeltmp));
    tick = 0;
    tick0_ifr = 0;
    
    figure;
    h1 = axes('position', [0.07 0.52 0.87 0.45]);
    h2 = axes('position', [0.07 0.04 0.91 0.45]);
    while (tick + winsize*Fs0 < size(ts,2))
        % subplot 211; % raster-plot
        axes(h1); cla;
        spy(ts(:,tick+1:tick+winsize*Fs0), 12); axis normal;
        tmp = tick + get(gca,'XTick');
        tmp = cellstr(num2str(roundn(tmp/Fs0,-2)'));
        set(gca, 'XTickLabel', tmp);
        set(gca, 'YTick', 1:nchann);
        set(gca, 'YTicklabel', labeltmp);
        ylabel('Spiking neurons ID');
        tick = tick + winsize*Fs0;
        
        % subplot 212 % IFR
        tick_ifr = tick0_ifr + (1:(winsize*Fs));
        tt = tick_ifr/Fs;
        axes(h2); cla;
        imagesc(tt, 1:size(ifr,1), ifr(:,tick_ifr)*timebin/1000);
        set(gca, 'YTick', 1:nchann);
        set(gca, 'YTicklabel', labeltmp);
        xlabel('Seconds');
        ylabel('Spiking neurons ID');
        tick0_ifr = tick_ifr(end);
        colorbar;
        
        pause
    end
end
tsec = (1:size(ifr,2))/Fs;
%ifr=mean(ifr);
%ifr = spikecount(NUM, labels, Fs0, timebin, type);
% Fsnew = 1000/timebin; t0 = timebin/2 + (0:size(ifr,2)-1)/Fsnew;
% figure; plot(t0, ifr(end,:))
% flag = (t0 >= 0) & (t0 <= 50);
% figure; plot(t0(flag), ifr(end,flag))
% figure; plot(t0(flag), ifr(end,flag)*timebin/1000)

% figure;
% tick0 = 0; ich = [24:30]; wins = 20;
% while (tick0 + wins*Fstmp < size(ifr,2))
%     Fstmp = 1000/timebin;
%     tick = tick0+(1:(wins*Fstmp));
%     tt = tick/Fs;
%     imagesc(tt, ich, ifr(ich,tick));
%     pause;
%     tick0 = tick(end);
% end


% flag = (tsec>2000 & tsec < 2001); figure; plot(tsec(flag), spike_count([26 28],flag))
% figure; plot(tsec, spike_count([26 28],:))

% flag = (tsec>2000 & tsec < 2001); figure; plot(tsec(flag), spike_count(:,flag))