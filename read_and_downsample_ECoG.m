wd = 'C:\Users\ck1be\OneDrive - Ulster University\Documents\all data\Desktop\current_works\oxford works\oxford project2\ruairi2nd\chandan';
basedir = pwd;
addpath(basedir);

STATUS_IFR_ECoG = true;
% STATUS_IFR_ECoG = false;

STATUS_ECoG = true;
%STATUS_ECoG = false;

STATUS_IFR = true;
% STATUS_IFR = false;

%% Read ECoG
if STATUS_ECoG || STATUS_IFR_ECoG
    fid = fopen(fullfile(wd,'analog_signals.csv'));
    varname_ECoG = textscan(fid, '%s', 1, 'delimiter', '\n');%names of every variable
    varname_ECoG = varname_ECoG{1};
    C_ECoG = textscan(fid, '%d%f%s%s%f', inf, 'delimiter', ',');% values of every column of every variable
    fclose(fid);
    
    % unique(C{3})
    % unique(C{4})
end

%% prepare ECoG data
if STATUS_ECoG || STATUS_IFR_ECoG
    tmp = unique(C_ECoG{3});
    N_mouse_ECoG = length(tmp);
    mouse_ECoG = [];
    for it = 1:N_mouse_ECoG
        mouse_ECoG(it).strcode = tmp{it}; %#ok<*SAGROW> session name
        flag = strcmp(C_ECoG{3}, mouse_ECoG(it).strcode);%
        mouse_ECoG(it).channel_id = unique(C_ECoG{4}(flag));% diff channel ids
    end
    
    for it = 1:N_mouse_ECoG
        disp([it N_mouse_ECoG]);
        mouse_ECoG(it).channel_samples = cell(1,length(mouse_ECoG(it).channel_id));%sample data
        mouse_ECoG(it).channel_num_samples = zeros(1,length(mouse_ECoG(it).channel_id));%no. of samples
        mouse_ECoG(it).channel_Fs = zeros(1,length(mouse_ECoG(it).channel_id));%sampling freq of each channel
        maxsample = 0;
        flag = strcmp(C_ECoG{3}, mouse_ECoG(it).strcode); % row entries for this mouse
        for k = 1:length(mouse_ECoG(it).channel_id)
            flag2 = strcmp(C_ECoG{4}, mouse_ECoG(it).channel_id{k}); % row entries for this channel
            mouse_ECoG(it).channel_samples{k} = C_ECoG{2}(flag & flag2);
            mouse_ECoG(it).channel_num_samples(k) = length(mouse_ECoG(it).channel_samples{k});
            tmp = C_ECoG{5}(flag & flag2);
            mouse_ECoG(it).channel_Fs(k) = 1/(tmp(2) - tmp(1));%sampling freq of ecog
        end
    end
    
    for it = 1:N_mouse_ECoG
        if any(diff(mouse_ECoG(it).channel_num_samples) > 0) || any(diff(mouse_ECoG(it).channel_Fs) > 0)
            error('These channels'' data was not acquired simultaneously.');
        end
        mouse_ECoG(it).recording_time = mouse_ECoG(it).channel_num_samples(1)/mouse_ECoG(it).channel_Fs(1)/60; %recording time in minutes
    end
end

%% Analyze ECoG data
if STATUS_ECoG
    Nepochs = 180;%
    flowcut = 0.05; % Hz
    fhighcut = 40; % Hz
    for it = 1:N_mouse_ECoG
        fprintf('MOUSE #%d of %d\n', it, N_mouse_ECoG);
        workdir = ['ECoG_' mouse_ECoG(it).strcode];
        mkdir(workdir);
        cd(workdir);
         
        Fs = mouse_ECoG(it).channel_Fs(1);
        X = cell2mat(mouse_ECoG(it).channel_samples);
        labels = mouse_ECoG(it).channel_id;
        ECoG_FC(X, Fs, labels, Nepochs, flowcut, fhighcut);
        
        cd(basedir);
    end
end
clear fid;
%% Read spiking data
if STATUS_IFR || STATUS_IFR_ECoG
    fid = fopen(fullfile(wd,'spiketimes.csv'));
    varname_spiketimes = textscan(fid, '%s', 1, 'delimiter', '\n');
    varname_spiketimes = varname_spiketimes{1};
    C_spiketimes = textscan(fid, '%d%d%d', inf, 'delimiter', ',');%data

    fclose(fid);
end

%% Read spiking cluster features
if STATUS_IFR || STATUS_IFR_ECoG
    fid = fopen(fullfile(wd,'manual_c.csv'));
%     fid = fopen(fullfile(wd,'clusters.csv'));
%     fid = fopen(fullfile(wd,'clusters - Copy.csv'));
    
    varname_clusters = textscan(fid, '%s', 1, 'delimiter', '\n');
    varname_clusters = varname_clusters{1};
    
%     C_clusters = textscan(fid, '%d%f%f%s', inf, 'delimiter', ',');
    C_clusters = textscan(fid, '%d%s%f%f', inf, 'delimiter', ',');
%     C_clusters = textscan(fid, '%d%s%s', inf, 'delimiter', ',');
    
    fclose(fid);
    
    % There is something suspicious because neuron_id 1753 appeared in the
    % clusters but how if this id is not in the spike data (spiketimes.csv)
    % setdiff(unique(C_clusters{1}), unique(C_spiketimes{2}))
end
clear fid;
%% Read neurons classes
if STATUS_IFR || STATUS_IFR_ECoG
    fid = fopen(fullfile(wd,'neurons.csv'));
    varname_neuron_class = textscan(fid, '%s', 1, 'delimiter', '\n');
    varname_neuron_class = varname_neuron_class{1};
   % C_neuron_class = textscan(fid, '%d%d%d', inf, 'delimiter', ',');
     C_neuron_class = textscan(fid, '%d%d%d%s', inf, 'delimiter', ',');
    fclose(fid);
end

%% intersect neuron class with cluster data
[flag,iloc] = ismember(C_neuron_class{2}, C_clusters{1});
iloc = iloc(flag);
all_neurons_session_name = C_neuron_class{4}(flag);
all_neurons_id = C_clusters{1}(iloc);
all_neurons_cluster = C_clusters{2}(iloc);
all_neurons_cv2_isi = C_clusters{3}(iloc);
all_neurons_mean_firing_rate = C_clusters{4}(iloc);

%% prepare IFR data
if STATUS_IFR || STATUS_IFR_ECoG
    tmp = unique(all_neurons_session_name);
    N_mouse_IFR = length(tmp);
    mouse_IFR = [];
    for it = 1:N_mouse_IFR
        mouse_IFR(it).strcode = tmp{it}; %#ok<*SAGROW>
    end
    for it = 1:N_mouse_IFR
        flag = strcmp(all_neurons_session_name, mouse_IFR(it).strcode);
        mouse_IFR(it).neuron_id = all_neurons_id(flag);
        mouse_IFR(it).neuron_cluster = all_neurons_cluster(flag);
        mouse_IFR(it).neuron_cv2_isi = all_neurons_cv2_isi(flag);
        mouse_IFR(it).neuron_mean_firing_rate = all_neurons_mean_firing_rate(flag);
        tmp = unique(mouse_IFR(it).neuron_id);
        if (length(tmp) < length(mouse_IFR(it).neuron_id))
            error('Some weird duplicated neuron_id in neurons.csv dataset');
        end
        
        mouse_IFR(it).neuron_samples = cell(1,length(mouse_IFR(it).neuron_id));
        mouse_IFR(it).neuron_num_samples = zeros(1,length(mouse_IFR(it).neuron_id));
        maxsample = 0;
        for k = 1:length(mouse_IFR(it).neuron_id)
            flag = (C_spiketimes{2} == mouse_IFR(it).neuron_id(k));
            %flag=~flag;
           % mouse_IFR(it).neuron_samples{k} = C_spiketimes{2}(flag);
            mouse_IFR(it).neuron_samples{k} = C_spiketimes{3}(flag);
            mouse_IFR(it).neuron_num_samples(k) = length(mouse_IFR(it).neuron_samples{k});
            if (mouse_IFR(it).neuron_num_samples(k) == 0)
                warning('Mouse %s: neuron id #%d has zero spikes.', mouse_IFR(it).strcode, mouse_IFR(it).neuron_id(k));
            else
                maxsample = max(maxsample, max(mouse_IFR(it).neuron_samples{k}));
                 
            end
        end
        mouse_IFR(it).recording_time = maxsample/3e4/60; % time in minutes
        
    end
end

%% Analyze spiking data
if STATUS_IFR 
    Fs0 = 3e4;
    timebin = 5; % bin size in ms
    
    type = 'default'; % change to 'num_spike_sort' when you have the fast/slow, regular/irregular classification
    % type = 'num_spike_sort'; % sort by the number of spikes (ascending)
    
    file_Type = 2;
    
%     fplot = true;
%     winsize = 0.2; % in seconds 
    
    for it = 1:N_mouse_IFR
        %%close all, clc
        fprintf('MOUSE #%d of %d\n', it, N_mouse_IFR);
        workdir = ['IFR_' mouse_IFR(it).strcode];
        mkdir(workdir);
        cd(workdir);
        
        NUM = zeros(sum(mouse_IFR(it).neuron_num_samples), 2);
        cnt = 0;
        for k = 1:length(mouse_IFR(it).neuron_id)
            ind = cnt + (1:mouse_IFR(it).neuron_num_samples(k));
            NUM(ind,1) = mouse_IFR(it).neuron_id(k);
            NUM(ind,2) = mouse_IFR(it).neuron_samples{k};
            cnt = cnt + mouse_IFR(it).neuron_num_samples(k);
        end
       [ifr,labels,tsec] = spikecount(NUM, mouse_IFR(it).neuron_cluster, Fs0, timebin, type, file_Type);
      %  [ifr,labels,tsec] = spikecount(NUM, labels, Fs0, timebin, type, file_Type, fplot, winsize);
%         [ifr,labels,tsec] = spikecount(NUM, mouse_IFR(it).neuron_cluster, Fs0, timebin, type, file_Type, fplot, winsize);
        Fs = 1000/timebin;
        Nepochs = 600;
        flowpass = 45; % Hz
        flowcut = 0; % Hz
        fhighcut = 45; % Hz
        fhighcut2 = 15; % Hz
        spikecount_FC(ifr', labels, Fs, Nepochs, flowpass, flowcut, fhighcut, fhighcut2, file_Type);
        
        cd(basedir);
    end
end

%% Analyze cross-coherence between IFR and ECoG
if STATUS_IFR_ECoG
    Fs0 = 3e4;
    timebin = 5 ; % bin size in ms
    type = 'default'; % change to 'num_spike_sort' when you have the fast/slow, regular/irregular classification
    % type = 'num_spike_sort'; % sort by the number of spikes (ascending)
    file_Type = 2;
    [exist_mouse,iloc] = ismember({mouse_IFR.strcode}, {mouse_ECoG.strcode});
    
    for it = 1:N_mouse_IFR
        if ~exist_mouse(it), continue; end
        
        fprintf('MOUSE #%d of %d\n', it, N_mouse_IFR);
        workdir = ['IFR_ECoG_' mouse_IFR(it).strcode];
        mkdir(workdir);
        cd(workdir);
        
        NUM = zeros(sum(mouse_IFR(it).neuron_num_samples), 2);
        cnt = 0;
        for k = 1:length(mouse_IFR(it).neuron_id)
            ind = cnt + (1:mouse_IFR(it).neuron_num_samples(k));
            NUM(ind,1) = mouse_IFR(it).neuron_id(k);
            NUM(ind,2) = mouse_IFR(it).neuron_samples{k};
            cnt = cnt + mouse_IFR(it).neuron_num_samples(k);
        end
        [X_IFR,labels_IFR,tsec] = spikecount(NUM, mouse_IFR(it).neuron_cluster, Fs0, timebin, type, file_Type);
        Fs_IFR = 1000/timebin;
        
        Fs_ECoG = mouse_ECoG(iloc(it)).channel_Fs(1);
        X_ECoG = cell2mat(mouse_ECoG(iloc(it)).channel_samples);
        labels_ECoG = mouse_ECoG(iloc(it)).channel_id;
        
        Nepochs = 180;
        flowpass = 25;
        flowcut = 0.05; % Hz
        fhighcut = 15; % Hz
        
        %spikecount_ECoG_FC_temp(X_IFR', labels_IFR, Fs_IFR, X_ECoG, labels_ECoG, Fs_ECoG, Nepochs, flowpass, fhighcut, file_Type);
          
        spikecount_ECoG_FC(X_IFR', labels_IFR, Fs_IFR, X_ECoG, labels_ECoG, Fs_ECoG, Nepochs, flowpass, fhighcut, file_Type);
        
        cd(basedir);
    end
end