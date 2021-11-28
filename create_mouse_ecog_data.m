wd = 'C:\Users\chandan';

%% Read ECoG
fid = fopen(fullfile(wd,'analog_signals.csv'));
varname_ECoG = textscan(fid, '%s', 1, 'delimiter', '\n');
varname_ECoG = varname_ECoG{1};
C_ECoG = textscan(fid, '%d%f%s%s%f', inf, 'delimiter', ',');
fclose(fid);

%% prepare ECoG data
tmp = unique(C_ECoG{3});
N_mouse_ECoG = length(tmp);
mouse_ECoG = [];
for it = 1:N_mouse_ECoG
    mouse_ECoG(it).strcode = tmp{it}; %#ok<*SAGROW>
    flag = strcmp(C_ECoG{3}, mouse_ECoG(it).strcode);
    mouse_ECoG(it).channel_id = unique(C_ECoG{4}(flag));
end

for it = 1:N_mouse_ECoG
    disp([it N_mouse_ECoG]);
    mouse_ECoG(it).channel_samples = cell(1,length(mouse_ECoG(it).channel_id));
    mouse_ECoG(it).channel_num_samples = zeros(1,length(mouse_ECoG(it).channel_id));
    mouse_ECoG(it).channel_Fs = zeros(1,length(mouse_ECoG(it).channel_id));
    maxsample = 0;
    flag = strcmp(C_ECoG{3}, mouse_ECoG(it).strcode); % row entries for this mouse
    for k = 1:length(mouse_ECoG(it).channel_id)
        flag2 = strcmp(C_ECoG{4}, mouse_ECoG(it).channel_id{k}); % row entries for this channel
        mouse_ECoG(it).channel_samples{k} = C_ECoG{2}(flag & flag2);
        mouse_ECoG(it).channel_num_samples(k) = length(mouse_ECoG(it).channel_samples{k});
        tmp = C_ECoG{5}(flag & flag2);
        mouse_ECoG(it).channel_Fs(k) = 1/(tmp(2) - tmp(1));
    end
end

for it = 1:N_mouse_ECoG
    if any(diff(mouse_ECoG(it).channel_num_samples) > 0) || any(diff(mouse_ECoG(it).channel_Fs) > 0)
        error('These channels'' data was not acquired simultaneously.');
    end
    mouse_ECoG(it).recording_time = mouse_ECoG(it).channel_num_samples(1)/mouse_ECoG(it).channel_Fs(1)/60; %recording time in minute
end
