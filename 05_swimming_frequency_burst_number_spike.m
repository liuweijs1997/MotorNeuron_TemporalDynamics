% analyze_swimming_burstFrequency_spikeCount
%
% Detect bursts from an extracellular ABF trace using threshold crossing.
%
% Definitions:
%   burstStartTimes  = burst onset times from threshold crossing
%   burstEndTimes    = burst offset times from threshold crossing
%   cycleDurations   = time from one burst onset to the next burst onset
%   burstFrequencies = 1 / cycleDurations
%   spikeCount       = number of spikes within each burst
%
% Main output:
%   x-axis: burst frequency (Hz)
%   y-axis: spike count per burst
%   one data point per burst
%
% Filters:
%   keep only bursts with:
%   2 < burst frequency < 10.5 Hz
%   burst duration / cycle > 0.1
%
% Output files:
%   - bursts.fig
%   - recording.fig
%   - spikeCount_vs_burstFreq.fig
%   - BurstFreqSpikeCountOverTime.fig
%   - spikeCount_vs_burstFreq.csv / .txt / .mat
%   - analysis.txt
%   - analysis.mat
%   - analysis_parameters.mat

clearvars;
close all;
clc;

%% choose file
[fileName, pathName, filterIndex] = uigetfile('*.abf');

if isequal(fileName, 0)
    error('No file selected.');
end

filePath = fullfile(pathName, fileName);

%% load abf
[data, si, h] = abfload(filePath);

%% time vector
samplingInterval = si / 1000000;               % seconds
samplingFrequency = 1 / samplingInterval;      % Hz
lengthOfData = size(data, 1);
time = (0:lengthOfData-1)' * samplingInterval; % column vector
maxTime = max(time);

%% preview channels
channelsToPlot = 1:min(size(data,2), 2); % preview first 2 channels if available

figurePreview = figure;
channelNames = h.recChNames;

for i = 1:length(channelsToPlot)
    ax(i) = subplot(length(channelsToPlot), 1, i);
    plot(time, data(:, channelsToPlot(i)));
    title(strcat(channelNames{channelsToPlot(i)}, ' channel# ', num2str(channelsToPlot(i))));
    ylabel('mV');
    xlabel('time (s)');
end

linkaxes(ax, 'x');

%% ask which channel to analyze
prompt = {'Which channel should be analyzed?'};
dlg_title = 'Enter channel number';
num_lines = 1;
def = {num2str(channelsToPlot(1))};
options.WindowStyle = 'normal';

answer = inputdlg(prompt, dlg_title, num_lines, def, options);

if isempty(answer)
    error('No channel selected.');
end

traceNumber = str2double(answer{1});

if isnan(traceNumber) || traceNumber < 1 || traceNumber > size(data,2)
    error('Invalid channel number.');
end

close(figurePreview);

%% choose analyzed trace
data = data(:, traceNumber);
data = data(:);

%% create output folder
baseFileName = fileName(1:end-4); % remove .abf
outputFolder = fullfile(pathName, baseFileName);

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% detect spikes
% NOTE:
% detect_spikes was originally designed for intracellular recordings.
% For extracellular data, always visually inspect whether spike detection is reasonable.
clear global phsc_data;
try
    spikeTimes = detect_spikes(data, si);
    spikeTimes = spikeTimes(:);
catch ME
    warning('Spike detection failed: %s. Continuing without spike analysis.', ME.message);
    spikeTimes = [];
end

%% plot raw + filtered trace
figureBursts = figure;
h1 = plot(time, data);
set(h1, 'DisplayName', 'recording');
xlabel('time (seconds)');
ylabel('signal (mV)');
hold on;

%% preprocess extracellular signal
dcWindow = round(0.025 / samplingInterval);
smoothWindow = round(0.025 / samplingInterval);

dcWindow = max(dcWindow, 1);
smoothWindow = max(smoothWindow, 1);

% sliding DC removal
filteredData = data - slidingavg(data, dcWindow);

% rectify
filteredData = abs(filteredData);

% smooth
filteredData = slidingavg(filteredData, smoothWindow);

% force column vector
filteredData = filteredData(:);

h2 = plot(time, filteredData, 'g');
set(h2, 'DisplayName', 'filtered data');

%% load or ask burst threshold
thresholdFile = fullfile(outputFolder, [baseFileName 'burstThreshold.mat']);

if exist(thresholdFile, 'file') == 2
    savedData = load(thresholdFile);
    burstThreshold = savedData.burstThreshold;
    clear savedData;
else
    prompt = {'Burst threshold (on filteredData):'};
    dlg_title = 'Enter burst detection threshold';
    num_lines = 1;
    def = {num2str(0.3 * max(filteredData), '%.4f')};
    options.WindowStyle = 'normal';

    answer = inputdlg(prompt, dlg_title, num_lines, def, options);

    if isempty(answer)
        error('No burst threshold provided.');
    end

    burstThreshold = str2double(answer{1});

    if isnan(burstThreshold) || burstThreshold <= 0
        error('Invalid burst threshold.');
    end

    save(thresholdFile, 'burstThreshold');
end

%% load or ask minimum burst duration
minDurFile = fullfile(outputFolder, [baseFileName 'minBurstDuration.mat']);

if exist(minDurFile, 'file') == 2
    savedData = load(minDurFile);
    minBurstDuration = savedData.minBurstDuration;
    clear savedData;
else
    prompt = {'Minimum burst duration (seconds):'};
    dlg_title = 'Enter minimum burst duration';
    num_lines = 1;
    def = {'0.02'};
    options.WindowStyle = 'normal';

    answer = inputdlg(prompt, dlg_title, num_lines, def, options);

    if isempty(answer)
        error('No minimum burst duration provided.');
    end

    minBurstDuration = str2double(answer{1});

    if isnan(minBurstDuration) || minBurstDuration <= 0
        error('Invalid minimum burst duration.');
    end

    save(minDurFile, 'minBurstDuration');
end

%% detect bursts
isBurst = filteredData > burstThreshold;

dBurst = diff([0; isBurst; 0]);
burstStartIdx = find(dBurst == 1);
burstEndIdx   = find(dBurst == -1) - 1;

if isempty(burstStartIdx) || isempty(burstEndIdx)
    error('No bursts detected. Try lowering the threshold.');
end

burstStartTimes = time(burstStartIdx);
burstEndTimes   = time(burstEndIdx);
burstDurations  = burstEndTimes - burstStartTimes;

%% remove short bursts
validBursts = burstDurations >= minBurstDuration;

burstStartIdx   = burstStartIdx(validBursts);
burstEndIdx     = burstEndIdx(validBursts);
burstStartTimes = burstStartTimes(validBursts);
burstEndTimes   = burstEndTimes(validBursts);
burstDurations  = burstDurations(validBursts);

if length(burstStartTimes) < 2
    error('Fewer than 2 valid bursts detected. Cannot calculate cycle durations/frequencies.');
end

%% calculate cycle, frequency, duty cycle
cycleDurations = burstStartTimes(2:end) - burstStartTimes(1:end-1);
burstFrequencies = 1 ./ cycleDurations;

% for burst i, pair its own duration with the cycle from burst i to burst i+1
burstDurationsForPlot = burstDurations(1:end-1);
burstDutyCycle = burstDurationsForPlot ./ cycleDurations;

cycleTimes = burstStartTimes(1:end-1);
burstStartForPlot = burstStartTimes(1:end-1);
burstEndForPlot   = burstEndTimes(1:end-1);
burstIndexForPlot = (1:length(burstFrequencies))';

%% count spikes in each burst
spikeCount_all = nan(length(burstFrequencies), 1);

if ~isempty(spikeTimes)
    for i = 1:length(burstFrequencies)
        spikesInBurst = spikeTimes(spikeTimes >= burstStartForPlot(i) & spikeTimes <= burstEndForPlot(i));
        spikeCount_all(i) = length(spikesInBurst);
    end
else
    warning('No spikes detected, so spike counts will be NaN.');
end

%% save raw copies before filtering
cycleDurations_raw      = cycleDurations;
burstFrequencies_raw    = burstFrequencies;
burstDurations_rawPlot  = burstDurationsForPlot;
burstDutyCycle_raw      = burstDutyCycle;
cycleTimes_raw          = cycleTimes;
burstStart_rawPlot      = burstStartForPlot;
burstEnd_rawPlot        = burstEndForPlot;
burstIndex_rawPlot      = burstIndexForPlot;
spikeCount_raw          = spikeCount_all;

%% raw table
T_spikeCount_raw = table( ...
    burstIndex_rawPlot, ...
    burstStart_rawPlot, ...
    burstEnd_rawPlot, ...
    cycleDurations_raw, ...
    burstFrequencies_raw, ...
    burstDurations_rawPlot, ...
    burstDutyCycle_raw, ...
    spikeCount_raw, ...
    'VariableNames', { ...
    'BurstIndex', ...
    'BurstStartTime_s', ...
    'BurstEndTime_s', ...
    'CycleDuration_s', ...
    'BurstFrequency_Hz', ...
    'BurstDuration_s', ...
    'BurstDuration_over_Cycle', ...
    'SpikeCount'});

%% apply filters
% keep only:
% 2 < burst frequency < 10.5 Hz
% burst duration / cycle > 0.1
validIdx = burstFrequencies > 2 & burstFrequencies < 10.5 & burstDutyCycle > 0.1;

cycleDurations      = cycleDurations(validIdx);
burstFrequencies    = burstFrequencies(validIdx);
burstDurationsForPlot = burstDurationsForPlot(validIdx);
burstDutyCycle      = burstDutyCycle(validIdx);
cycleTimes          = cycleTimes(validIdx);
burstStartForPlot   = burstStartForPlot(validIdx);
burstEndForPlot     = burstEndForPlot(validIdx);
burstIndexForPlot   = burstIndexForPlot(validIdx);
spikeCount_all      = spikeCount_all(validIdx);

if isempty(burstFrequencies)
    error('No bursts remained after applying the filters.');
end

%% filtered table
T_spikeCount = table( ...
    burstIndexForPlot, ...
    burstStartForPlot, ...
    burstEndForPlot, ...
    cycleDurations, ...
    burstFrequencies, ...
    burstDurationsForPlot, ...
    burstDutyCycle, ...
    spikeCount_all, ...
    'VariableNames', { ...
    'BurstIndex', ...
    'BurstStartTime_s', ...
    'BurstEndTime_s', ...
    'CycleDuration_s', ...
    'BurstFrequency_Hz', ...
    'BurstDuration_s', ...
    'BurstDuration_over_Cycle', ...
    'SpikeCount'});

%% plot burst onset/offset on filtered trace
figure(figureBursts);
plot(burstStartTimes, filteredData(burstStartIdx), 'bo', ...
    'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', 'burst start');
plot(burstEndTimes, filteredData(burstEndIdx), 'mo', ...
    'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', 'burst end');
yline(burstThreshold, '--b', 'burst threshold', 'LineWidth', 1.5);

%% raw recording figure
figureRecording = figure;
h3 = plot(time, data);
set(h3, 'DisplayName', 'recording');
axis([0 maxTime -Inf Inf]);
xlabel('time (seconds)');
ylabel('signal (mV)');

%% mean spike count at 4-5 Hz
indexFourToFiveHz = find(4 <= burstFrequencies & burstFrequencies <= 5);

if isempty(indexFourToFiveHz)
    meanSpikeCountFourToFiveHz = NaN;
else
    meanSpikeCountFourToFiveHz = mean(spikeCount_all(indexFourToFiveHz), 'omitnan');
end

%% plot RAW spike count vs frequency
figureSpikeCountRaw = figure;
scatter(burstFrequencies_raw, spikeCount_raw, ...
    'MarkerEdgeColor', [0.5 0.5 0.5], 'LineWidth', 1.0);
xlabel('burst frequency (Hz)');
ylabel('spike count per burst');
title('Raw burst-level spike count');

axRaw = gca;
set(axRaw, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.015 0.03], ...
    'LineWidth', 2.5, 'GridLineStyle', 'none', 'FontSize', 14, ...
    'FontName', 'Arial');

%% plot FILTERED spike count vs frequency
figureSpikeCount = figure;
scatter(burstFrequencies, spikeCount_all, ...
    'MarkerEdgeColor', [0 0 0], 'LineWidth', 1.5);
xlabel('burst frequency (Hz, 2<f<10.5)');
ylabel('spike count per burst');
title('Burst frequency vs spike count');

axSC = gca;
set(axSC, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.015 0.03], ...
    'LineWidth', 2.5, 'GridLineStyle', 'none', 'FontSize', 14, ...
    'FontName', 'Arial');

%% plot burst frequency and spike count over time
figureBurstFreqSpikeCountOverTime = figure;

yyaxis left
scatter(cycleTimes, burstFrequencies, 20, 'k', 'filled');
ylabel('burst frequency (Hz)', 'FontSize', 14, 'FontName', 'Arial');

yyaxis right
scatter(cycleTimes, spikeCount_all, 20, ...
    'MarkerEdgeColor', [0.5 0.5 0.5]);
ylabel('spike count per burst', 'FontSize', 14, 'FontName', 'Arial', ...
    'Color', [0.5 0.5 0.5]);

xlabel('time (s)', 'FontSize', 14, 'FontName', 'Arial');

axTime = gca;
set(axTime, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.015 0.03], ...
    'LineWidth', 2.5, 'GridLineStyle', 'none', 'FontSize', 14, ...
    'FontName', 'Arial');

%% save txt data
savefileTxt = fullfile(outputFolder, [baseFileName 'analysis.txt']);

saveMatrix = [ ...
    burstStartForPlot, ...
    burstEndForPlot, ...
    cycleDurations, ...
    burstFrequencies, ...
    burstDurationsForPlot, ...
    burstDutyCycle, ...
    spikeCount_all];

myfile = fopen(savefileTxt, 'wt');

fprintf(myfile, 'burst threshold:\t%4.6f\r\n', burstThreshold);
fprintf(myfile, 'minimum burst duration (s):\t%4.6f\r\n', minBurstDuration);
fprintf(myfile, 'filters:\t2 < burst frequency < 10.5 Hz; burst duration/cycle > 0.0\r\n');
fprintf(myfile, 'mean spike count at 4-5 Hz:\t%4.6f\r\n', meanSpikeCountFourToFiveHz);
fprintf(myfile, ['burst start time (s)\t burst end time (s)\t' ...
    'cycle duration (s)\t burst frequency (Hz)\t burst duration (s)\t' ...
    'burst duration/cycle\t spike count\r\n']);

fclose(myfile);

save(savefileTxt, 'saveMatrix', '-ascii', '-tabs', '-append');

%% save csv/txt/mat for RAW table
writetable(T_spikeCount_raw, fullfile(outputFolder, [baseFileName 'spikeCount_vs_burstFreq_raw.csv']));
writetable(T_spikeCount_raw, fullfile(outputFolder, [baseFileName 'spikeCount_vs_burstFreq_raw.txt']), ...
    'Delimiter', '\t');
save(fullfile(outputFolder, [baseFileName 'spikeCount_vs_burstFreq_raw.mat']), 'T_spikeCount_raw');

%% save csv/txt/mat for FILTERED table
writetable(T_spikeCount, fullfile(outputFolder, [baseFileName 'spikeCount_vs_burstFreq_filtered.csv']));
writetable(T_spikeCount, fullfile(outputFolder, [baseFileName 'spikeCount_vs_burstFreq_filtered.txt']), ...
    'Delimiter', '\t');
save(fullfile(outputFolder, [baseFileName 'spikeCount_vs_burstFreq_filtered.mat']), 'T_spikeCount');

%% save matlab variables
analysisParametersFile = fullfile(outputFolder, 'analysis_parameters.mat');
save(analysisParametersFile, ...
    'burstStartTimes', 'burstEndTimes', 'burstDurations', ...
    'cycleDurations', 'burstFrequencies', 'burstDurationsForPlot', ...
    'burstDutyCycle', 'cycleTimes', ...
    'burstThreshold', 'minBurstDuration', ...
    'spikeTimes', 'spikeCount_all', ...
    'T_spikeCount', 'T_spikeCount_raw');

%% save workspace
savefileMat = fullfile(outputFolder, [baseFileName 'analysis.mat']);
save(savefileMat);

%% format and save figures
figure(figureBursts);
axes1 = gca;
set(axes1, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.015 0.03], ...
    'LineWidth', 2.5, 'GridLineStyle', 'none', 'FontSize', 14, ...
    'FontName', 'Arial');
xlabel('time (seconds)', 'FontSize', 14, 'FontName', 'Arial');
ylabel('signal (mV)', 'FontSize', 14, 'FontName', 'Arial');
savefigure = fullfile(outputFolder, [baseFileName 'bursts.fig']);
saveas(figureBursts, savefigure, 'fig');

figure(figureRecording);
axes1 = gca;
set(axes1, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.015 0.03], ...
    'LineWidth', 2.5, 'GridLineStyle', 'none', 'FontSize', 14, ...
    'FontName', 'Arial');
xlabel('time (seconds)', 'FontSize', 14, 'FontName', 'Arial');
ylabel('signal (mV)', 'FontSize', 14, 'FontName', 'Arial');
savefigure = fullfile(outputFolder, [baseFileName 'recording.fig']);
saveas(figureRecording, savefigure, 'fig');

figure(figureSpikeCountRaw);
axes1 = gca;
set(axes1, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.015 0.03], ...
    'LineWidth', 2.5, 'GridLineStyle', 'none', 'FontSize', 14, ...
    'FontName', 'Arial');
xlabel('burst frequency (Hz)', 'FontSize', 14, 'FontName', 'Arial');
ylabel('spike count per burst', 'FontSize', 14, 'FontName', 'Arial');
savefigure = fullfile(outputFolder, [baseFileName 'spikeCount_vs_burstFreq_raw.fig']);
saveas(figureSpikeCountRaw, savefigure, 'fig');

figure(figureSpikeCount);
axes1 = gca;
set(axes1, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.015 0.03], ...
    'LineWidth', 2.5, 'GridLineStyle', 'none', 'FontSize', 14, ...
    'FontName', 'Arial');
xlabel('burst frequency (Hz, 2<f<10.5)', 'FontSize', 14, 'FontName', 'Arial');
ylabel('spike count per burst', 'FontSize', 14, 'FontName', 'Arial');
savefigure = fullfile(outputFolder, [baseFileName 'spikeCount_vs_burstFreq_filtered.fig']);
saveas(figureSpikeCount, savefigure, 'fig');

figure(figureBurstFreqSpikeCountOverTime);
savefigure = fullfile(outputFolder, [baseFileName 'BurstFreqSpikeCountOverTime.fig']);
saveas(figureBurstFreqSpikeCountOverTime, savefigure, 'fig');

%% done
disp('Analysis completed successfully.');
disp(['Results saved to: ' outputFolder]);