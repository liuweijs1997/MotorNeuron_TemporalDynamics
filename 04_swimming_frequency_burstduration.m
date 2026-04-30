% analyze_swimming_extra_burst_dutycycle
% Detect bursts from an extracellular ABF trace.
%
% Definitions:
% cycleDurations = time from one burst onset to the next burst onset
% burstDurations = firing duration within one burst
% burstDutyCycle = burstDurations / cycleDurations
% burstFrequencies = 1 / cycleDurations
%
% Filters:
% keep only bursts with:
% 2 < burst frequency < 10.5 Hz
% burst duration / cycle > 0.1

clear all;
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

figure;
channelNames = h.recChNames(channelsToPlot);

for i = 1:length(channelsToPlot)
    ax(i) = subplot(length(channelsToPlot), 1, i);
    plot(time, data(:, channelsToPlot(i)));
    title(strcat(channelNames{i}, ' channel# ', num2str(channelsToPlot(i))));
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

close gcf; % close preview figure

%% choose analyzed trace
data = data(:, traceNumber);

%% create output folder
baseFileName = fileName(1:end-4); % remove .abf
outputFolder = fullfile(pathName, baseFileName);

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% plot raw + filtered trace
figureBursts = figure;
h1 = plot(time, data);
set(h1, 'DisplayName', 'recording');
xlabel('time (seconds)');
ylabel('signal (mV)');
hold on;

%% preprocess extracellular signal
% sliding DC removal
filteredData = data - slidingavg(data, round(0.025 / samplingInterval));

% rectify
filteredData = abs(filteredData);

% smooth
filteredData = slidingavg(filteredData, round(0.025 / samplingInterval));

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

burstDurationsForPlot = burstDurations(1:end-1);
burstDutyCycle = burstDurationsForPlot ./ cycleDurations;

cycleTimes = burstStartTimes(1:end-1);

%% apply filters
% keep only:
% 2 < burst frequency < 10.5 Hz
% burst duration / cycle > 0.1
validIdx = burstFrequencies > 2 & burstFrequencies < 10.5 & burstDutyCycle > 0.1;

cycleDurations = cycleDurations(validIdx);
burstFrequencies = burstFrequencies(validIdx);
burstDurationsForPlot = burstDurationsForPlot(validIdx);
burstDutyCycle = burstDutyCycle(validIdx);
cycleTimes = cycleTimes(validIdx);

if isempty(burstFrequencies)
    error('No bursts remained after applying the filters (2 < frequency < 10 Hz and duty cycle > 0.1).');
end

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

%% mean duty cycle at 4-5 Hz
indexFourToFiveHz = find(4 <= burstFrequencies & burstFrequencies <= 5);

if isempty(indexFourToFiveHz)
    meanDutyCycleFourToFiveHz = NaN;
else
    meanDutyCycleFourToFiveHz = mean(burstDutyCycle(indexFourToFiveHz));
end

%% plot duty cycle vs frequency
figureDutyCycle = figure;
scatter(burstFrequencies, burstDutyCycle, ...
    'MarkerEdgeColor', [0 0 0], 'LineWidth', 2);
xlabel('burst frequency (Hz, 2<f<10)');
ylabel('burst duration / cycle');
title('Burst duty cycle vs burst frequency (2<f<10, duty>0.1)');
% ylim([0 1]); % uncomment if needed

%% plot frequency and duty cycle over time
figureBurstFreqDutyOverTime = figure;

yyaxis left
scatter(cycleTimes, burstFrequencies, 20, 'k', 'filled');
ylabel('burst frequency (Hz)', 'FontSize', 14, 'FontName', 'Arial');

yyaxis right
scatter(cycleTimes, burstDutyCycle, 20, ...
    'MarkerEdgeColor', [0.5 0.5 0.5]);
ylabel('burst duration / cycle', 'FontSize', 14, 'FontName', 'Arial', ...
    'Color', [0.5 0.5 0.5]);

xlabel('time (s)', 'FontSize', 14, 'FontName', 'Arial');

ax = gca;
set(ax, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.015 0.03], ...
    'LineWidth', 2.5, 'GridLineStyle', 'none', 'FontSize', 14, ...
    'FontName', 'Arial');

%% save txt data
savefileTxt = fullfile(outputFolder, [baseFileName 'analysis.txt']);

saveMatrix = [ ...
    cycleTimes, ...
    cycleTimes + cycleDurations, ...
    cycleDurations, ...
    burstFrequencies, ...
    burstDurationsForPlot, ...
    burstDutyCycle];

myfile = fopen(savefileTxt, 'wt');

fprintf(myfile, 'burst threshold:\t%4.6f\r\n', burstThreshold);
fprintf(myfile, 'minimum burst duration (s):\t%4.6f\r\n', minBurstDuration);
fprintf(myfile, 'filters:\t2 < burst frequency < 10 Hz; burst duration/cycle > 0.1\r\n');
fprintf(myfile, 'mean duty cycle at 4-5 Hz:\t%4.6f\r\n', meanDutyCycleFourToFiveHz);
fprintf(myfile, ['burst start time (s)\t next burst start time (s)\t' ...
    'cycle duration (s)\t burst frequency (Hz)\t burst duration (s)\t burst duration/cycle\r\n']);

fclose(myfile);

save(savefileTxt, 'saveMatrix', '-ascii', '-tabs', '-append');

%% save matlab variables
analysisParametersFile = fullfile(outputFolder, 'analysis_parameters.mat');
save(analysisParametersFile, ...
    'burstStartTimes', 'burstEndTimes', 'burstDurations', ...
    'cycleDurations', 'burstFrequencies', 'burstDurationsForPlot', ...
    'burstDutyCycle', 'cycleTimes', ...
    'burstThreshold', 'minBurstDuration');
%% save burst frequency and burst duty cycle to csv
csvFile = fullfile(outputFolder, [baseFileName '_burstFrequency_DutyCycle.csv']);

T = table(burstFrequencies(:), burstDutyCycle(:), ...
    'VariableNames', {'burstFrequency_Hz', 'burstDuration_over_cycle'});

writetable(T, csvFile);

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

figure(figureDutyCycle);
axes1 = gca;
set(axes1, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.015 0.03], ...
    'LineWidth', 2.5, 'GridLineStyle', 'none', 'FontSize', 14, ...
    'FontName', 'Arial');
xlabel('burst frequency (Hz, 2<f<10)', 'FontSize', 14, 'FontName', 'Arial');
ylabel('burst duration / cycle', 'FontSize', 14, 'FontName', 'Arial');

savefigure = fullfile(outputFolder, [baseFileName 'burstDutyCycleVsFrequency.fig']);
saveas(figureDutyCycle, savefigure, 'fig');

figure(figureBurstFreqDutyOverTime);
savefigure = fullfile(outputFolder, [baseFileName 'BurstFreqDutyOverTime.fig']);
saveas(figureBurstFreqDutyOverTime, savefigure, 'fig');

%% clean up
clear('AX','H1','H2','a','ans','axes1','ax','channelNames','dBurst', ...
    'figureBursts','figureDutyCycle','figureBurstFreqDutyOverTime', ...
    'figureRecording','filterIndex','h1','h2','h3','i','indexFourToFiveHz', ...
    'isBurst','lengthOfData','maxTime','myfile','num_lines','options', ...
    'prompt','savefigure','savefileMat','savefileTxt','saveMatrix', ...
    'samplingFrequency','samplingInterval','thresholdFile','minDurFile');