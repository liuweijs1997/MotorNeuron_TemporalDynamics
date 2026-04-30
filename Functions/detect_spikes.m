function[spikeTimes]= detect_spikes(data,si)

% detect_spikes - Spike detection from intracellular/extracellular recordings
%
% This function detects spike events from a continuous voltage trace and
% returns the corresponding spike times. The detection pipeline includes:
%   (1) High-pass filtering to remove low-frequency components
%   (2) Removal of line noise (e.g., 50 Hz interference)
%   (3) Detection of candidate spike events
%   (4) Manual or semi-automatic selection of true spikes
%
% Inputs:
%   data - Continuous voltage trace (mV)
%   si   - Sampling interval (microseconds)
%
% Output:
%   spikeTimes - Detected spike times (seconds)
%
% -------------------------------------------------------------------------
% Signal processing steps:
%
%   1. High-pass filtering:
%      Removes low-frequency components below 50 Hz.
%
%   2. Line noise removal:
%      Suppresses periodic noise (e.g., 50 Hz electrical interference).
%
%   3. Spike candidate detection:
%      Identifies potential spike-like events using threshold-based methods.
%
%   4. Spike selection:
%      Refines candidate events to retain biologically relevant spikes.
%
% -------------------------------------------------------------------------
% Dependencies:
%
%   This function relies on MATLAB built-in functions and additional
%   third-party functions, including:
%
%   - butterhigh1        (high-pass filter design)
%   - templatefilter     (line noise removal)
%   - suc2spike          (spike candidate detection)
%   - selectspike        (spike selection/refinement)
%
%   These functions originate from the Wagenaar Lab toolbox and/or
%   Adam Taylor's analysis code.
%
% -------------------------------------------------------------------------
% Usage in this study:
%
%   This function was used to extract spike timing information from
%   electrophysiological recordings, which were subsequently used for
%   burst detection, frequency analysis, and temporal precision evaluation.
%
%   Study:
%   "Agility Training Enhances Motor Temporal Precision by Reweighting
%    Spinal Phase-Locked Commissural Inhibition"
%
% -------------------------------------------------------------------------
% Adaptation and integration:
%   Integrated into analysis pipeline by Wei Liu, 2026
%
% -------------------------------------------------------------------------
% Note:
%   This function includes third-party components and is used for
%   research purposes only. All rights belong to the original authors.



% calculate the time vector
samplingInterval = si/1000000;                         % Sampling interval in s
samplingFrequency = 1/samplingInterval;                 % Sampling frequency
lengthOfData = length(data);                           % Length of signal
time = (0:lengthOfData-1)*samplingInterval;            % Time vector in s


% plot voltage trace for visual control
figure;
%subplot(2,1,1);
plot(time,data)
xlabel('time (seconds)')
ylabel('Vm amp (mV)')


% % ask for the min and max values to zoom into the data
% prompt = {'lower margin (s)','upper margin (s)'};
% dlg_title = 'enter margins for magnified trace';
% num_lines = 1;
% def = {'0','120'};
% answer = inputdlg(prompt,dlg_title,num_lines,def);
% 
% xMin=str2double(answer{1});
% xMax=str2double(answer{2});
% 
% subplot(2,1,2);
% plot(time,data)
% xlabel('time (seconds)')
% ylabel('Vm amp (mV)')
% axis([xMin xMax -Inf Inf]);



%removing all frequencies below 50Hz
cutoffFrequencyHz = 50;
[b,a]=butterhigh1(cutoffFrequencyHz/samplingFrequency);
filteredData=filtfilt(b,a,data);




% % plot filtered data for visual control
% figure;
% subplot(2,1,1);
% plot(time,filteredData)
% xlabel('time (seconds)')
% ylabel('Vm filtered (mV)')
% subplot(2,1,2);
% plot(time,filteredData)
% xlabel('time (seconds)')
% ylabel('Vm amp (mV)')
% axis([xMin xMax -Inf Inf]);


% remove 50Hz noise
filteredData2 = templatefilter(filteredData,samplingFrequency/50,200/samplingFrequency,50);


% % plot filtered data without line pickup
% figure;
% subplot(2,1,1);
% plot(time,filteredData2)
% xlabel('time (seconds)')
% ylabel('Vm amp (mV)')
% subplot(2,1,2);
% plot(time,filteredData2)
% xlabel('time (seconds)')
% ylabel('Vm filtered (mV)')
% axis([xMin xMax -Inf Inf]);

%detect all spike-like events
potentialSpikes=suc2spike(filteredData2,time,2,2,1);


% plot preselected spikes
figureSpikes=figure;
plot(time,filteredData2);
hold on
plot(potentialSpikes.tms,potentialSpikes.amp,'k.');
hold on
xlabel('time (seconds)')
ylabel('Vm amp (mV)')


%select the real spikes out of the potential ones
spikes=selectspike(potentialSpikes);

% plot selected spikes in same figure
figure(figureSpikes);
hold on;
plot(spikes.tms,spikes.amp,'r.');

% output-------------------------------------------------------------------
spikeTimes=spikes.tms;











