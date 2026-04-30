clc;
clear all;
close all;
% =========================================================================
%                  Parameter Configuration Area (Main Script)
% =========================================================================
base_folder = '/Users/weiliu/Documents';

% --- Train Group ---
group1.name = 'Train Group';
group1.folder = fullfile(base_folder, 'same T');
group1.color = [1.0, 0.2, 0.0];             
group1.smoothing_window = ;                
group1.min_peak_amp   = ;                
group1.min_valley_amp = ;              
group1.peak_time_range   = [0, 100];        
group1.valley_time_range = [0, 100];        

% --- Control Group ---
group2.name = 'Control Group';
group2.folder = fullfile(base_folder, 'same C');
group2.color = [0.6392, 0.6471, 0.6549];    % Dark gray
group2.smoothing_window = ;
group2.min_peak_amp   = ;
group2.min_valley_amp = ;
group2.peak_time_range   = [0, 100];
group2.valley_time_range = [0, 100];

file_pattern = '*DLC_mobnet_100_SCI zebrafishMay17shuffle1_200000.csv';

% --- General Algorithm Parameters ---
FRAME_RATE = 240;
APPLY_BASELINE_CORRECTION     = true;    % Apply before finding peak/valley/zero-crossing
APPLY_AMPLITUDE_NORMALIZATION = false;   % Whether to normalize each cycle's amplitude to [-1, 1]
NUM_RESAMPLE_POINTS = 101;               % Resample each cycle to 101 points, i.e., 0-100%

% =========================================================================
%                         Figure Initialization
% =========================================================================
figure_trends = figure( ...
    'Name','Trendline Comparison: Train vs Control', ...
    'NumberTitle','off', ...
    'Position',[100,100,1000,600]);
ax_trends = axes('Parent', figure_trends);
hold(ax_trends,'on');
grid(ax_trends,'on');
title(ax_trends,'Average Trend Comparison: Train vs Control');
xlabel(ax_trends,'Normalized Cycle Progress (%)');
ylabel(ax_trends,'Average Angle (Degrees)');
yline(ax_trends,0,'r--','LineWidth',1.5,'DisplayName','Baseline (0 deg)');

figure_individuals = figure( ...
    'Name','File Average Cycle Comparison for All Samples', ...
    'NumberTitle','off', ...
    'Position',[110,110,1000,600]);
ax_individuals = axes('Parent', figure_individuals);
hold(ax_individuals,'on');
grid(ax_individuals,'on');
title(ax_individuals,'File Average Cycle Comparison for All Samples');
xlabel(ax_individuals,'Normalized Cycle Progress (%)');
ylabel(ax_individuals,'Average Angle (Degrees)');
yline(ax_individuals,0,'r--','LineWidth',1.5,'DisplayName','Baseline (0 deg)');

all_abnormal_files = {};
all_results = {};

% =========================================================================
%                               Main Loop
% =========================================================================
groups = {group1, group2};
for g = 1:numel(groups)
    cg = groups{g};
    fprintf('\n========== Start Processing: %s ==========\n', cg.name);
    [group_abnormal_files, group_results, ~] = process_group_data( ...
        cg, file_pattern, FRAME_RATE, ...
        APPLY_BASELINE_CORRECTION, APPLY_AMPLITUDE_NORMALIZATION, ...
        ax_trends, ax_individuals, NUM_RESAMPLE_POINTS);
    all_abnormal_files = [all_abnormal_files; group_abnormal_files];
    all_results = [all_results; group_results];
end

% =========================================================================
%                         Save Plots/Results
% =========================================================================
legend(ax_trends,'Location','best','Interpreter','none');
hold(ax_trends,'off');
trendlines_plot_path = fullfile(base_folder,'T_vs_C_trendlines_plot.png');
saveas(figure_trends, trendlines_plot_path);
fprintf('\nTrendline comparison plot saved to: %s\n', trendlines_plot_path);

legend(ax_individuals,'Location','best','Interpreter','none');
hold(ax_individuals,'off');
individuals_plot_path = fullfile(base_folder,'T_vs_C_individual_file_mean_cycles_plot.png');
saveas(figure_individuals, individuals_plot_path);
fprintf('All individual file average cycle plots saved to: %s\n', individuals_plot_path);

% Save detailed cycle parameters
results_path = fullfile(base_folder,'T_vs_C_cycle_parameters_detailed.csv');
if ~isempty(all_results)
    results_table = cell2table(all_results, 'VariableNames', ...
        {'Group','Filename','n_cycles', ...
         'mean_p4_sec','mean_p1_sec','mean_p2_sec','mean_p3_sec','mean_T_sec', ...
         'mean_p4_div_T','mean_p1_div_T','mean_p2_div_T','mean_p3_div_T'});
    writetable(results_table, results_path);
    fprintf('Detailed cycle parameters saved to: %s\n', results_path);
else
    results_table = cell2table(cell(0,12), 'VariableNames', ...
        {'Group','Filename','n_cycles', ...
         'mean_p4_sec','mean_p1_sec','mean_p2_sec','mean_p3_sec','mean_T_sec', ...
         'mean_p4_div_T','mean_p1_div_T','mean_p2_div_T','mean_p3_div_T'});
    writetable(results_table, results_path);
    fprintf('No valid cycles extracted, but empty result table generated: %s\n', results_path);
end

% Print abnormal files
if ~isempty(all_abnormal_files)
    fprintf('\n----------------- List of All Abnormal Files -----------------\n');
    fprintf('Detected %d abnormal files:\n', size(all_abnormal_files,1));
    for i = 1:size(all_abnormal_files,1)
        fprintf('Group: %s, Filename: %s\nReason: %s\n\n', ...
            all_abnormal_files{i,1}, ...
            all_abnormal_files{i,2}, ...
            all_abnormal_files{i,3});
    end
    fprintf('--------------------------------------------------\n');
else
    fprintf('\nValid cycles successfully extracted from all files.\n');
end

% =========================================================================
%                         Statistical Summary
% =========================================================================
if ~isempty(all_results) && exist('ttest2','file')==2
    % Columns 9-12 of all_results are now:
    % mean_p4_div_T, mean_p1_div_T, mean_p2_div_T, mean_p3_div_T
    Ratios = cell2mat(all_results(:, 9:12));
    isTrain = strcmp(all_results(:,1),'Train Group');
    metrics = {'mean_p4/T', 'mean_p1/T', 'mean_p2/T', 'mean_p3/T'};
    
    fprintf('\n===== Inter-group Statistical Results: Based on Average Cycle Parameters per File (Mean ± SD) =====\n');
    for k = 1:numel(metrics)
        data_all = Ratios(:,k);
        data_train = data_all(isTrain & ~isnan(data_all));
        data_ctrl  = data_all(~isTrain & ~isnan(data_all));
        p_val = NaN;
        
        if numel(data_train) >= 2 && numel(data_ctrl) >= 2
            [~, p_val] = ttest2(data_train, data_ctrl, 'Vartype','unequal');
        else
            warning('%s: Insufficient sample size, requires >= 2 per group, skipping test.', metrics{k});
        end
        
        fprintf('%s  Train: %.3f ± %.3f (n=%d) | Control: %.3f ± %.3f (n=%d) | p=%.4f\n', ...
            metrics{k}, ...
            mean(data_train,'omitnan'), std(data_train,'omitnan'), numel(data_train), ...
            mean(data_ctrl,'omitnan'), std(data_ctrl,'omitnan'), numel(data_ctrl), ...
            p_val);
    end
else
    fprintf('\nNo results available for statistics, or MATLAB lacks ttest2 function.\n');
end

fprintf('\nScript execution completed.\n');

% =========================================================================
%                       Core Processing Function
% =========================================================================
function [group_abnormal_files, group_results, file_mean_cycles_matrix] = process_group_data( ...
    group, file_pattern, frame_rate, ...
    apply_baseline, apply_amp_norm, ...
    ax_trends_handle, ax_individuals_handle, num_resample_points)
    
    folder_path = group.folder;
    group_name = group.name;
    smoothing_window = group.smoothing_window;
    min_peak_amp = group.min_peak_amp;
    min_valley_amp = group.min_valley_amp;
    peak_time_range = group.peak_time_range;
    valley_time_range = group.valley_time_range;
    
    plot_folder = fullfile(folder_path,'analysis_plots');
    if ~exist(plot_folder,'dir')
        mkdir(plot_folder);
    end
    
    files = dir(fullfile(folder_path, file_pattern));
    if isempty(files)
        fprintf('Warning: No matching files found in path "%s".\n', folder_path);
        group_abnormal_files = {};
        group_results = {};
        file_mean_cycles_matrix = [];
        return;
    end
    
    group_abnormal_files = {};
    group_results = {};
    file_mean_cycles_matrix = [];
    
    for i = 1:numel(files)
        fprintf('Processing file %d/%d: %s\n', i, numel(files), files(i).name);
        file_path = fullfile(folder_path, files(i).name);
        
        try
            raw = readmatrix(file_path);
        catch
            warning('Read failed: %s', files(i).name);
            group_abnormal_files(end+1,:) = {group_name, files(i).name, 'File read failed'};
            continue;
        end
        
        % -----------------------------------------------------------------
        % Extract DLC Coordinates
        % -----------------------------------------------------------------
        data = raw(:,2:end);
        p1  = data(:,1:2);
        p2  = data(:,4:5);
        p3  = data(:,7:8);
        p4  = data(:,10:11);
        p5  = data(:,13:14);
        p6  = data(:,16:17);
        p7  = data(:,19:20);
        p8  = data(:,22:23);
        p9  = data(:,25:26);
        p10 = data(:,28:29);
        p11 = data(:,31:32);
        p12 = data(:,34:35);
        
        all_points = [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12];
        valid_rows = all(~isnan(all_points),2);
        
        if sum(valid_rows) < 20
            group_abnormal_files(end+1,:) = {group_name, files(i).name, 'Too few data rows'};
            continue;
        end
        
        p1  = p1(valid_rows,:);
        p2  = p2(valid_rows,:);
        p3  = p3(valid_rows,:);
        p4  = p4(valid_rows,:);
        p5  = p5(valid_rows,:);
        p6  = p6(valid_rows,:);
        p7  = p7(valid_rows,:);
        p8  = p8(valid_rows,:);
        p9  = p9(valid_rows,:);
        p10 = p10(valid_rows,:);
        p11 = p11(valid_rows,:);
        p12 = p12(valid_rows,:);
        
        % -----------------------------------------------------------------
        % Calculate Angle Curve
        % -----------------------------------------------------------------
        angles_deg = zeros(size(p1,1),6);
        points_for_angle = {p2, p3, p4, p5, p6, p7};
        for k = 1:6
            angles_deg(:,k) = calculate_angle(p1, p12, p1, points_for_angle{k});
        end
        mean_theta = mean(angles_deg,2);
        
        % -----------------------------------------------------------------
        % SG Smoothing
        % -----------------------------------------------------------------
        sg_order = 3;
        sg_frame = 2 * smoothing_window + 1;
        
        if sg_frame > numel(mean_theta)
            sg_frame = 2 * floor((numel(mean_theta)-1)/2) + 1;
        end
        if sg_frame < sg_order + 2
            sg_frame = sg_order + 2;
        end
        if mod(sg_frame,2) == 0
            sg_frame = sg_frame + 1;
        end
        if sg_frame > numel(mean_theta)
            sg_frame = numel(mean_theta);
            if mod(sg_frame,2) == 0
                sg_frame = sg_frame - 1;
            end
        end
        if sg_frame <= sg_order
            group_abnormal_files(end+1,:) = {group_name, files(i).name, 'Valid data too short for SG smoothing'};
            continue;
        end
        
        theta_proc = sgolayfilt(mean_theta, sg_order, sg_frame);
        
        % -----------------------------------------------------------------
        % Baseline correction
        % -----------------------------------------------------------------
        if apply_baseline
            base = linspace(theta_proc(1), theta_proc(end), numel(theta_proc))';
            theta_proc = theta_proc - base;
        end
        
        % -----------------------------------------------------------------
        % Extract all complete cycles
        % -----------------------------------------------------------------
        cycles = extract_all_aligned_cycles(theta_proc);
        
        if isempty(cycles)
            group_abnormal_files(end+1,:) = {group_name, files(i).name, 'Failed to extract complete cycles'};
            continue;
        end
        
        valid_cycle_rows = [];
        valid_resampled_cycles = [];
        
        % Each row of cycles:
        % [zs, zmid, ze, p, v, peak_val, valley_val]
        for c = 1:size(cycles,1)
            zs   = cycles(c,1);
            zmid = cycles(c,2);
            ze   = cycles(c,3);
            p    = cycles(c,4);
            v    = cycles(c,5);
            peak_val   = cycles(c,6);
            valley_val = cycles(c,7);
            
            % Amplitude threshold filtering
            if ~(peak_val > min_peak_amp && valley_val < min_valley_amp)
                continue;
            end
            
            peak_time_percent   = ((p - zs) / (ze - zs)) * 100;
            valley_time_percent = ((v - zs) / (ze - zs)) * 100;
            
            % Peak/Valley occurrence time range filtering
            if ~(peak_time_percent >= peak_time_range(1) && peak_time_percent <= peak_time_range(2) && ...
                 valley_time_percent >= valley_time_range(1) && valley_time_percent <= valley_time_range(2))
                continue;
            end
            
            T  = (ze - zs) / frame_rate;
            p4 = (p - zs) / frame_rate;
            p1 = (zmid - p) / frame_rate;
            p2 = (v - zmid) / frame_rate;
            p3 = (ze - v) / frame_rate;
            
            p4_div_T = p4 / T;
            p1_div_T = p1 / T;
            p2_div_T = p2 / T;
            p3_div_T = p3 / T;
            
            valid_cycle_rows = [valid_cycle_rows; ...
                p4, p1, p2, p3, T, ...
                p4_div_T, p1_div_T, p2_div_T, p3_div_T];
            
            % -------------------------------------------------------------
            % Resample this cycle for trend plotting
            % -------------------------------------------------------------
            cycle_data = theta_proc(floor(zs):ceil(ze));
            if apply_amp_norm
                ma = max(abs(cycle_data));
                if ma > 1e-6
                    cycle_data = cycle_data / ma;
                end
            end
            
            % Bring start and end of cycle back to the same baseline to avoid unclosed ends
            cycle_data = cycle_data - linspace(cycle_data(1), cycle_data(end), numel(cycle_data))';
            
            original_x = linspace(0,100,numel(cycle_data));
            standard_x = linspace(0,100,num_resample_points);
            
            resampled_cycle = interp1(original_x, cycle_data, standard_x, 'spline');
            valid_resampled_cycles = [valid_resampled_cycles; resampled_cycle];
        end
        
        if isempty(valid_cycle_rows)
            group_abnormal_files(end+1,:) = {group_name, files(i).name, 'All cycles failed amplitude or time window filtering'};
            continue;
        end
        
        % -----------------------------------------------------------------
        % Average all valid cycles within each file
        % -----------------------------------------------------------------
        mean_cycle_values = mean(valid_cycle_rows, 1, 'omitnan');
        n_cycles = size(valid_cycle_rows, 1);
        
        group_results = [group_results; ...
            {group_name, files(i).name, n_cycles, ...
             mean_cycle_values(1), mean_cycle_values(2), mean_cycle_values(3), mean_cycle_values(4), mean_cycle_values(5), ...
             mean_cycle_values(6), mean_cycle_values(7), mean_cycle_values(8), mean_cycle_values(9)}];
         
        % Average cycle curve for each file
        file_mean_cycle = mean(valid_resampled_cycles, 1, 'omitnan');
        file_mean_cycles_matrix(end+1,:) = file_mean_cycle;
        
        % -----------------------------------------------------------------
        % Save single file analysis plot
        % -----------------------------------------------------------------
        fig_single = figure( ...
            'Name', files(i).name, ...
            'NumberTitle','off', ...
            'Visible','off');
        ax = axes('Parent', fig_single);
        hold(ax,'on');
        grid(ax,'on');
        
        plot(ax, theta_proc, 'b-', 'DisplayName','Preprocessed Angle');
        yline(ax,0,'r--','DisplayName','Baseline');
        
        % All original zero crossings
        zc_idx = find(theta_proc(1:end-1).*theta_proc(2:end) <= 0);
        if ~isempty(zc_idx)
            zc_pos = zc_idx + abs(theta_proc(zc_idx)) ./ ...
                (abs(theta_proc(zc_idx)) + abs(theta_proc(zc_idx+1)));
            scatter(ax, zc_pos, zeros(size(zc_pos)), 18, 'g', 'filled', ...
                'DisplayName','Zero Crossings');
        end
        
        % Mark valid cycle peaks and valleys
        valid_peaks = [];
        valid_valleys = [];
        for c = 1:size(cycles,1)
            zs   = cycles(c,1);
            zmid = cycles(c,2);
            ze   = cycles(c,3);
            p    = cycles(c,4);
            v    = cycles(c,5);
            peak_val   = cycles(c,6);
            valley_val = cycles(c,7);
            
            if ~(peak_val > min_peak_amp && valley_val < min_valley_amp)
                continue;
            end
            
            peak_time_percent   = ((p - zs) / (ze - zs)) * 100;
            valley_time_percent = ((v - zs) / (ze - zs)) * 100;
            if ~(peak_time_percent >= peak_time_range(1) && peak_time_percent <= peak_time_range(2) && ...
                 valley_time_percent >= valley_time_range(1) && valley_time_percent <= valley_time_range(2))
                continue;
            end
            
            valid_peaks = [valid_peaks; p, peak_val];
            valid_valleys = [valid_valleys; v, valley_val];
        end
        
        if ~isempty(valid_peaks)
            plot(ax, valid_peaks(:,1), valid_peaks(:,2), ...
                'k^','MarkerFaceColor','k', ...
                'LineStyle','none', ...
                'DisplayName','Valid Peak');
        end
        if ~isempty(valid_valleys)
            plot(ax, valid_valleys(:,1), valid_valleys(:,2), ...
                'rv','MarkerFaceColor','r', ...
                'LineStyle','none', ...
                'DisplayName','Valid Valley');
        end
        
        title(ax, ['File: ', files(i).name, ...
            sprintf(' | Valid Cycles: %d', n_cycles)], ...
            'Interpreter','none');
        xlabel(ax,'Valid Frames');
        ylabel(ax,'Average Angle (Degrees)');
        legend(ax,'Location','best','Interpreter','none');
        
        [~, baseFileName, ~] = fileparts(files(i).name);
        saveas(fig_single, fullfile(plot_folder,[baseFileName,'.png']));
        close(fig_single);
    end
    
    % ---------------------------------------------------------------------
    % Plot group average trend
    % ---------------------------------------------------------------------
    if ~isempty(file_mean_cycles_matrix)
        standard_x = linspace(0,100,num_resample_points);
        mean_trend = mean(file_mean_cycles_matrix,1,'omitnan');
        sem_trend  = std(file_mean_cycles_matrix,0,1,'omitnan') / sqrt(size(file_mean_cycles_matrix,1));
        
        upper_bound = mean_trend + sem_trend;
        lower_bound = mean_trend - sem_trend;
        
        fill(ax_trends_handle, ...
            [standard_x, fliplr(standard_x)], ...
            [upper_bound, fliplr(lower_bound)], ...
            group.color, ...
            'FaceAlpha',0.2, ...
            'EdgeColor','none', ...
            'HandleVisibility','off');
            
        plot(ax_trends_handle, ...
            standard_x, mean_trend, ...
            'Color', group.color, ...
            'LineWidth',3.0, ...
            'DisplayName', [group.name, ' Average Trend']);
            
        % Average cycle curve for each file
        for k = 1:size(file_mean_cycles_matrix,1)
            if k == 1
                plot(ax_individuals_handle, ...
                    standard_x, file_mean_cycles_matrix(k,:), ...
                    'Color', group.color, ...
                    'LineWidth',0.8, ...
                    'DisplayName', group.name);
            else
                plot(ax_individuals_handle, ...
                    standard_x, file_mean_cycles_matrix(k,:), ...
                    'Color', group.color, ...
                    'LineWidth',0.8, ...
                    'HandleVisibility','off');
            end
        end
    end
end

% =========================================================================
%       Extract all complete cycles "0 -> Peak -> 0 -> Valley -> 0"
% =========================================================================
function cycles = extract_all_aligned_cycles(theta_proc)
    % Each row of cycles:
    % [zs, zmid, ze, p, v, peak_val, valley_val]
    %
    % zs   = Cycle start zero crossing
    % zmid = Middle zero crossing
    % ze   = Cycle end zero crossing
    % p    = Positive peak position
    % v    = Negative valley position
    
    cycles = [];
    if numel(theta_proc) < 4
        return;
    end
    
    zc_idx = find(theta_proc(1:end-1).*theta_proc(2:end) <= 0);
    if numel(zc_idx) < 3
        return;
    end
    
    zc_pos = zc_idx + abs(theta_proc(zc_idx)) ./ ...
        (abs(theta_proc(zc_idx)) + abs(theta_proc(zc_idx+1)));
        
    [~, peak_locs_all]   = findpeaks(theta_proc);
    [~, valley_locs_all] = findpeaks(-theta_proc);
    
    if isempty(peak_locs_all) || isempty(valley_locs_all)
        return;
    end
    
    % Iterate through all three consecutive zero crossings:
    % zs -> zmid -> ze
    % Require positive peak between zs-zmid, negative valley between zmid-ze
    for zi = 1:(numel(zc_pos)-2)
        zs   = zc_pos(zi);
        zmid = zc_pos(zi+1);
        ze   = zc_pos(zi+2);
        
        peak_candidates = peak_locs_all( ...
            peak_locs_all > zs & peak_locs_all < zmid);
        valley_candidates = valley_locs_all( ...
            valley_locs_all > zmid & valley_locs_all < ze);
            
        if isempty(peak_candidates) || isempty(valley_candidates)
            continue;
        end
        
        % Select the highest peak in this positive half-cycle
        [peak_val, peak_idx] = max(theta_proc(peak_candidates));
        p = peak_candidates(peak_idx);
        
        % Select the lowest valley in this negative half-cycle
        [valley_val, valley_idx] = min(theta_proc(valley_candidates));
        v = valley_candidates(valley_idx);
        
        if peak_val > 0 && valley_val < 0
            cycles = [cycles; zs, zmid, ze, p, v, peak_val, valley_val];
        end
    end
end

% =========================================================================
%                            Auxiliary Functions
% =========================================================================
function angle = calculate_angle(p1, p12, p1_other, p2_other)
    vec_ref  = p12 - p1;
    vec_body = p2_other - p1_other;
    
    dotp = sum(vec_ref .* vec_body, 2);
    mag_ref  = sqrt(sum(vec_ref.^2, 2));
    mag_body = sqrt(sum(vec_body.^2, 2));
    
    denom = mag_ref .* mag_body;
    denom(denom == 0) = 1e-6;
    
    c = dotp ./ denom;
    c = max(min(c,1),-1);
    ang = acosd(c);
    
    cross_val = vec_ref(:,1).*vec_body(:,2) - vec_ref(:,2).*vec_body(:,1);
    sgn = sign(cross_val);
    sgn(sgn == 0) = 1;
    
    angle = ang .* sgn;
end