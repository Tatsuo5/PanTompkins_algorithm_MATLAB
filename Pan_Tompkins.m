function [peaks_idx, peaks_data, THRESHOLD_I1, smoothed_data] = Pan_Tompkins(fs, data, ecg_column, k, plot_flag)
% Pan_Tompkins
% Input:
%   fs         - Sampling frequency (Hz)
%   data       - Raw data matrix
%   ecg_column - Column number containing the ECG signal in the data
%   k          - Condition to exclude peaks from threshold update if they exceed k times the threshold
%   plot_flag  - 1 to generate plots, 0 to not generate plots
% Output:
%   peaks_idx - Detected R-wave indices
%   peaks_data - Amplitudes (magnitudes) of the detected R-waves
%   smoothed_data - Preprocessed ECG signal
%   THRESHOLD_I1 - Moving threshold signal (first threshold, with the second threshold being half of this)
% When block processing is desired, peaks_idx is output so that R-waves can be determined for a given time (fs*time) interval

%% Data Preprocessing
data = data(:, ecg_column); % Extract ECG data from the specified column
data = data(isfinite(data)); % Remove NaN values
data = data(:);
% Apply bandpass filter
[b, a] = butter(5, [5 / (fs / 2), 15 / (fs / 2)], 'bandpass'); % Filter coefficients
filtered_data = filtfilt(b, a, data);
% Differentiation filter
derivative_filter = [-1 -2 0 2 1] .* (1 / 8) * fs;% Modified to match the original paper notation [-1 -2 0 2 1]; the previous implementation [1 2 0 -2 -1] was equally correct (only sign inverted).
differentiated_data = filtfilt(derivative_filter, 1, filtered_data);
% Squaring and moving average
squared_data = differentiated_data .^ 2;
smoothed_data = movmean(squared_data, round(0.15 * fs));

%% Initial Threshold Setting
refractory_period = round(0.2 * fs); % Refractory period
initial_segment = smoothed_data(1:2 * fs); % Learning phase
[pks, ~] = findpeaks(initial_segment, 'MinPeakDistance', refractory_period);
initial_SPKI = mean(pks);
initial_NPKI = mean(pks(pks < prctile(pks, 50)));

%% Peak Detection and Threshold Setting
% Peak detection and real-time re-detection
signal_peaks_idx = [];
SPKI = initial_SPKI;
NPKI = initial_NPKI;
retrieved_peaks_idx = []; % Array to store re-detected R-wave peaks
signal_pks_data = [];
t_wave_detected = [];
t_wave_amplitudes = [];
THRESHOLD_I1 = zeros(1, length(smoothed_data)); % Signal for threshold
threshold_I1 = NPKI + 0.25 * (SPKI - NPKI);
threshold_I2 = 0.5 * threshold_I1;
% Peak detection (R-wave candidates)
[pks, locs] = findpeaks(smoothed_data, 'MinPeakDistance', refractory_period); 
% Real-time re-detection and T-wave identification with threshold update for each index
for i = 1:length(locs)
    % Execute re-detection
    if ~isempty(signal_peaks_idx)
        % Calculate RR intervals
        RR_interval = diff(signal_peaks_idx);
        % RRaverage1: Average RR interval considering regular heart rhythm
        RRaverage1 = mean(RR_interval(max(1, end-7):end));
        % RRaverage2: Average RR interval within the range (92% to 116%)
        valid_RR = RR_interval(RR_interval >= 0.92 * RRaverage1 & RR_interval <= 1.16 * RRaverage1);
        if ~isempty(valid_RR)
            RRaverage2 = mean(valid_RR);
        else
            RRaverage2 = RRaverage1; % Substitute with the regular average
        end
        % Missed beat check: if no R-wave within 166%, re-detect
        if (locs(i) - signal_peaks_idx(end)) > 1.66 * RRaverage2
            search_start = signal_peaks_idx(end) + refractory_period;
            search_end = locs(i) - refractory_period;
            range_length = search_end - search_start + 1; % Calculate the length of the re-detection range
            re_locs = [];
            if range_length > refractory_period
                % If the re-detection range is larger than the refractory period
                [re_pks, re_locs] = findpeaks(smoothed_data(search_start:search_end), 'MinPeakDistance', refractory_period);
            end
            % Detect peaks in the re-detection range
            if ~isempty(re_locs)
                % Determine peaks above the second threshold
                above_threshold_idx = find(re_pks >= threshold_I2);
                % If there are peaks above the second threshold, process them all
                for j = above_threshold_idx'
                    current_peak = re_pks(j);
                    current_loc = re_locs(j) + search_start;
                    is_t_wave = false;
                    if ~isempty(signal_peaks_idx)
                        last_r_peak = signal_peaks_idx(end);
                        if (current_loc - last_r_peak) < round(0.36 * fs)
                            window_size = round(0.075 * fs);
                            if last_r_peak > window_size
                                prev_slope = max(diff(smoothed_data(last_r_peak-window_size:last_r_peak)));
                            else
                                prev_slope = max(diff(smoothed_data(1:last_r_peak)));
                            end
                            if current_loc > window_size
                                curr_slope = max(diff(smoothed_data(current_loc-window_size:current_loc)));
                            else
                                curr_slope = max(diff(smoothed_data(1:current_loc)));
                            end
                            if curr_slope < 0.5 * prev_slope
                                t_wave_detected = [t_wave_detected; current_loc];
                                t_wave_amplitudes = [t_wave_amplitudes; current_peak];
                                is_t_wave = true;
                            end
                        end
                    end
                    if is_t_wave
                        NPKI = 0.125 * current_peak + 0.875 * NPKI;
                    else
                        retrieved_peaks_idx = [retrieved_peaks_idx; current_loc];
                        signal_peaks_idx = [signal_peaks_idx; current_loc];
                        signal_pks_data = [signal_pks_data; current_peak];
                        SPKI = 0.125 * current_peak + 0.875 * SPKI;
                    end
                    threshold_I1 = NPKI + 0.25 * (SPKI - NPKI);
                    threshold_I2 = 0.5 * threshold_I1;
                end
            end
        end
    end
    % T-wave flag
    is_t_wave = false;
    % T-wave identification (comparing the current peak with the previous peak)
    if ~isempty(signal_peaks_idx)
        last_r_peak = signal_peaks_idx(end); % Last detected R-wave peak
        if (locs(i) - last_r_peak) < round(0.36 * fs)
            % Compare slopes to determine T-wave
            window_size = round(0.075 * fs); % 75ms window
            % Calculate the slope of the previous peak (preventing out-of-bound errors)
            if last_r_peak > window_size
                prev_slope = max(diff(smoothed_data(last_r_peak-window_size:last_r_peak)));
            else
                prev_slope = max(diff(smoothed_data(1:last_r_peak))); % Use from the start of data to current position
            end
            % Calculate the slope of the current peak (preventing out-of-bound errors)
            if locs(i) > window_size
                curr_slope = max(diff(smoothed_data(locs(i)-window_size:locs(i))));
            else
                curr_slope = max(diff(smoothed_data(1:locs(i)))); % Use from the start of data to current position
            end
            % T-wave determination
            if curr_slope < 0.5 * prev_slope
                t_wave_detected = [t_wave_detected; locs(i)]; % Peak identified as T-wave
                t_wave_amplitudes = [t_wave_amplitudes; pks(i)]; % Amplitude of the peak identified as T-wave
                is_t_wave = true;
            end
        end
    end
    % Peak classification (if T-wave, treat as noise)
    if is_t_wave
        NPKI = 0.125 * pks(i) + 0.875 * NPKI;
    elseif pks(i) > k * threshold_I1
        % Peak exceeding k times: detected as peak but not used for threshold update
        signal_peaks_idx = [signal_peaks_idx; locs(i)];
        signal_pks_data = [signal_pks_data; pks(i)];
    elseif pks(i) > threshold_I1
        % Normal R-wave peak
        SPKI = 0.125 * pks(i) + 0.875 * SPKI;
        signal_peaks_idx = [signal_peaks_idx; locs(i)];
        signal_pks_data = [signal_pks_data; pks(i)];
    else
        % Noise peak
        NPKI = 0.125 * pks(i) + 0.875 * NPKI;
    end 
    % Update thresholds
    threshold_I1 = NPKI + 0.25 * (SPKI - NPKI);
    threshold_I2 = 0.5 * threshold_I1;
    % Reflect threshold on signal
    if i == 1
        THRESHOLD_I1(1:locs(i)) = threshold_I1;
    elseif i == length(locs)
        THRESHOLD_I1(locs(i-1)+1:end) = threshold_I1;
    else
        THRESHOLD_I1(locs(i-1)+1:locs(i)) = threshold_I1;
    end
end
% Save data after re-detection and T-wave identification
peaks_idx = signal_peaks_idx;
peaks_data = smoothed_data(peaks_idx);

%% Plot
if plot_flag == 1
    time = (0:length(smoothed_data) - 1) / fs; % Create time axis
    figure;
    plot(time, smoothed_data, 'b', 'DisplayName', 'Smoothed Signal'); % Plot smoothed signal
    hold on;
    % Plot re-detected R-waves
    if ~isempty(retrieved_peaks_idx)
        valid_indices = retrieved_peaks_idx(retrieved_peaks_idx > 0 & retrieved_peaks_idx <= length(smoothed_data)); 
        retrieved_amplitudes = smoothed_data(valid_indices-1); % Adjusted for time axis (since indexing starts at 0)
        plot((valid_indices-2) / fs, retrieved_amplitudes, 'go', 'MarkerSize', 6, 'LineWidth', 1, 'DisplayName', 'R-peaks (Re-detected)');
    end
    % Plot R-waves detected by the first threshold
    if ~isempty(peaks_idx)
        initial_peaks = setdiff(peaks_idx, retrieved_peaks_idx); % Exclude re-detected peaks
        initial_peaks = initial_peaks(initial_peaks > 0 & initial_peaks <= length(smoothed_data)); 
        if ~isempty(initial_peaks)
            plot((initial_peaks-1) / fs, smoothed_data(initial_peaks), 'ro', 'MarkerSize', 6, 'LineWidth', 1, 'DisplayName', 'R-peaks (Initial)');
        end
    end
    % Plot T-waves above the second threshold
    if ~isempty(t_wave_detected)
        valid_t_wave_times = [];
        valid_t_wave_amplitudes = [];
        for i = 1:length(t_wave_detected)
            if t_wave_amplitudes(i) > THRESHOLD_I1(t_wave_detected(i)) / 2
                valid_t_wave_times = [valid_t_wave_times, t_wave_detected(i)];
                valid_t_wave_amplitudes = [valid_t_wave_amplitudes, t_wave_amplitudes(i)];
            end
        end
        if ~isempty(valid_t_wave_times)
            plot((valid_t_wave_times-1) / fs, valid_t_wave_amplitudes, 'kx', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'T-wave (Excluded)');
        end
    end
    % Plot thresholds
    plot(time, THRESHOLD_I1, 'm--', 'LineWidth', 1, 'DisplayName', 'Threshold I1');
    plot(time, THRESHOLD_I1 / 2, 'g--', 'LineWidth', 1, 'DisplayName', 'Threshold I2');
    title('ECG Signal with Detected R-peaks and T-wave Exclusion');
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend('show', 'Location', 'best');
    grid on;
    hold off;
end
end
