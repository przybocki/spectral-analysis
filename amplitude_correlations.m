function [tau_max] = amplitude_correlations(x,t,band1,band2,bw,filt_order,lp_cutoff,trunc, N, M, max_lag,displ)

% MATLAB function amplitude_correlations determines the direction of energy flow between frequency bands using the amplitude correlation technique
%
% Inputs:
%
% x: time series data from anode segment
% t: time vector
% band1: center frequency of lower frequency band 
% band2: center frequency of higher frequency band 
% bw: bandwidth 
% filt_order: bandpass filter order (default = 10)
% lp_cutoff: low pass filter cutoff (default = bandwidth)
% trunc: number of samples to truncate from signal (default = 0)
% N: number of samples in each window
% M: fraction of overlap for adjacent blocks (default = 0.9)
% max_lag: maximum number of samples in cross-correlation (default = N)
% displ: if true, plots correlation function
%
% Output:
%
% tau_max = time lag of maximum correlation K(tau) between frequency bands
% Syntax K(tau) = < [f1(t)^2] [f2(t+tau)^2] >
%
% 2025 Ryan Przybocki

clc;
close all;

dt = t(2)-t(1);
fs = 1/dt;

% Filter the time series into two separate frequency bands
low_band = [band1-bw/2 band1+bw/2];     % Low frequency band
high_band = [band2-bw/2 band2+bw/2];    % High frequency band

low_band_str = string(band1/1e6) + '+/-' + string(bw/2/1e6) + ' MHz';
high_band_str = string(band2/1e6) + '+/-' + string(bw/2/1e6) + ' MHz';


% Design bandpass filters using Butterworth filters:
bpFilt_low = designfilt('bandpassiir', 'FilterOrder', filt_order, ...
                        'HalfPowerFrequency1', low_band(1), 'HalfPowerFrequency2', low_band(2), ...
                        'SampleRate', fs, 'DesignMethod', 'butter');

bpFilt_high = designfilt('bandpassiir', 'FilterOrder', filt_order, ...
                         'HalfPowerFrequency1', high_band(1), 'HalfPowerFrequency2', high_band(2), ...
                         'SampleRate', fs, 'DesignMethod', 'butter');

% Apply the filters to the time series
x_low = filtfilt(bpFilt_low, x);  % Low-frequency filtered data
x_high = filtfilt(bpFilt_high, x); % High-frequency filtered data

% Square the filtered signals
x_low_squared = x_low .^ 2;   
x_high_squared = x_high .^ 2; 

% Design low-pass filter:
lpFilt = designfilt('lowpassiir', 'FilterOrder', filt_order, ...
                    'HalfPowerFrequency', lp_cutoff, ...
                    'SampleRate', fs, 'DesignMethod', 'butter');

% Apply the low-pass filter to the squared signals
x_low_envelope = filtfilt(lpFilt, x_low_squared);
x_high_envelope = filtfilt(lpFilt, x_high_squared);

% Truncate signal if desired
x_low_envelope = x_low_envelope(trunc+1:end-trunc);
x_high_envelope = x_high_envelope(trunc+1:end-trunc);
t = t(trunc+1:end-trunc);
L = length(t);

% Break time series into overlapping blocks
overlap = round(N * M);                 % Number of overlapping points
step = N - overlap;                     % Step size between blocks
N_block = floor((L - overlap) / step);  % Number of blocks for averaging

% Initialize variable to accumulate cross-correlations:
avg_cross_corr = zeros(2*max_lag+1,1); 

% Obtain correlations for each block:
% Note the MATLAB syntax: xcov(x1,x2) is the expected value of x1(t+tau)*x2(t), so the first argument is the higher frequency band 

for i = 1:N_block

    start_idx = trunc + (i - 1) * step + 1;
    end_idx = start_idx + N - 1;
    
    % Check if the end index exceeds the length
    if end_idx > length(x_low_envelope)
        break; % Exit if there are not enough samples for another window
    end
    
    low_envelope_windowed = x_low_envelope(start_idx:end_idx);
    high_envelope_windowed = x_high_envelope(start_idx:end_idx);
    
    % Compute cross-correlation for the current window
    [cross_corr, ~] = xcov(high_envelope_windowed, low_envelope_windowed, max_lag, 'coeff');
    
    % Accumulate the cross-correlation results
    avg_cross_corr = avg_cross_corr + cross_corr; 

end

% Average the cross-correlation results
avg_cross_corr = avg_cross_corr / N_block;

% Identify the lag with the maximum correlation
[~, max_corr_idx] = max(avg_cross_corr);        % Index of the maximum correlation
max_corr_lag = max_corr_idx - (max_lag + 1);    % Convert index to sample lag
tau_max = max_corr_lag*dt;                      % Convert sample lag to time lag

% Step 6: Determine the direction of energy flow from the sign of tau
if max_corr_lag > 0
    direction = 'Energy flows from low frequency band to high frequency band';
elseif max_corr_lag < 0
    direction = 'Energy flows from high frequency band to low frequency band';
else
    direction = 'No clear directional energy flow';
end


% Display the result
if displ
    fprintf('Lag: %d samples\n', mean(max_corr_lag));
    fprintf('Time lag: %d s\n', mean(tau_max));
    fprintf('Direction: %s\n', direction);
end

% Plot the average cross-correlation result
if displ
    figure;
    box on;
    lags = [-max_lag:max_lag].*dt; % Time lag vector
    plot(1e6.* lags, avg_cross_corr, 'LineWidth', 1.5);
    title('Average Cross-correlation K(\tau): (' + low_band_str + ') \rightarrow (' + high_band_str + ')');
    xlabel(['\tau (' 956 's)']);
    ylabel('K(\tau)');
    set(gca, 'FontSize', 14);
    grid on;
end

