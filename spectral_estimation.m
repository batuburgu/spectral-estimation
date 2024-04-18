time_domain = linspace(-5, 5, 4001); % Time domain
fs = 1 / (time_domain(2) - time_domain(1)); % Sampling frequency

hamming_window_length = 4000; 
blackman_window_length = 4000;
MA_window_length = 4000;
ar_order = 1000;

sampled_sinusoidal_1 = sin(2 * pi * 100 * time_domain); % 100 Hz sinusoidal signal
sampled_sinusoidal_2 = sin(2 * pi * 150 * time_domain); % 150 Hz sinusoidal signal 

% Creation of the noisy signals
noisy_signal_1_10dB = awgn(sampled_sinusoidal_1, 10);
noisy_signal_1_0dB = awgn(sampled_sinusoidal_1, 0);
noisy_signal_1_m10dB = awgn(sampled_sinusoidal_1, -10);

noisy_signal_2_10dB = awgn(sampled_sinusoidal_2, 10);
noisy_signal_2_0dB = awgn(sampled_sinusoidal_2, 0);
noisy_signal_2_m10dB = awgn(sampled_sinusoidal_2, -10);

% Periodogram Calculations
[perio_freq_domain_1, periodogram_1_10dB] = periodogram_function(noisy_signal_1_10dB, fs);
[perio_freq_domain_2, periodogram_1_0dB] = periodogram_function(noisy_signal_1_0dB, fs);
[perio_freq_domain_3, periodogram_1_m10dB] = periodogram_function(noisy_signal_1_m10dB, fs);

[perio_freq_domain_4, periodogram_2_10dB] = periodogram_function(noisy_signal_2_10dB, fs);
[perio_freq_domain_5, periodogram_2_0dB] = periodogram_function(noisy_signal_2_0dB, fs);
[perio_freq_domain_6, periodogram_2_m10dB] = periodogram_function(noisy_signal_2_m10dB, fs);

% Plot for Signal 1 with SNR = +10 dB
figure(1);
plot(perio_freq_domain_1, periodogram_1_10dB,'Color','blue');
hold on
plot(-perio_freq_domain_1, periodogram_1_10dB,'Color','blue');
title('Periodogram Result for Signal 1 ( sin(2 * pi * 100 * t) ) with SNR = 10dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Signal 1 with SNR = 10dB')

% Plot for Signal 1 with SNR = 0 dB
figure(2);
plot(perio_freq_domain_2, periodogram_1_0dB,'Color','blue');
hold on
plot(-perio_freq_domain_2, periodogram_1_0dB,'Color','blue');
title('Periodogram Result for Signal 1 ( sin(2 * pi * 100 * t) ) with SNR = 0dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Signal 1 with SNR = 0dB')

% Plot for Signal 1 with SNR = -10 dB
figure(3);
plot(perio_freq_domain_3, periodogram_1_m10dB,'Color','blue');
hold on
plot(-perio_freq_domain_3, periodogram_1_m10dB,'Color','blue');
title('Periodogram Result for Signal 1 ( sin(2 * pi * 100 * t) ) with SNR = -10dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Signal 1 with SNR = -10dB')


% Plot for Signal 2 with SNR = +10 dB
figure(4);
plot(perio_freq_domain_4, periodogram_2_10dB,'Color','blue');
hold on
plot(-perio_freq_domain_4, periodogram_2_10dB,'Color','blue');
title('Periodogram Result for Signal 2 ( sin(2 * pi * 150 * t) ) with SNR = 10dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Signal 2 with SNR = 10dB')

% Plot for Signal 2 with SNR = 0 dB
figure(5);
plot(perio_freq_domain_5, periodogram_2_0dB,'Color','red');
hold on
plot(-perio_freq_domain_5, periodogram_2_0dB,'Color','red');
title('Periodogram Result for Signal 2 ( sin(2 * pi * 150 * t) ) with SNR = 0dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Signal 2 with SNR = 0dB')

% Plot for Signal 1 with SNR = -10 dB
figure(6);
plot(perio_freq_domain_6, periodogram_2_m10dB,'Color','black');
hold on
plot(-perio_freq_domain_6, periodogram_2_m10dB,'Color','black');
title('Periodogram Result for Signal 2 ( sin(2 * pi * 150 * t) ) with SNR = -10dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Signal 2 with SNR = -10dB')

%Welch's Method Calculations
[welchs_freq_domain_1, welchs_power_1_10dB] = welchs_function(noisy_signal_1_10dB, fs, hamming_window_length, hamming_window_length/2);
[welchs_freq_domain_2, welchs_power_1_0dB] = welchs_function(noisy_signal_1_0dB, fs, hamming_window_length, hamming_window_length/2);
[welchs_freq_domain_3, welchs_power_1_m10dB] = welchs_function(noisy_signal_1_m10dB, fs, hamming_window_length, hamming_window_length/2);

[welchs_freq_domain_4, welchs_power_2_10dB] = welchs_function(noisy_signal_2_10dB, fs, hamming_window_length, hamming_window_length/2);
[welchs_freq_domain_5, welchs_power_2_0dB] = welchs_function(noisy_signal_2_0dB, fs, hamming_window_length, hamming_window_length/2);
[welchs_freq_domain_6, welchs_power_2_m10dB] = welchs_function(noisy_signal_2_m10dB, fs, hamming_window_length, hamming_window_length/2);

% Plot for Signal 1 with SNR = +10 dB
figure(7);
plot(welchs_freq_domain_1, welchs_power_1_10dB,'Color','red');
hold on
plot(-welchs_freq_domain_1, welchs_power_1_10dB,'Color','red');
title('Welchs Method Result for Signal 1 ( sin(2 * pi * 100 * t) ) with SNR = 10dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Signal 1 with SNR = 10dB')

% Plot for Signal 1 with SNR = 0 dB
figure(8);
plot(welchs_freq_domain_2, welchs_power_1_0dB,'Color','red');
hold on
plot(-welchs_freq_domain_2, welchs_power_1_0dB,'Color','red');
title('Welchs Method Result for Signal 1 ( sin(2 * pi * 100 * t) ) with SNR = 0dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Signal 1 with SNR = 0dB')

% Plot for Signal 1 with SNR = -10 dB
figure(9);
plot(welchs_freq_domain_3, welchs_power_1_m10dB,'Color','red');
hold on
plot(-welchs_freq_domain_3, welchs_power_1_m10dB,'Color','red');
title('Welchs Method Result for Signal 1 ( sin(2 * pi * 100 * t) ) with SNR = -10dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Signal 1 with SNR = -10dB')

% Plot for Signal 2 with SNR = +10 dB
figure(10);
plot(welchs_freq_domain_4, welchs_power_2_10dB,'Color','red');
hold on
plot(-welchs_freq_domain_4, welchs_power_2_10dB,'Color','red');
title('Welchs Method Result for Signal 2 ( sin(2 * pi * 150 * t) ) with SNR = 10dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Signal 2 with SNR = 10dB')

% Plot for Signal 2 with SNR = 0 dB
figure(11);
plot(welchs_freq_domain_5, welchs_power_2_0dB,'Color','red');
hold on
plot(-welchs_freq_domain_5, welchs_power_2_0dB,'Color','red');
title('Welchs Method Result for Signal 2 ( sin(2 * pi * 150 * t) ) with SNR = 0dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Signal 2 with SNR = 0dB')

% Plot for Signal 2 with SNR = -10 dB
figure(12);
plot(welchs_freq_domain_6, welchs_power_2_m10dB,'Color','red');
hold on
plot(-welchs_freq_domain_6, welchs_power_2_m10dB,'Color','red');
title('Welchs Method Result for Signal 2 ( sin(2 * pi * 150 * t) ) with SNR = -10dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Signal 2 with SNR = -10dB')

Blackman-Tukey Calculations
[blackman_freq_domain_1, blackman_power_1_10dB] = blackman_tukey_function(noisy_signal_1_10dB, fs, blackman_window_length, blackman_window_length/2);
[blackman_freq_domain_2, blackman_power_1_0dB] = blackman_tukey_function(noisy_signal_1_0dB, fs, blackman_window_length, blackman_window_length/2);
[blackman_freq_domain_3, blackman_power_1_m10dB] = blackman_tukey_function(noisy_signal_1_m10dB, fs, blackman_window_length, blackman_window_length/2);

[blackman_freq_domain_4, blackman_power_2_10dB] = blackman_tukey_function(noisy_signal_2_10dB, fs, blackman_window_length, blackman_window_length/2);
[blackman_freq_domain_5, blackman_power_2_0dB] = blackman_tukey_function(noisy_signal_2_0dB, fs, blackman_window_length, blackman_window_length/2);
[blackman_freq_domain_6, blackman_power_2_m10dB] = blackman_tukey_function(noisy_signal_2_m10dB, fs, blackman_window_length, blackman_window_length/2);

% Plot for Signal 1 with SNR = +10 dB
figure(13);
plot(blackman_freq_domain_1, blackman_power_1_10dB,'Color','black');
hold on
plot(-blackman_freq_domain_1, blackman_power_1_10dB,'Color','black');
title('Blackman-Tukey Method Result for Signal 1 ( sin(2 * pi * 100 * t) ) with SNR = 10dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Signal 1 with SNR = 10dB')

% Plot for Signal 1 with SNR = 0 dB
figure(14);
plot(blackman_freq_domain_2, blackman_power_1_0dB,'Color','black');
hold on
plot(-blackman_freq_domain_2, blackman_power_1_0dB,'Color','black');
title('Blackman-Tukey Method Result for Signal 1 ( sin(2 * pi * 100 * t) ) with SNR = 0dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Signal 1 with SNR = 0dB')

% Plot for Signal 1 with SNR = -10 dB
figure(15);
plot(blackman_freq_domain_3, blackman_power_1_m10dB,'Color','black');
hold on
plot(-blackman_freq_domain_3, blackman_power_1_m10dB,'Color','black');
title('Blackman-Tukey Method Result for Signal 1 ( sin(2 * pi * 100 * t) ) with SNR = -10dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Signal 1 with SNR = -10dB')

% Plot for Signal 2 with SNR = 10 dB
figure(16);
plot(blackman_freq_domain_4, blackman_power_2_10dB,'Color','black');
hold on
plot(-blackman_freq_domain_4, blackman_power_2_10dB,'Color','black');
title('Blackman-Tukey Method Result for Signal 2 ( sin(2 * pi * 150 * t) ) with SNR = 10dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Signal 2 with SNR = 10dB')

% Plot for Signal 2 with SNR = 0 dB
figure(17);
plot(blackman_freq_domain_5, blackman_power_2_0dB,'Color','black');
hold on
plot(-blackman_freq_domain_5, blackman_power_2_0dB,'Color','black');
title('Blackman-Tukey Method Result for Signal 2 ( sin(2 * pi * 150 * t) ) with SNR = 0dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Signal 1 with SNR = 0dB')

% Plot for Signal 2 with SNR = -10 dB
figure(18);
plot(blackman_freq_domain_6, blackman_power_2_m10dB,'Color','black');
hold on
plot(-blackman_freq_domain_6, blackman_power_2_m10dB,'Color','black');
title('Blackman-Tukey Method Result for Signal 2 ( sin(2 * pi * 150 * t) ) with SNR = -10dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Signal 1 with SNR = -10dB')

AR (Autoregressive) Model Calculations
[ar_freq_domain_1, ar_power_1_10dB] = AR_function(noisy_signal_1_10dB, fs, ar_order);
[ar_freq_domain_2, ar_power_1_0dB] = AR_function(noisy_signal_1_0dB, fs, ar_order);
[ar_freq_domain_3, ar_power_1_m10dB] = AR_function(noisy_signal_1_m10dB, fs, ar_order);

[ar_freq_domain_4, ar_power_2_10dB] = AR_function(noisy_signal_2_10dB, fs, ar_order);
[ar_freq_domain_5, ar_power_2_0dB] = AR_function(noisy_signal_2_0dB, fs, ar_order);
[ar_freq_domain_6, ar_power_2_m10dB] = AR_function(noisy_signal_2_m10dB, fs, ar_order);

% Plot for Signal 1 with SNR = 10dB
figure(19);
plot(ar_freq_domain_1, ar_power_1_10dB);
hold on
plot(-ar_freq_domain_1, ar_power_1_10dB);
title('Autoregressive Method Result for Signal 1 ( sin(2 * pi * 100 * t) ) with SNR = 10dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)')
legend('Signal 1 with SNR = 10dB')

% Plot for Signal 1 with SNR = 0dB
figure(20);
plot(ar_freq_domain_2, ar_power_1_0dB, 'Color', 'magenta');
hold on
plot(-ar_freq_domain_2, ar_power_1_0dB, 'Color', 'magenta');
title('Autoregressive Method Result for Signal 1 ( sin(2 * pi * 100 * t) ) with SNR = 0dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)')
legend('Signal 1 with SNR = 0dB')

% Plot for Signal 1 with SNR = -10dB
figure(21);
plot(ar_freq_domain_3, ar_power_1_m10dB, 'Color', 'magenta');
hold on
plot(-ar_freq_domain_3, ar_power_1_m10dB, 'Color', 'magenta');
title('Autoregressive Method Result for Signal 1 ( sin(2 * pi * 100 * t) ) with SNR = -10dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)')
legend('Signal 1 with SNR = -10dB')

% Plot for Signal 2 with SNR = 10dB
figure(22);
plot(ar_freq_domain_4, ar_power_2_10dB, 'Color', 'magenta');
hold on
plot(-ar_freq_domain_4, ar_power_2_10dB, 'Color', 'magenta');
title('Autoregressive Method Result for Signal 2 ( sin(2 * pi * 150 * t) ) with SNR = 10dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)')
legend('Signal 2 with SNR = 10dB')

% Plot for Signal 2 with SNR = 0dB
figure(23);
plot(ar_freq_domain_5, ar_power_2_0dB, 'Color', 'magenta');
hold on
plot(-ar_freq_domain_5, ar_power_2_0dB, 'Color', 'magenta');
title('Autoregressive Method Result for Signal 2 ( sin(2 * pi * 150 * t) ) with SNR = 0dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)')
legend('Signal 2 with SNR = 0dB')

% Plot for Signal 2 with SNR = -10dB
figure(24);
plot(ar_freq_domain_6, ar_power_2_m10dB, 'Color', 'magenta');
hold on
plot(-ar_freq_domain_6, ar_power_2_m10dB, 'Color', 'magenta');
title('Autoregressive Method Result for Signal 2 ( sin(2 * pi * 150 * t) ) with SNR = -10dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)')
legend('Signal 2 with SNR = -10dB')

% MA (Moving Average) Model Calculations
[ma_freq_domain_1, ma_power_1_10dB] = MA_function(noisy_signal_1_10dB, fs, MA_window_length);
[ma_freq_domain_2, ma_power_1_0dB] = MA_function(noisy_signal_1_0dB, fs, MA_window_length);
[ma_freq_domain_3, ma_power_1_m10dB] = MA_function(noisy_signal_1_m10dB, fs, MA_window_length);

[ma_freq_domain_4, ma_power_2_10dB] = MA_function(noisy_signal_2_10dB, fs, MA_window_length);
[ma_freq_domain_5, ma_power_2_0dB] = MA_function(noisy_signal_2_0dB, fs, MA_window_length);
[ma_freq_domain_6, ma_power_2_m10dB] = MA_function(noisy_signal_2_m10dB, fs, MA_window_length);

% Plot for Signal 1 with SNR = 10dB
figure(25);
plot(ma_freq_domain_1, ma_power_1_10dB, 'Color', 'cyan');
hold on
plot(-ma_freq_domain_1, ma_power_1_10dB, 'Color', 'cyan');
title('Moving Average Method Result for Signal 1 ( sin(2 * pi * 100 * t) ) with SNR = 10dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)')
legend('Signal 1 with SNR = 10dB')

% Plot for Signal 1 with SNR = 0dB
figure(26);
plot(ma_freq_domain_2, ma_power_1_0dB, 'Color', 'cyan');
hold on
plot(-ma_freq_domain_2, ma_power_1_0dB, 'Color', 'cyan');
title('Moving Average Method Result for Signal 1 ( sin(2 * pi * 100 * t) ) with SNR = 0dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)')
legend('Signal 1 with SNR = 0dB')

% Plot for Signal 1 with SNR = -10dB
figure(27);
plot(ma_freq_domain_3, ma_power_1_m10dB, 'Color', 'cyan');
hold on
plot(-ma_freq_domain_3, ma_power_1_m10dB, 'Color', 'cyan');
title('Moving Average Method Result for Signal 1 ( sin(2 * pi * 100 * t) ) with SNR = -10dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)')
legend('Signal 1 with SNR = -10dB')

% Plot for Signal 2 with SNR = 10dB
figure(28);
plot(ma_freq_domain_4, ma_power_2_10dB, 'Color', 'cyan');
hold on
plot(-ma_freq_domain_4, ma_power_2_10dB, 'Color', 'cyan');
title('Moving Average Method Result for Signal 2 ( sin(2 * pi * 150 * t) ) with SNR = 10dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)')
legend('Signal 2 with SNR = 10dB')

% Plot for Signal 2 with SNR = 0dB
figure(29);
plot(ma_freq_domain_5, ma_power_2_0dB, 'Color', 'cyan');
hold on
plot(-ma_freq_domain_5, ma_power_2_0dB, 'Color', 'cyan');
title('Moving Average Method Result for Signal 2 ( sin(2 * pi * 150 * t) ) with SNR = 0dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)')
legend('Signal 2 with SNR = 0dB')

% Plot for Signal 2 with SNR = -10dB
figure(30);
plot(ma_freq_domain_6, ma_power_2_m10dB, 'Color', 'cyan');
hold on
plot(-ma_freq_domain_6, ma_power_2_m10dB, 'Color', 'cyan');
title('Moving Average Method Result for Signal 2 ( sin(2 * pi * 150 * t) ) with SNR = -10dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)')
legend('Signal 2 with SNR = -10dB')

function [freq_domain, MA_power_dB] = MA_function(signal, fs, window_length)
    % Length of the signal
    N = length(signal);

    % FFT length
    fft_length = 2^nextpow2(N);

    plain_PSD_estimation = abs(fft(signal, fft_length)).^2 / N;

    % MA filter
    MA_filter_freq_domain = rectwin(window_length);
    MA_filter_time_domain = ifft(MA_filter_freq_domain)*fs;

    % Filtering the signal with MA filter
    filtered_signal = conv(signal, MA_filter_time_domain, 'same');

    % PSD estimation for the filtered signal
    filtered_PSD_estimation = abs(fft(filtered_signal, fft_length) / fs).^2 / N;

    % Taking the average of the plain PSD and the filtered PSD
    PSD = (plain_PSD_estimation + filtered_PSD_estimation) / 2;

    freq_domain =  (0:(fft_length/2)-1) * (fs / (fft_length));
    
    % Normalized PSD estimation in decibels
    MA_power_dB = 10 * log10(PSD / max(PSD));
    MA_power_dB = MA_power_dB(1:length(MA_power_dB) / 2);
end

function [freq_domain, AR_power_dB] = AR_function(signal, fs, order)
    % Length of the signal
    N = length(signal);
    
    % Calculating the autocorrelation matrix
    temp=signal(1:order)' .* signal(1:order);
    expected_value=mean(temp);
    Rxx=toeplitz(expected_value);

    % Calculating the variance of the signal
    variance = var(signal(1:order));

    % Creating a vector with all zeros except the first element
    variance_matrix = zeros(order,1);
    % Assigning the variance of the signal to the first element 
    variance_matrix(1) = variance;
    % Calculating the AR_coefficients with: variance_matrix * inv(Rxx)
    AR_coefficients = variance_matrix\Rxx;
    
    freq_domain = (0:order-1) * (fs / order);

    % Finding the transfer function to obtain the PSD estimation
    transfer_function = zeros(size(freq_domain));
    for i = 1:length(freq_domain)
        exponential = exp(-1i * 2 * pi * freq_domain(i) / fs);
        transfer_function(i) = 1 / (1 - sum(AR_coefficients .* exponential.^(-1 * (1:order))));

    end
    % Estimating the PSD and taking the positive side
    PSD_estimation = variance ./ abs(transfer_function).^2;
    PSD_estimation = PSD_estimation(1 : length(PSD_estimation) / 2);
    
    AR_power_dB = 10 * log10(PSD_estimation / max(PSD_estimation));
    % Half of the frequency domain
    freq_domain = freq_domain(1:length(freq_domain) / 2);
end

function [freq_domain, blackman_power_dB] = blackman_tukey_function(signal,fs,window_length,nooverlap)
    % Length of the signal
    N = length(signal);
    
    % Calculating the segment length and the number of windows
    segment_length = window_length - nooverlap;
    number_of_windows = floor((N - nooverlap) / segment_length);

    freq_domain = (0:window_length/2 - 1) * fs / window_length;

    blackman_power = zeros(1, window_length);
    
    % Taking the autocorrelation of the signal
    signal_correlation = xcorr(signal);
    
    for i = 1:number_of_windows
        % The first and the last index of the corresponding signal segment
        start_index = (i - 1) * segment_length + 1;
        end_index = start_index + window_length - 1;
        
        % Autocorrelation of the corresponding signal segment
        segment_correlation = signal_correlation(start_index:end_index);
        
        % Windowing the autocorrelation
        windowed_correlation = segment_correlation .* blackman(window_length)';
        
        % Estimating PSD
        blackman_power = blackman_power + abs(fft(windowed_correlation)); 
    end
    blackman_power = blackman_power / number_of_windows;
    
    % Normalized decibel values of the PSD
    blackman_power_dB = 10 * log10(blackman_power / max(blackman_power));
    blackman_power_dB = blackman_power_dB(1:length(blackman_power_dB)/2);
end

function [freq_domain, welch_power_dB] = welchs_function(signal, fs, window_length, nooverlap)
    % Length of the signal
    N = length(signal);
    
    window = transpose(hamming(window_length));
    
    % Calculating the number of segments and starting points of them
    number_of_segments = floor((N - nooverlap) / (window_length - nooverlap));
    segment_starting_points = 1 + (0:number_of_segments - 1) * (window_length - nooverlap);
    
    freq_domain = (0:window_length/2-1) * fs / window_length;

    total_power = zeros(size(freq_domain));

    for i = 1:number_of_segments
        % Taking the corresponding segment
        segment = signal(segment_starting_points(i):segment_starting_points(i) + window_length - 1);
        % Applying windowing to the segment taken
        segment = segment .* window;
        
        % Taking the corresponding segment's fourier transform and
        % estimating the power spectrum
        segment_fourier = fft(segment);
        segment_power_spectrum = (1 / window_length) * abs(segment_fourier(1:window_length/2)) .^2;
    
        % Adding the segment's power spectrum to the power spectrum of the signal
        total_power = total_power + segment_power_spectrum;
        
    end
    
    % Normalized decibel values of the total power
    welch_power_dB = 10 * log10(total_power / max(total_power));    

end

function [freq_domain, perio_power_dB] = periodogram_function(signal, fs)
    % Length of the signal
    N = length(signal);
    
    % Frequencies of the corresponding signal
    freq_domain = (0:N/2-1) * fs / N;
    
    % Fourier transform of the signal
    signal_fourier = fft(signal);
    
    % Power spectrum of the signal
    power_spectrum = (1 / N) * abs(signal_fourier(1:N/2)).^2;
    
    % Normalized decibel values of the power spectrum 
    perio_power_dB = 10 * log10(power_spectrum / max(power_spectrum));
end

