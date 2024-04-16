time_domain = linspace(-5, 5, 4000); % Time domain
fs = 1 / (time_domain(2) - time_domain(1)); % Sampling frequency

hamming_window_length = 4000; 

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
plot(perio_freq_domain_2, periodogram_1_0dB,'Color','red');
hold on
plot(-perio_freq_domain_2, periodogram_1_0dB,'Color','red');
title('Periodogram Result for Signal 1 ( sin(2 * pi * 100 * t) ) with SNR = 0dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Signal 1 with SNR = 0dB')

% Plot for Signal 1 with SNR = -10 dB
figure(3);
plot(perio_freq_domain_3, periodogram_1_m10dB,'Color','black');
hold on
plot(-perio_freq_domain_3, periodogram_1_m10dB,'Color','black');
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
plot(welchs_freq_domain_1, welchs_power_1_10dB,'Color','blue');
hold on
plot(-welchs_freq_domain_1, welchs_power_1_10dB,'Color','blue');
title('Welchs Method Result for Signal 1 ( sin(2 * pi * 100 * t) ) with SNR = 10dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Signal 1 with SNR = 10dB')

% Plot for Signal 1 with SNR = 0 dB
figure(8);
plot(welchs_freq_domain_2, welchs_power_1_0dB,'Color','blue');
hold on
plot(-welchs_freq_domain_2, welchs_power_1_0dB,'Color','blue');
title('Welchs Method Result for Signal 1 ( sin(2 * pi * 100 * t) ) with SNR = 0dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Signal 1 with SNR = 0dB')

% Plot for Signal 1 with SNR = -10 dB
figure(9);
plot(welchs_freq_domain_3, welchs_power_1_m10dB,'Color','blue');
hold on
plot(-welchs_freq_domain_3, welchs_power_1_m10dB,'Color','blue');
title('Welchs Method Result for Signal 1 ( sin(2 * pi * 100 * t) ) with SNR = -10dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Signal 1 with SNR = -10dB')

% Plot for Signal 2 with SNR = +10 dB
figure(10);
plot(welchs_freq_domain_4, welchs_power_2_10dB,'Color','blue');
hold on
plot(-welchs_freq_domain_4, welchs_power_2_10dB,'Color','blue');
title('Welchs Method Result for Signal 2 ( sin(2 * pi * 150 * t) ) with SNR = 10dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Signal 2 with SNR = 10dB')

% Plot for Signal 2 with SNR = 0 dB
figure(11);
plot(welchs_freq_domain_5, welchs_power_2_0dB,'Color','blue');
hold on
plot(-welchs_freq_domain_5, welchs_power_2_0dB,'Color','blue');
title('Welchs Method Result for Signal 2 ( sin(2 * pi * 150 * t) ) with SNR = 0dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Signal 2 with SNR = 0dB')

% Plot for Signal 2 with SNR = -10 dB
figure(12);
plot(welchs_freq_domain_6, welchs_power_2_m10dB,'Color','blue');
hold on
plot(-welchs_freq_domain_6, welchs_power_2_m10dB,'Color','blue');
title('Welchs Method Result for Signal 2 ( sin(2 * pi * 150 * t) ) with SNR = -10dB');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Signal 2 with SNR = -10dB')


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

