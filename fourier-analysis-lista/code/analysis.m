%Author: Larry Parker
%Project: Fourier analysis of lista wind power data


%% Step 1: read in the data and create a time series object

% make sure to be in the right path
filepath = strcat(pwd,'/data/Lista.csv')

data_table = readtable(filepath, 'HeaderLines', 0, 'Delimiter', ',');
data_table.POWER = str2double(data_table.POWER);, 

time_table = table2timetable(data_table);


time_range = timerange('2013-01-01', '2019-01-01');
time_series = time_table(time_range,:);

N = length(time_series.POWER);


%% Step 2: Define gausian filter kernel and plot results
% normalization to remove the DC offsetq
time_series.POWER = time_series.POWER - mean(time_series.POWER); % get rid of dc offset


hz = linspace(0,1,N);

fft_wind_POWER = fft(time_series.POWER);

pow = 2*(abs(fft_wind_POWER/N));

[MaxPOWER, Index] = maxk(pow,1)

freq_max = hz(Index) % frequency with highest amplitude / strongest cycle

% Create in frequency domain with a gaussian kernel

fwhm = 0.05; % modifying this parameter determines the width of the kernel

s = fwhm*(2*pi-1)/(4*pi);
x = hz-freq_max;
fx = exp(-.5*(x/s).^2); % define gaussain kernel
filtered_fft = fft_wind_POWER.*fx';
filtered_pow = 2*(abs(filtered_fft/N)); % filtered signal in frequency domain
filtered_signal = 2*real(ifft(filtered_fft)); % filgered signal in time domain


%% Step 3: Produce plots
clf, figure(1)
plot(time_series.TIMESTAMPS, time_series.POWER,'b', 'linew',1)
xlabel(''), ylabel('Wind power (MW)')
title('Time domain of original signal')
saveas(gcf,'graphs/Time_domain_original_data.png')

figure(2)
stem(hz, pow, 'b','ks-','markerfacecolor', 'w', 'linew',3, 'markersize', 10)
xlabel('Frequency (norm.)'), ylabel('Wind power (MW)')
%set(gca,'xlim',[max(0, (Index- 30)*(1/N)) (Index+ 30)*(1/N)])
title('Frequency domain of original signal')
saveas(gcf,'graphs/Frequency_domain_original_data.png')

figure(3)
plot(time_series.TIMESTAMPS,filtered_signal,'r','markerfacecolor', 'w','linew',1)
xlabel(''), ylabel('Wind power (MW)')
title('Time domain of filtered signal')
saveas(gcf,'graphs/Time_domain_filtered_data.png')

figure(4)
stem(hz, filtered_pow, 'r','ms-','markerfacecolor', 'w', 'linew',3, 'markersize', 10)
xlabel('Frequency (norm.)'), ylabel('Wind power (MW)')
%set(gca,'xlim',[max(0, (Index- 30)*(1/N)) (Index+ 30)*(1/N)])
title('Frequency domain of filtered signal')
saveas(gcf,'graphs/Frequency_domain_filtered_data.png')

%% Step 5: compare

figure(5)
plot(time_series.TIMESTAMPS, time_series.POWER,'b', 'linew',2), hold on
plot(time_series.TIMESTAMPS, filtered_signal,'r', 'linew',2)
title('Original data (blue) vs filered data (red)')
ylabel('Wind power (MW)')
saveas(gcf,'graphs/Original_data_VS_filtered_data.png')
