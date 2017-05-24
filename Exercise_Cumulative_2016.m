    %% Working with the ENG rat data recorded from the experiments with cuff electrodes

    clear all
    close all
    clear
    clc

    %% Data Import 

    cd '/Users/marianafalcao/Desktop/erasmus cadeiras/Fundamentals of Neuroengineering/Exercices/Exercise 4/Ex1'

    d=importdata('Brush.txt');        %nocio



    %% Raw data and labels
    signal=d(:,1);
    labels=d(:,2);
    labels=round(labels); % Because there could be some noise in the analogical signal. 
    Fs = 20000;
    time=(1:size(signal,1))/Fs;

    figure(1)
    subplot(3,1,3)
    plot(time,zscore(signal));
    hold on;
    plot(time,zscore(labels),'r','LineWidth',1.5);
    title('ENG original signal and stimuli application labels - Nociceptive stimulus')
    xlim([0,time(1,end)])
    ylim([-6,6])


    %% Examination of the spectral context of the signal-PSD

    [rows_act,cols_act,values_act] = find(labels>0);
    [rows_rest,cols_rest,values_rest] = find(labels==0);

    figure(2)
    subplot(1,3,3)
    signalOfInterest=signal(rows_act);
    notOfInterest=signal(rows_rest);
    h = spectrum.welch; % creates the Welch spectrum estimator
    SOIf=psd(h,signalOfInterest,'Fs',Fs); % calculates and plot the one sided PSD
    NOIf=psd(h,notOfInterest,'Fs',Fs);
    plot(SOIf); % Plot the one-sided PSD. 
    hold on;
    plot(NOIf);
    title('Spectrum of SOI, and NOI for the stimuli application - Nociceptive stimulus')
    legend('stimuli activity','rest activity')

    %% Filtering 

    %Determine the Low pass cutoff frequency for your filter
    Fl=800;
    %Determine the High pass cutoff frequency for your filter
    Fh=2200;

    [Signal_filtered]=filtra(signal',Fs,Fl,Fh);

    figure(3)
    subplot(1,3,3)
    soi_f=Signal_filtered(rows_act);
    noi_f=Signal_filtered(rows_rest);
    h1 = spectrum.welch; % creates the Welch spectrum estimator
    SOI_f=psd(h,soi_f,'Fs',Fs); % calculates and plot the one sided PSD
    NOI_f=psd(h,noi_f,'Fs',Fs);
    plot(SOI_f); % Plot the one-sided PSD. 
    hold on;
    plot(NOI_f);
    title('Spectrum of SOI, and NOI post-filtering - Nociceptive stimulus')
    legend('stimuli activity','rest activity')
    
    %% Feature extraction 

    % Determine the step size                                          
    step=Fs * 100e-3; 

    % Compute the MAV and Zero-Crossing features for each time window                      
    % Hint: MAV and ZeroCross can be functions in separate files         

    pos=1;
    L2=length(signal);
    for i=1:step:(L2-step)
        window=signal(i:i+step-1,1);
        vec_mav(pos,1)=mav(window);
        vec_zCross(pos,1)=zeroCross(window);
        vec_wave(pos,1)=waveLen(window);
        pos=pos+1;
    end

    % Resizing the label vector to the size of the features vector
    for i=1:step:(L2-step);  
        z=ceil(i/step);
        labels_resized(z)=labels(i);
    end

    % plot the zscored features along with the corresponding labels
    figure(4)
    subplot(3,1,1)
    plot(zscore(vec_mav));
    hold on;
    plot(zscore(labels_resized),'r');
    xlim([0,size(labels_resized,2)])
    title('MAV features - Nociceptive stimulus')
    % title('ENG reduced signal and stimuli application labels')
    subplot(3,1,2)
    plot(zscore(vec_zCross));
    hold on;
    plot(zscore(labels_resized),'r');
    xlim([0,size(labels_resized,2)])
    title('Zero Cross features - Nociceptive stimulus')
    subplot(3,1,3)
    plot(zscore(vec_wave));
    hold on;
    plot(zscore(labels_resized),'r');
    xlim([0,size(labels_resized,2)])
    title('Wavelength features - Nociceptive stimulus')


%% Signal-to-Noise Estimation

% Extract the features for the non zero labels
[rows_act,cols_act,values_act] = find(labels_resized'>0);
non_zero_mav=vec_mav(rows_act);
non_zero_zcross=vec_zCross(rows_act);
non_zero_wl=vec_wave(rows_act);

% Compute Activation ??


% Extract the features for the labels = zero
[rows_rest,cols_rest,values_rest] = find(labels_resized'==0);
zero_mav=vec_mav(rows_rest);
zero_zcross=vec_zCross(rows_rest);
zero_wl=vec_wave(rows_rest);

% Compute the noise, SNR and SNR in dB
snr_mav1=mean(non_zero_mav,1)/mean(zero_mav,1);
snr_mav1_dB=20 * log10(snr_mav1);
snr_zcross1=mean(non_zero_zcross,1)/mean(zero_zcross,1);
snr_zcross1_dB=20 * log10(snr_zcross1);
snr_wl1=mean(non_zero_wl,1)/mean(zero_wl,1);
snr_wl1_dB=20 * log10(snr_wl1);


%% %% %% %%
%% perform the same for the other data sets:

c=importdata('Touch.txt');          %mechanical stimulus - VF
b=importdata('FootFlexion.txt');    %prioperceptive
Fs=20000;
Fl=800;
Fh=2200;

%% Raw data and labels (same Fs)
signal_2=c(:,1);
labels_2=c(:,2);
labels_2=round(labels_2);
time_2=(size(signal,1):size(signal,1)+size(signal_2,1)-1)/Fs;
signal_3=b(:,1);
labels_3=b(:,2);
labels_3=round(labels_3);
time_3=(size(signal,1)+size(signal_2,1):size(signal,1)+size(signal_2,1)+size(signal_3,1)-1)/Fs;

    figure(1)
    subplot(3,1,1)
    plot(time_2,zscore(signal_2));
    hold on;
    plot(time_2,zscore(labels_2),'y','LineWidth',1.5);
    title('ENG original signal and stimuli application labels - Mechanical Stimulus')
    xlim([time(1,end),time_2(1,end)])
    ylim([-6,6])
    subplot(3,1,2)
    plot(time_3,zscore(signal_3));
    hold on;
    plot(time_3,zscore(labels_3),'g','LineWidth',1.5);
    title('ENG original signal and stimuli application labels - Propioceptive Stimulus')
    xlim([time_2(1,end),time_3(1,end)])
    ylim([-6,6])
    
%% Examination of the spectral context of the signal-PSD
[rows_act_2,cols_act,values_act] = find(labels_2>0);
[rows_rest_2,cols_rest,values_rest] = find(labels_2==0);
signalOfInterest_2=signal_2(rows_act_2);
notOfInterest_2=signal_2(rows_rest_2);

[rows_act_3,cols_act,values_act] = find(labels_3>0);
[rows_rest_3,cols_rest,values_rest] = find(labels_3==0);
signalOfInterest_3=signal_3(rows_act_3);
notOfInterest_3=signal_3(rows_rest_3);

    figure(2)
    subplot(1,3,1)
    h = spectrum.welch; % creates the Welch spectrum estimator
    SOIf=psd(h,signalOfInterest_2,'Fs',Fs); % calculates and plot the one sided PSD
    NOIf=psd(h,notOfInterest_2,'Fs',Fs);
    plot(SOIf); % Plot the one-sided PSD. 
    hold on;
    plot(NOIf);
    title('Spectrum of SOI, and NOI for the stimuli application - Mechanical Stimulus')
    legend('stimuli activity','rest activity')
    subplot(1,3,2)
    h = spectrum.welch; % creates the Welch spectrum estimator
    SOIf=psd(h,signalOfInterest_3,'Fs',Fs); % calculates and plot the one sided PSD
    NOIf=psd(h,notOfInterest_3,'Fs',Fs);
    plot(SOIf); % Plot the one-sided PSD. 
    hold on;
    plot(NOIf);
    title('Spectrum of SOI, and NOI for the stimuli application - Propioceptive Stimulus')
    legend('stimuli activity','rest activity')

%% Filtering (same cutoff frequencies)
[Signal_filtered_2]=filtra(signal_2',Fs,Fl,Fh);
[Signal_filtered_3]=filtra(signal_3',Fs,Fl,Fh);

    figure(3)
    subplot(1,3,1)
    soi_f=Signal_filtered_2(rows_act_2);
    noi_f=Signal_filtered_2(rows_rest_2);
    h1 = spectrum.welch; % creates the Welch spectrum estimator
    SOI_f=psd(h,soi_f,'Fs',Fs); % calculates and plot the one sided PSD
    NOI_f=psd(h,noi_f,'Fs',Fs);
    plot(SOI_f); % Plot the one-sided PSD. 
    hold on;
    plot(NOI_f);
    title('Spectrum of SOI, and NOI post-filtering - Mechanical Stimulus')
    legend('stimuli activity','rest activity')
    subplot(1,3,2)
    soi_f=Signal_filtered_3(rows_act_3);
    noi_f=Signal_filtered_3(rows_rest_3);
    h1 = spectrum.welch; % creates the Welch spectrum estimator
    SOI_f=psd(h,soi_f,'Fs',Fs); % calculates and plot the one sided PSD
    NOI_f=psd(h,noi_f,'Fs',Fs);
    plot(SOI_f); % Plot the one-sided PSD. 
    hold on;
    plot(NOI_f);
    title('Spectrum of SOI, and NOI post-filtering - Proprioceptive Stimulus')
    legend('stimuli activity','rest activity')

%% Feature extraction (same step size)
step=Fs * 100e-3; 

pos=1;
L2=length(signal_2);
for i=1:step:(L2-step)
    window_2=signal_2(i:i+step-1,1);
    vec_mav_2(pos,1)=mav(window_2);
    vec_zCross_2(pos,1)=zeroCross(window_2);
    vec_wave_2(pos,1)=waveLen(window_2);
    pos=pos+1;
end
for i=1:step:(L2-step);  
    z=ceil(i/step);
    labels_resized_2(z)=labels_2(i);
end
figure(5)
subplot(3,1,1)
plot(zscore(vec_mav_2));
hold on;
plot(zscore(labels_resized_2),'r');
xlim([0,size(labels_resized_2,2)])
title('MAV features - Mechanical Stimulus')
% title('For VF situation')
subplot(3,1,2)
plot(zscore(vec_zCross_2));
hold on;
plot(zscore(labels_resized_2),'r');
xlim([0,size(labels_resized_2,2)])
title('Zero Cross features - Mechanical Stimulus')
subplot(3,1,3)
plot(zscore(vec_wave_2));
hold on;
plot(zscore(labels_resized_2),'r');
xlim([0,size(labels_resized_2,2)])
title('Wavelength features - Mechanical Stimulus')


pos=1;
L2=length(signal_3);
for i=1:step:(L2-step)
    window_3=signal_3(i:i+step-1,1);
    vec_mav_3(pos,1)=mav(window_3);
    vec_zCross_3(pos,1)=zeroCross(window_3);
    vec_wave_3(pos,1)=waveLen(window_3);
    pos=pos+1;
end
for i=1:step:(L2-step);  
    z=ceil(i/step);
    labels_resized_3(z)=labels_3(i);
end
figure(6)
subplot(3,1,1)
plot(zscore(vec_mav_3));
hold on;
plot(zscore(labels_resized_3),'r');
xlim([0,size(labels_resized_3,2)])
title('MAV features - Proprioceptive Stimulus')
subplot(3,1,2)
plot(zscore(vec_zCross_3));
hold on;
plot(zscore(labels_resized_3),'r');
xlim([0,size(labels_resized_3,2)])
title('Zero Cross features - Proprioceptive Stimulus')
subplot(3,1,3)
plot(zscore(vec_wave_3));
hold on;
plot(zscore(labels_resized_3),'r');
xlim([0,size(labels_resized_3,2)])
title('Wavelength features - Proprioceptive Stimulus')

%% Signal-to-Noise Estimation
[rows_act_2,cols_act,values_act] = find(labels_resized_2'>0);
non_zero_mav_2=vec_mav_2(rows_act_2);
non_zero_zcross2=vec_zCross_2(rows_act_2);
non_zero_wl2=vec_wave_2(rows_act_2);
[rows_rest_2,cols_rest,values_rest] = find(labels_resized_2'==0);
zero_mav_2=vec_mav_2(rows_rest_2);
zero_zcross_2=vec_zCross_2(rows_rest_2);
zero_wl_2=vec_wave_2(rows_rest_2);
snr_mav2=mean(non_zero_mav_2,1)/mean(zero_mav_2,1);
snr_mav2_dB=20 * log10(snr_mav2);
snr_zcross2=mean(non_zero_zcross2,1)/mean(zero_zcross_2,1);
snr_zcross2_dB=20 * log10(snr_zcross2);
snr_wl2=mean(non_zero_wl2,1)/mean(zero_wl_2,1);
snr_wl2_dB=20 * log10(snr_wl2);

[rows_act_3,cols_act,values_act] = find(labels_resized_3'>0);
non_zero_mav_3=vec_mav_3(rows_act_3);
non_zero_zcross3=vec_zCross_3(rows_act_3);
non_zero_wl3=vec_wave_3(rows_act_3);
[rows_rest_3,cols_rest,values_rest] = find(labels_resized_3'==0);
zero_mav_3=vec_mav_3(rows_rest_3);
zero_zcross_3=vec_zCross_3(rows_rest_3);
zero_wl_3=vec_wave_3(rows_rest_3);
snr_mav3=mean(non_zero_mav_3,1)/mean(zero_mav_3,1);
snr_mav3_dB=20 * log10(snr_mav3);
snr_zcross3=mean(non_zero_zcross3,1)/mean(zero_zcross_3,1);
snr_zcross3_dB=20 * log10(snr_zcross3);
snr_wl3=mean(non_zero_wl3,1)/mean(zero_wl_3,1);
snr_wl3_dB=20 * log10(snr_wl3);


