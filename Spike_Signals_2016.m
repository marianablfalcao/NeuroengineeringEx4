clear all
close all
neuralData=load ('Proprio_1.txt');
Fs=24414.06;

% Compute the time vector from the frequendy value
time=;

%% Pre-processing of the signal 
figure(1)
plot(time,neuralData')
title('Original Signal'); 

%low pass 
Fb=5000;
%High Pass
Fh=500;

neuralData = filtra(neuralData,Fs,Fb,Fh);


%% Spike Extraction 

%%% Vary the threshold from 2 to 6  of the spike extracting algorithm
%%% (spike_extract.m)  and observe how this affects the detection of spike
%%% events and the computation time 

timeWindow = 0.002;
[filteredSpikes, spikesIndex]=spike_extract(neuralData,threshold,Fs,timeWindow);

% Plot the waveform of the extracted spikes (filteredSpikes) for each
% threshold and choose the optimal threshold value



%% Principal component analysis

% Extract the filteredSpikes using the optimal threshold (found previously)

threshold= ;
[filteredSpikes, spikesIndex]=spike_extract_solution(neuralData,threshold,Fs,timeWindow);

% Compute the principal components of the filteredSpikes (size n= number of spikes)
% Hint: use Matlab's function pca


% Store the representation of the filteredSpikes on the 2 first PC -  in a nx2 matrix 



%% Spike sorting

% Cluster the spikes into 2 groups using the kmeans algorithm (already
% implemented in matlab) and store the cluster indices in a vector "idx".


% Plot the two clusters (two distinct colours) of data and the centroids of each cluster 


% Plot the  waveforms contained in the two clusters in two distinct colours 




%% Calculation of the Firing Rate 
% Compute the firing rate using a time window of 200 ms (and the optimal
% threshold found previously)

% Plot the original neural signal overlapped with the detected spikes

% Plot the original neural signal overlapped with the calculated firing
% rate