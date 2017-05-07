clear all
close all
neuralData=load ('Proprio_1.txt');
Fs=24414.06;

% Compute the time vector from the frequendy value
time=(1:length(neuralData))./Fs;

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
spiks ={};
ns = [];
for threshold = 2:6
    [filteredSpikes, spikesIndex]=spike_extract(neuralData,threshold,Fs,timeWindow);
    spiks{threshold-1} = filteredSpikes;
    ns(threshold-1) = (length(filteredSpikes));
end

% Plot the waveform of the extracted spikes (filteredSpikes) for each
% threshold and choose the optimal threshold value

figure(1)
plot(ns)

figure(2)
for j = 1:5
    subplot(2,3,j)
for i = 1:size(spiks{j},1)
    
    plot(spiks{j}(i,:))
    hold on
end
end


%% Principal component analysis

% Extract the filteredSpikes using the optimal threshold (found previously)

threshold= 3;
[filteredSpikes, spikesIndex]=spike_extract(neuralData,threshold,Fs,timeWindow);

% Compute the principal components of the filteredSpikes (size n= number of spikes)
% Hint: use Matlab's function pca

[coeff,score] = pca(filteredSpikes);


% Store the representation of the filteredSpikes on the 2 first PC -  in a nx2 matrix 
X = score(:,1:2);


%% Spike sorting

% Cluster the spikes into 2 groups using the kmeans algorithm (already
% implemented in matlab) and store the cluster indices in a vector "idx".
[idx, ctrs] = kmeans(X, 2);

% Plot the two clusters (two distinct colours) of data and the centroids of each cluster 
figure(3)
plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',12)
hold on
plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',12)

% Plot the  waveforms contained in the two clusters in two distinct colours 
figure(4)
hold on
for i = idx==2
plot(filteredSpikes(i,:)', 'b')
end
for i = idx==1
plot(filteredSpikes(i,:)', 'r')
end

figure(5)
hold on
a = mean(filteredSpikes(idx==1,:));
b = mean(filteredSpikes(idx==2,:));
for i = idx==2
plot(a', 'b')
end
for i = idx==1
plot(b', 'r')
end




%% Calculation of the Firing Rate 
% Compute the firing rate using a time window of 200 ms (and the optimal
% threshold found previously)
window = 0.2*Fs;
rateData = zeros(1,length(neuralData));
for i = 1:length(spikesIndex)
    start = spikesIndex(i);
    for j = start:start+window
        if j <= length(neuralData)
            rateData(j) = rateData(j) + 5;
        end
    end
end

%% h
% Plot the original neural signal overlapped with the detected spikes
figure(6)
plot(neuralData)
hold on
scatter(spikesIndex, 10^-5*ones(1,length(spikesIndex)), 5) 
% for i = 1:length(spikesIndex)
% scatter(spikesIndex(i):spikesIndex(i)+48, filteredSpikes(i,:)*0.2*10^-5,5) 
% end

% Plot the original neural signal overlapped with the calculated firing
% rate
figure(7)
yyaxis left
hold on
scatter(spikesIndex, 10^-5*ones(1,length(spikesIndex)), 5,'m') 
plot(neuralData)

yyaxis right
ylim([-120 70])
plot(rateData)