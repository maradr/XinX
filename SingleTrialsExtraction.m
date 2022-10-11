clear; close all;
totcond = 4;
OutPath = fullfile(pwd, 'Outputs');


counter = 1;
%% Extract & Process trials per subject
for cnd = 1:totcond % For each condition
    path = fullfile(pwd, ['Cnd', num2str(cnd)]); %set the right path
    ff = dir(fullfile(path,'*.set')); %read all .set files
    
    BC_all    = cell(1,1,length(ff)); % preallocate 
    
    for s = 1:length(ff)
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab; %open EEGlab
        EEG = pop_loadset('filename',ff(s).name,'filepath',path); %read each sbj file
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        
        Fs = EEG.srate; % Get Sampling rate
        L  = size(EEG.data,2); % Length of the signal (one trial)
        electrodes = size(EEG.data,1); % get number of electrodes (128)
        n_trials = size(EEG.data, 3); % get retained number of trials
        
        for trial = 1:n_trials
            disp([num2str(s) ', ' num2str(trial)]); %let me know you're not crashing stuck somewhere
            temp = EEG.data(:,:,trial); %temporary file with that specific trial
            D=temp'; % transpose it
            
            %run the freqanalysis magic:
            Y  = fft(D); % Fast Fourier Transform            
            Y  = sqrt(real(Y).^2+imag(Y).^2)/L;% square root of the sum of squares of the real and imaginary parts divided by the number of data points            
            Y  = Y(1:L/2+1,:); % Take half
            f  = Fs*(0:(L/2))/L; % Vector of values for the freq bins
            BC     = zeros(size(Y)); %allocate empty matrix to speed up computations
            
            for e = 1:electrodes % For each electrode
                for n = 12:length(Y)-12
                Ys = Y([n-11:n-2,n+2:n+11],e); % surrounding bins
                                               % excluding adjacent and two
                                               % extremes (min and max)
                [~, idmin] = min(Ys);
                [~, idmax] = max(Ys);
                Ys([idmin,idmax]) = [];
                BC(n,e)     = Y(n,e) - mean(Ys); % Baseline corrected
                end
            end      
            
            %store it!
            Data(counter).Sbj = ff(s).name;
            Data(counter).Cond = str2double(regexp(ff(s).name,'^[0-9]*','match'));
            Data(counter).Trial = trial;
            Data(counter).Y     = (Y);
            %Data(trial).zScore_all = cell2mat(zScore_all);
            Data(counter).BC    = BC;
            
            save(fullfile(OutPath,'Data.mat'),'Data')
            counter = counter + 1;
        end
    end
end


save('Data_all.mat', 'Data', '-v7.3')
% 580 raws: 29 sbjs x 4 conditions x 5 trials

% data all contains also the raw signal; BC has only the BC amplitude

%% Extract BC info (L & R, harmonics & their average)
CarrierFreq = 6;
CalcActualF
load Chanlocs
Lr = [9,10,11,127,15];
Rr = [38,39,40,43,28];

fi     = Ft;
fI     = Ft*(1:4)';
[~,fw_har] = min(abs(repmat(f,[length(fI),1]) - repmat(fI,[1,length(f)])),[],2);
fI_1     = Ft*1';
[~,fw_1] = min(abs(repmat(f,[length(fI_1),1]) - repmat(fI_1,[1,length(f)])),[],2);
fI_2     = Ft*2';
[~,fw_2] = min(abs(repmat(f,[length(fI_2),1]) - repmat(fI_2,[1,length(f)])),[],2);
fI_3     = Ft*3';
[~,fw_3] = min(abs(repmat(f,[length(fI_3),1]) - repmat(fI_3,[1,length(f)])),[],2);
fI_4     = Ft*4';
[~,fw_4] = min(abs(repmat(f,[length(fI_4),1]) - repmat(fI_4,[1,length(f)])),[],2);


for i=1:length(tempBC)
    tempBC(i).Har_Left = squeeze(sum(mean(tempBC(i).BC(fw_har,Lr,:),2)));
    tempBC(i).first_Left = squeeze(mean(tempBC(i).BC(fw_1,Lr,:),2));
    tempBC(i).second_Left = squeeze(mean(tempBC(i).BC(fw_2,Lr,:),2));
    tempBC(i).third_Left = squeeze(mean(tempBC(i).BC(fw_3,Lr,:),2));
    tempBC(i).fourth_Left = squeeze(mean(tempBC(i).BC(fw_4,Lr,:),2));
    tempBC(i).Har_Right = squeeze(sum(mean(tempBC(i).BC(fw_har,Rr,:),2)));
    tempBC(i).first_Right = squeeze(mean(tempBC(i).BC(fw_1,Rr,:),2));
    tempBC(i).second_Right = squeeze(mean(tempBC(i).BC(fw_2,Rr,:),2));
    tempBC(i).third_Right = squeeze(mean(tempBC(i).BC(fw_3,Rr,:),2));
    tempBC(i).fourth_Right = squeeze(mean(tempBC(i).BC(fw_4,Rr,:),2));
end

% Save only the necessary info to run analyses in R
Data_Short = rmfield(tempBC,'BC');
writetable(struct2table(Data_Short), 'SL30_singletrials.txt')


%%
summary = readtable('SL30_singletrials.txt');
tidy = table2struct(summary);


%%

