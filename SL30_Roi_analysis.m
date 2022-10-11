%% X in x - analysis on Regions of Interest
% L2R lab, Mara De Rosa
% Requires: Effect Size Toolbox (Harald Hentschke & Maik Stüttgen, 2015)
%           CalcActualF - from Yamil Vidal
% Condition labels, argument for "Cond"

clear all; close all; clc
Cond        = 'Words'; %Pseudofonts, Nonwords, Pseudowords, Words
CarrierFreq = 6;
Type        = 'T';  % NT (base stimulation, NonTarget)or T (oddball stimulation, Target)
ThePath = fullfile('Data'); %path to datasets

CalcActualF
load Chanlocs
load(fullfile(pwd,ThePath,'Freq',Cond))

%% Electrodes of interest
occROI = [9,10,11,12,...
    13,14,15,16,...
    22,23,24,25,...
    26,27,28,29,...
    40,41,39,38]; % Occipital ROI
if (Type == "NT")
    fOfInterest     = Fnt*(1:4)';
    [~,fw] = min(abs(repmat(f,[length(fOfInterest),1]) - repmat(fOfInterest,[1,length(f)])),[],2);
    
    tokeep = {'BC_all', 'Cond', 'Type', 'Ft', 'Fnt', 'fw', 'occROI'};
    removed_Variables = setdiff(who, tokeep);
    clear(removed_Variables{:})
    
    dataOcc = squeeze(sum(mean(BC_all(fw,occROI,:),2)));
    [~,Pocc,CIocc,STATSocc] = ttest(dataOcc);
    ESocc = mes(dataOcc,0,'g1');
    effectSize = num2str(ESocc.g1);
    display(['Base respons, Occipital ROI for ', Cond, ': t = ',num2str(STATSocc.tstat),', p = ',num2str(Pocc), ', effect size: = ',num2str(ESocc.g1),', CI: ',num2str(ESocc.g1Ci(1)),' to ',num2str(ESocc.g1Ci(2))])
    
elseif (Type == "T")
    fOfInterest     = Ft*(1:4)';
    [~,fw] = min(abs(repmat(f,[length(fOfInterest),1]) - repmat(fOfInterest,[1,length(f)])),[],2);
    leftROI = [9,10,11,127,15];
    rightROI = [38,39,40,43,28];
    tokeep = {'BC_all', 'Cond', 'Type', 'Ft', 'Fnt', 'fw', 'leftROI', 'rightROI'};
    removed_Variables = setdiff(who, tokeep);
    clear(removed_Variables{:})
    
    dataLeft = squeeze(sum(mean(BC_all(fw,leftROI,:),2)));
    [~,Pleft,CIleft,STATSleft] = ttest(dataLeft);
    ESleft = mes(dataLeft,0,'g1');
    effectSize = num2str(ESleft.g1);
    display(['Odd response, left ROI for: ', Cond,': t = ',num2str(STATSleft.tstat),', p = ',num2str(Pleft), ', effect size: = ',num2str(ESleft.g1),', CI: ',num2str(ESleft.g1Ci(1)),' to ',num2str(ESleft.g1Ci(2))])
    
    
    dataRight = squeeze(sum(mean(BC_all(fw,rightROI,:),2)));
    [~,Pright,CIright,STATSlright] = ttest(dataRight);
    ESright = mes(dataLeft,0,'g1');
    effectSize = num2str(ESright.g1);
    display(['Odd response, right ROI for ', Cond, ': t = ',num2str(STATSlright.tstat),', p = ',num2str(Pright), ', effect size: = ',num2str(ESright.g1),', CI: ',num2str(ESright.g1Ci(1)),' to ',num2str(ESright.g1Ci(2))])
    
end





