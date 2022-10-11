clear
close all
CarrierFreq = 6;
frequencies =  6./[6,5,4,3,2, 1];
looping = 1;
load Chanlocs
for condition= 1:4
    Cond = ['Cnd_' num2str(condition) '.mat'];
    %% Upload data & pick freq
    ThePath = fullfile('D:\tegolino-derosa\_Experiments\SL30\Clustering and Bayesian analysis\Data\');
    CalcActualF
    load(fullfile(ThePath,Cond))
    
    % Define the electrodes in the ROIs
    Lr = [9,10,11,127,15]; % Left ROI : 9:A9, 10:A10(PO7), 11:A11, 127:D31, 15:A15
    Rr = [38,39,40,43,28]; % Right ROI: 38:B6, 39:B7(PO8), 40:B8, 43:B11,28:A28
    
    for band =1:length(frequencies)
        
        freq = frequencies(band);
        [~,fw] = min(abs(repmat(f,[length(freq),1]) - repmat(freq,[1,length(f)])),[],2);
        %% Left ROI
        %This outputs the measure of interest (SNR or zScore) in the bin of
        %interest for each of the participants for the left ROI
        data(looping).condition = Cond;
        data(looping).frequency = freq;
        data(looping).left =     squeeze(mean(BC_all(fw,Lr,:),2));
        data(looping).right = squeeze(mean(BC_all(fw,Rr,:),2));
        looping = looping+1;
    end
end


clc
for har = 1:length(data)
tempdata = data(har).left;
[h4,P,CI,STATS] = ttest(tempdata,0,'Tail', 'right');
Effect = mes(tempdata,0,'g1');
effectSize = num2str(Effect.g1);
disp(['Mean response in the left ROI for ', data(har).condition, ', frequency:' num2str(data(har).frequency) ': t(28) = ',num2str(STATS.tstat), ', p = ',num2str(P), ', effect size: = ',num2str(Effect.g1),', [',num2str(Effect.g1Ci(1)) ' ' num2str(Effect.g1Ci(2)) ']']);
disp(h4)
tempdata = data(har).right;
[h4,P,CI,STATS] = ttest(tempdata,0,'Tail', 'right');
Effect = mes(tempdata,0,'g1');
effectSize = num2str(Effect.g1);
disp(['Mean response in the right ROI for ', data(har).condition, ', frequency:' num2str(data(har).frequency) ': t(28) = ',num2str(STATS.tstat), ', p = ',num2str(P), ', effect size: = ',num2str(Effect.g1),', [',num2str(Effect.g1Ci(1)) ' ' num2str(Effect.g1Ci(2)) ']']);
disp(h4)

end

left = nan(24,29);


% l =1;
% for i=1:2:48
% all(i,1:29) = data(l).left;
% all(i+1,1:29) = data(l).right;
% l = l+1;
% end


for i=1:24
left(i,1:29) = data(i).left;
end

for i=1:24
right(i,1:29) = data(i).right;
end
Cond='C3';cond=12;
h4 = notBoxPlot(left(1+cond:5+cond,:)','jitter',.5,'interval','tInterval');
%h4 = notBoxPlot(right(1+cond:5+cond,:)','jitter',.5,'interval','tInterval');
%xticklabels({'1 Hz', '1.2 Hz', '1.5 Hz', '2 Hz', '3 Hz', '6Hz'})
xticklabels({'1 Hz', '1.2 Hz', '1.5 Hz', '2 Hz', '3 Hz'})

ylim([-0.1,0.25]);
yline(0,'k:','LineWidth',2)
hold on
d=[h4.data];
set(d,'markerfacecolor',[.8,.8,.8],'color',[0,0,0])
set(d,'MarkerSize', 4)

set([h4.mu],'Color','k'); set([h4.semPtch], 'FaceColor',[1,1,1]*0.75,'EdgeColor','none'); set([h4.sdPtch],'FaceAlpha',0);set([h4.sdPtch],'LineStyle','none')

if (Cond == "C1")
    hue = [0.4660 0.6740 0.1880];
elseif (Cond == "C2")
    hue = [0, 0.4470, 0.7410];
elseif  (Cond == "C3")
    hue = [0.9290 0.6940 0.1250];
elseif (Cond == "C4")
    hue = [0.8500 0.3250 0.0980];
end
%set(d(2),'marker', 'o') %,
set([h4(2).semPtch],'FaceColor', hue,'EdgeColor','none')
set([h4(2).semPtch],'FaceColor', hue,'EdgeColor','none')
set(d(2),'markerfacecolor',hue)%'color', 'markerfacecolor',[0,0,0])

saveas(gcf,[Cond 'left.png'])
%saveas(gcf,[Cond 'right.png'])









cond = 0;
for condition= 1:4
    Cond = ['C' num2str(condition)];    
        if (Cond == "C1")
        hue = [0.4660 0.6740 0.1880];
    elseif (Cond == "C2")
        hue = [0, 0.4470, 0.7410];
    elseif  (Cond == "C3")
        hue = [0.9290 0.6940 0.1250];
    elseif (Cond == "C4")
        hue = [0.8500 0.3250 0.0980];
    end

    h4 = notBoxPlot(left(1+cond:5+cond,:)','jitter',.5,'interval','tInterval');
    xticklabels({'1 Hz', '1.2 Hz', '1.5 Hz', '2 Hz', '3 Hz'})   
    ylim([-0.1,0.25]);
    hold on;    d=[h4.data]; set(d,'markerfacecolor',[.8,.8,.8],'color',[0,0,0]); set(d,'MarkerSize', 4);    
    set([h4.mu],'Color','k'); set([h4.semPtch], 'FaceColor',[1,1,1]*0.75,'EdgeColor','none'); set([h4.sdPtch],'FaceAlpha',0);set([h4.sdPtch],'LineStyle','none')
    set([h4(2).semPtch],'FaceColor', hue,'EdgeColor','none'); set([h4(2).semPtch],'FaceColor', hue,'EdgeColor','none'); set(d(2),'markerfacecolor',hue);
    saveas(gcf,[Cond 'left_nodots.png']); close;
    
    h4 = notBoxPlot(right(1+cond:5+cond,:)','jitter',.5,'interval','tInterval');
    xticklabels({'1 Hz', '1.2 Hz', '1.5 Hz', '2 Hz', '3 Hz'})   
    ylim([-0.1,0.25]);
    hold on;    d=[h4.data]; set(d,'markerfacecolor',[.8,.8,.8],'color',[0,0,0]); set(d,'MarkerSize', 4);    
    set([h4.mu],'Color','k'); set([h4.semPtch], 'FaceColor',[1,1,1]*0.75,'EdgeColor','none'); set([h4.sdPtch],'FaceAlpha',0);set([h4.sdPtch],'LineStyle','none')
    set([h4(2).semPtch],'FaceColor', hue,'EdgeColor','none'); set([h4(2).semPtch],'FaceColor', hue,'EdgeColor','none'); set(d(2),'markerfacecolor',hue);
    saveas(gcf,[Cond 'right_nodots.png']); close;
    cond = cond +6;
end