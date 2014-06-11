function plot(obj)

figure
    plot(flipTimes,ones(size(flipTimes)),'.');
    hold on;
    colors = hsv(length(stimfiles));
    for iStim = 1:length(stimfiles)
        correctedSwaps =  swapsTest{iStim};
        plot(correctedSwaps,(iStim+1)*ones(size(correctedSwaps)),'.','Color',colors(iStim,:))
    end
    set(gca,'YLim',[0 100])
    title(tprname);
% end