function plot(fc,varargin)

% function VMPlot(fc,varargin)
%
% Plots  the tuning curve with the standar errorbars and the raw response
% matrix for a single cell. It shows also the significance of tuning for
% orientation and direction with tuning index and dprime
%
% MF 2011-06-28

params.trials = 'ste';

params = getParams(params,varargin);

keys = fetch(fc);

for ikey = 1:length(keys)
    key = keys(ikey);
    
    [ori Poti Pdoti] = fetch1(OriTraces(key),'ori','Poti','Pdoti');
    stim = fetch1(VisStims(key),'stim_file');
    AreaMean = mean(ori,1);
    oriNum = length(AreaMean);
    uOri = 0:360/oriNum:360 - 360/oriNum;
    
    if strcmp(params.trials,'ste')
        errorbar(uOri,AreaMean,std(ori,[],1)/sqrt(size(ori,1)),'.k')
    elseif  strcmp(params.trials,'std')
        errorbar(uOri,AreaMean,std(ori,[],1),'.k')
    else
        plot(repmat(uOri,size(ori,1),1),ori,'.k')
    end
    
    hold on
    m = double(AreaMean - min(AreaMean));
    orientations = stim.params.constants.orientation/360 * 2 * pi;
    v = fitVonMises(m,orientations);
    x = 0:.2:2*pi;
    plot(x/(2*pi)*360,vonMises(v,x)+min(AreaMean),'r')
    
    AxisPro = axis;
    Yscale = AxisPro(4)-AxisPro(3);
    text(300,(AxisPro(3)+(Yscale*10)/12),['Pdpr: ' num2str(Pdoti)],'Color',[0.5 0.5 0.5]);
    text(300,(AxisPro(3)+(Yscale*9)/12),['Poti: ' num2str(Poti)],'Color',[0.5 0.5 0.5]);
    title(['Cell:',num2str(key.masknum),', Site:',num2str(key.scan_idx),', Day:',num2str(key.exp_date)])
    
    set(gcf,'Color',[1 1 1])
    set(gca,'Box','Off')
    set(gca,'XLim',[0 360])
    
    if length(keys) > 1
        key %#ok<NOPRT>
        pause
        clf
    end
    
end
function y = vonMises(b,x)
y = b(3) * exp(b(1)*(cos((x-(b(2))))-1)) + b(4)' + b(5) * exp(b(1)*(cos((pi+x-(b(2))))-1));
