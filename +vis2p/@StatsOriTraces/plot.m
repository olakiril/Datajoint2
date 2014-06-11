function plot(fc,varargin)

% function plot(fc,varargin)
%
% Plots  the tuning curve with the standar errorbars and the raw response
% matrix for a single cell. It shows also the significance of tuning for
% orientation and direction with tuning index and dprime
%
% MF 2011-06-28

params.trace_opt = 6;
params.movie_type = 'natural';

params = getParams(params,varargin);

keys = fetch(fc);

for ikey = 1:length(keys)
    key = keys(ikey);
    key = catstruct(key,params);
    
    [ori Poti Pdoti] = fetch1(StatsOriTraces(key),'ori','Poti','Pdoti');
    AreaMean = mean(ori,1);
    oriNum = length(AreaMean);
    uOri = 0:180/oriNum:180 - 180/oriNum;
      
    errorbar(uOri,AreaMean,std(ori,[],1)/sqrt(size(ori,1)),'.k')
    
    hold on
    m = AreaMean - min(AreaMean);
    orientations = ((1:length(m))-1)*(pi/length(m));
    v = fitVonMisesOne(m,orientations);
    x = 0:.1:pi;
    plot(x/pi*180,vonMises(v,x)+min(AreaMean),'r')
    
    AxisPro = axis;
    Yscale = AxisPro(4)-AxisPro(3);
    text(140,(AxisPro(3)+(Yscale*10)/12),['Pdpr: ' num2str(Pdoti)],'Color',[0.5 0.5 0.5]);
    text(140,(AxisPro(3)+(Yscale*9)/12),['Poti: ' num2str(Poti)],'Color',[0.5 0.5 0.5]);
    title(['Cell:',num2str(key.masknum),', Site:',num2str(key.scan_idx),', Day:',num2str(key.exp_date)])
    
    set(gcf,'Color',[1 1 1])
    set(gca,'Box','Off')
    set(gca,'XLim',[0 180])
    
    if length(keys) > 1
        key %#ok<NOPRT>
        pause
        clf
    end
    
end
function y = vonMises(b,x)
y = b(3) * exp(b(1)*(cos(2*(x-(b(2))))-1)) + b(4)' ;
