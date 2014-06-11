function testplot(fc,varargin)

% function VMPlot(fc,varargin)
%
% Plots  the tuning curve with the standar errorbars and the raw response
% matrix for a single cell. It shows also the significance of tuning for
% orientation and direction with tuning index and dprime
%
% MF 2011-06-28

params.text = 1;

params = getParams(params,varargin);
keys = fetch(fc);

for ikey = 1:length(keys)
    clf
    key = keys(ikey);
    cases{1} = key;
    cases{1}.movie_type = 'natural';
    cases{1}.trace_opt = 3;
    cases{2} = key;
    cases{2}.movie_type = 'phase';
    cases{2}.trace_opt = 3;
    cases{3} = key;
    cases{3}.movie_type = 'natural';
    cases{3}.trace_opt = 6;
    cases{4} = key;
    cases{4}.movie_type = 'phase';
    cases{4}.trace_opt = 6;
    for icase = 1:4
        
        [ori Poti Pdoti] = fetch1(StatsOriTraces(cases{icase}),'ori','Poti','Pdoti');
        AreaMean = mean(ori,1);
        oriNum = length(AreaMean);
        uOri = 0:180/oriNum:180 - 180/oriNum;
        
        subplot(2,2,icase)
        errorbar(uOri,AreaMean,std(ori,[],1)/sqrt(size(ori,1)),'.k')
        
        hold on
        m = AreaMean - min(AreaMean);
        orientations = ((1:length(m))-1)*(pi/length(m));
        v = fitVonMisesOne(m,orientations);
        x = 0:.1:pi;
        plot(x/pi*180,vonMises(v,x)+min(AreaMean),'r')
        
        if params.text
            AxisPro = axis;
            Yscale = AxisPro(4)-AxisPro(3);
            text(140,(AxisPro(3)+(Yscale*10)/12),['Pdpr: ' num2str(Pdoti)]);
            text(140,(AxisPro(3)+(Yscale*9)/12),['Poti: ' num2str(Poti)]);
            title(['Cell:',num2str(key.masknum),', Site:',num2str(key.scan_idx),', Day:',num2str(key.exp_date)])
        end
        set(gcf,'Color',[1 1 1])
        set(gca,'Box','Off')
    end
    key %#ok<NOPRT>
    pause
end
function y = vonMises(b,x)
y = b(3) * exp(b(1)*(cos(2*(x-(b(2))))-1)) + b(4)' ;
