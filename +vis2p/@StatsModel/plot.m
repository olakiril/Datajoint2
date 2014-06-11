function plot(obj,key,varargin)

params.type = 'performance';
params.linestyle = '-';
params.colors = [];
params.chance = 1;

params = getParams(params,varargin);

[popSprsN date scan] = fetchn(StatsModel(key),params.type,'exp_date','scan_idx');

if isempty(params.colors)
    colors = hsv(length(popSprsN));
else
    colors = params.colors;
end
for i = 1:length(popSprsN)
    hold on;
    plot(popSprsN{i},'linestyle',params.linestyle,'LineWidth',3,'color',colors(i,:));
    ylabel('Classification Performance','FontSize',9)
    xlabel('# of neurons','FontSize',9)
    set(gcf,'Color',[1 1 1])
    set(gca,'Box','Off')
    set(gca,'FontSize',9)
end
title([date{1} '  scan:  ' num2str(scan(1))],'fontsize',11,'fontweight','bold')
set(gcf,'name',['Classification performance ' date{1} ' ' num2str(scan(1))])

key.class_opt = 9;
if ~isempty(StatsModel(key)) && params.chance
    params.linestyle = '-.';
    [popSprsN date scan] = fetchn(StatsModel(key),params.type,'exp_date','scan_idx');
    for i = 1:length(popSprsN)
        hold on;
        plot(popSprsN{i},'linestyle',params.linestyle,'LineWidth',3,'color',colors(i,:));
        ylabel('Classification performance','FontSize',9)
    end
end