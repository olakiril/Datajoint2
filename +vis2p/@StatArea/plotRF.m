
function plotRF( obj,varargin )

params.trace_opt = 17;

params = getParams(params,varargin);

keys = fetch(obj);
for ikey = 1:length(keys)
    key = keys(ikey)
    key.trace_opt = params.trace_opt;
    cluster =  fetch1(RFMap(key,'rf_opt_num=6 and masknum = 1'),'on_rf');
    [v1o,v2o,v1i,v2i] = fetch1(StatArea(key),'v1o','v2o','v1i','v2i');
    [snr, fitInfo,masknum] = fetchn(RFFit(key,'rf_opt_num=6 and masknum>0'),'snr','gauss_fit','masknum');
    p = fetchn(RFStats(key,'rf_opt_num = 6'),'onpoff_p');
    
    imagesc(zeros(size(cluster)));
    colormap gray
    axis image
    axis off
    hold on
    for icell = 1:length(masknum)
        if sum(masknum(icell)==v1o)
            color= [0 0 1];
            line = '-';
        elseif sum(masknum(icell)==v2o)
            color= [1 0 0];
            line = '-';
        elseif sum(masknum(icell)==v1i)
            color= [0 0.2 0.5];
            line = '-.';
        elseif sum(masknum(icell)==v2i)
            color= [0.5 0.2 0];
            line = '-.';
        end
        if snr(icell)>1.5
            a = fitInfo{icell};
            m=a(1:2); C=diag(a(3:4)); cc=a(5)*sqrt(prod(a(3:4))); C(1,2)=cc; C(2,1)=cc;
            plotGauss(m,C,2,'color',color,'line',line,'linewidth',1);
        end
    end
    title([key.exp_date ' ' num2str(key.scan_idx)])
    if ikey~=length(keys)
        pause
        clf
    end
end