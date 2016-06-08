lamdas = -(0.2:0.025:0.55);
sigmas = 1:0.1:2;
pcas = [1 30 60 90];
highs = 0.01:0.02:0.2;
gains = 0: 0.25: 2;

key = []; tpO = struct;
key.exp_date = '2012-07-24';
key.scan_idx = 13;
Pmask =fetch1(EphysTraces(key),'masknum');
traces = fetchn(MaskTraces(key,'masknum>0'),'calcium_trace');
fps = fetch1(Movies(key),'fps');
traces = double([traces{:}]);  % traces are in columns
gspikes = fetch1(Traces('trace_opt = 17 and masknum = 0',key),'trace');
[c1r11, c1r12, c1r21, c1r22 ,q1r1 ,q1r2] = initialize('cell',length(sigmas),1);
parfor isig = 1:length(sigmas);tp = tpO;   tp.sigma = sigmas(isig);
    for ilamda = 1:length(lamdas);         tp.lamda = lamdas(ilamda);
        for igain = 1:length(gains);       tp.gain = gains(igain);
            for ipca = 1:length(pcas);     tp.pca = pcas(ipca);
                for ihp= 1:length(highs);  tp.highPass = highs(ihp);
                    [tracesP2, q1r2{isig}(ilamda,ihp,igain,ipca)] = oopsi(traces(:,Pmask),fps,tp);
                    [tracesP,out,tp] = pcMinusHighPass(traces,fps,tp);
                    [tracesP, q1r1{isig}(ilamda,ihp,igain,ipca)] = oopsi(tracesP(:,Pmask),fps,tp);
                    c1r11{isig}(ilamda,ihp,igain,ipca) = corr(tracesP(2:end),gspikes(1:end-1));
                    c1r12{isig}(ilamda,ihp,igain,ipca) = corr(tracesP,gspikes);
                    c1r21{isig}(ilamda,ihp,igain,ipca) = corr(tracesP2(2:end),gspikes(1:end-1));
                    c1r22{isig}(ilamda,ihp,igain,ipca) = corr(tracesP2,gspikes);
                end
            end
        end
    end
end
save('data1','c1r11', 'c1r12', 'c1r21', 'c1r22' ,'q1r1' ,'q1r2')

key = [];
key.exp_date = '2012-06-29';
key.scan_idx = 26;
Pmask =fetch1(EphysTraces(key),'masknum');
traces = fetchn(MaskTraces(key,'masknum>0'),'calcium_trace');
fps = fetch1(Movies(key),'fps');
traces = double([traces{:}]);  % traces are in columns
gspikes = fetch1(Traces('trace_opt = 17 and masknum = 0',key),'trace');
[c2r11, c2r12, c2r21, c2r22 ,q2r1 ,q2r2] = initialize('cell',length(sigmas),1);
parfor isig = 1:length(sigmas);tp = tpO;   tp.sigma = sigmas(isig);
    for ilamda = 1:length(lamdas);         tp.lamda = lamdas(ilamda);
        for igain = 1:length(gains);       tp.gain = gains(igain);
            for ipca = 1:length(pcas);     tp.pca = pcas(ipca);
                for ihp= 1:length(highs);  tp.highPass = highs(ihp);
                    [tracesP2, q2r2{isig}(ilamda,ihp,igain,ipca)] = oopsi(traces(:,Pmask),fps,tp);
                    [tracesP,out,tp] = pcMinusHighPass(traces,fps,tp);
                    [tracesP, q2r1{isig}(ilamda,ihp,igain,ipca)] = oopsi(tracesP(:,Pmask),fps,tp);
                     c2r11{isig}(ilamda,ihp,igain,ipca) = corr(tracesP(2:end),gspikes(1:end-1));
                    c2r12{isig}(ilamda,ihp,igain,ipca) = corr(tracesP,gspikes);
                    c2r21{isig}(ilamda,ihp,igain,ipca) = corr(tracesP2(2:end),gspikes(1:end-1));
                    c2r22{isig}(ilamda,ihp,igain,ipca) = corr(tracesP2,gspikes);
                end
            end
        end
    end
end
save('data2','c2r11', 'c2r12', 'c2r21', 'c2r22' ,'q2r1' ,'q2r2')

key = []; tpO = struct;
key.exp_date = '2012-07-28';
key.scan_idx = 1;
Pmask =fetch1(EphysTraces(key),'masknum');
traces = fetchn(MaskTraces(key,'masknum>0'),'calcium_trace');
fps = fetch1(Movies(key),'fps');
traces = double([traces{:}]);  % traces are in columns
gspikes = fetch1(Traces('trace_opt = 17 and masknum = 0',key),'trace');
[c3r11, c3r12, c3r21, c3r22 ,q3r1 ,q3r2] = initialize('cell',length(sigmas),1);
parfor isig = 1:length(sigmas);tp = tpO;   tp.sigma = sigmas(isig);
    for ilamda = 1:length(lamdas);         tp.lamda = lamdas(ilamda);
        for igain = 1:length(gains);       tp.gain = gains(igain);
            for ipca = 1:length(pcas);     tp.pca = pcas(ipca);
                for ihp= 1:length(highs);  tp.highPass = highs(ihp);
                    [tracesP2, q3r2{isig}(ilamda,ihp,igain,ipca)] = oopsi(traces(:,Pmask),fps,tp);
                    [tracesP,out,tp] = pcMinusHighPass(traces,fps,tp);
                    [tracesP, q3r1{isig}(ilamda,ihp,igain,ipca)] = oopsi(tracesP(:,Pmask),fps,tp);
                    c3r11{isig}(ilamda,ihp,igain,ipca) = corr(tracesP(2:end),gspikes(1:end-1));
                    c3r12{isig}(ilamda,ihp,igain,ipca) = corr(tracesP,gspikes);
                    c3r21{isig}(ilamda,ihp,igain,ipca) = corr(tracesP2(2:end),gspikes(1:end-1));
                    c3r22{isig}(ilamda,ihp,igain,ipca) = corr(tracesP2,gspikes);
                end
            end
        end
    end
end
save('data3','c3r11', 'c3r12', 'c3r21', 'c3r22' ,'q3r1' ,'q3r2')


%%
figure
c = reshape(cell2mat(cr11),length(lamdas),length(sigmas),length(highs));
cc = C26{1,1};
for i = 1:size(c,1)
    subplot(3,3,i)
    imagesc(squeeze(c(i,:,:)).*squeeze(cc(i,:,:)))
    colorbar
end
%%
c = squeeze(cc(:,:,:,7,5));
for i = 1:size(cc,1)
    subplot(3,3,i)
    imagesc(squeeze(c(i,:,:)))
    colorbar
end
%%
c1 = reshape(cell2mat(c1r11),length(lamdas),length(sigmas),length(highs),length(gains),length(pcas));
c2 = reshape(cell2mat(c2r11),length(lamdas),length(sigmas),length(highs),length(gains),length(pcas));
c3 = reshape(cell2mat(c3r11),length(lamdas),length(sigmas),length(highs),length(gains),length(pcas));
cc = c1.*c2.*c3;
[~,ic] = max(cc(:));
[x, y, z, p, o] = ind2sub(size(cc),ic);
[lamdas(x) sigmas(y) highs(z) gains(p) pcas(o)]