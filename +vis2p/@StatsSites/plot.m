function plot( this )
close all;

f1 = figure;
%f2 = figure;
f3 = figure;
Z = zeros(0,2);
Y1 = [];
Y2 = [];
Y3 = [];

for binsize = unique(fetchn(this,'binsize'))'
    for key = fetch( CellTracesGroups*BrainStates.*this )'
        key.binsize = binsize;
        [cellnums,prefOri,oriP] = fetchn( BrainStateCellOriTuning(key),'cellnum','pref_ori','ori_p');
        x1=fetch1( BrainStateCellCovars(setfield(key,'stimulus_condition','spont' )), 'xcov');
        x2=fetch1( BrainStateCellCovars(setfield(key,'stimulus_condition','evoked')), 'xcov');
        x3=fetch1( BrainStateCellCovars(setfield(key,'stimulus_condition','noise_corr')),'xcov');

        % convert cross-covariance to cross-correlation
        c1 = xcov2xcorr(x1);
        c2 = xcov2xcorr(x2);
        c3 = xcov2xcorr(x3);

        L = (size(c1,1)+1)/2;
        N = sqrt(size(c1,2));
        u1 = reshape(c1(L,:),N,N);
        u2 = reshape(c2(L,:),N,N);
        u3 = reshape(c3(L,:),N,N);

        figure( f1 );
        subplot(231);  imagesc( u1 );  axis image;  title('spont');
        subplot(232);  imagesc( u2 );  axis image;  title('evoked');
        subplot(233);  imagesc( u3 );  axis image;  title('noise');
        suptitle( sprintf('mouse %d, %s Scan %d Brainstate "%s", %d-ms bins'...
            ,key.mouse_id, fetch1(Sessions(key),'sess_date'), key.scannum, key.brain_state, key.binsize) );

        [i,j] = ndgrid(1:N,1:N);
        ix = find(i<j);
        Z(end+1,:) = [corr(x1(L,ix)',x2(L,ix)'), corr(x1(L,ix)',x3(L,ix)')];

        d = oriDiff(prefOri(i),prefOri(j));
        i1 = find( i<j & oriP(i)<0.05 & oriP(j)<0.05 & d<10 );
        i2 = find( i<j & oriP(i)<0.05 & oriP(j)<0.05 & d>=10 & d <=30 );
        i3 = find( i<j & oriP(i)<0.05 & oriP(j)<0.05 & d>30 );

        Y1 = [Y1; u1(i1) u2(i1) u3(i1)];  % <10 degrees
        Y2 = [Y2; u1(i2) u2(i2) u3(i2)];  % 10 .. 30 degrees
        Y3 = [Y3; u1(i3) u2(i3) u3(i3)];  % > 30 degrees


        n = sqrt([length(i1) length(i2) length(i3)]);

        subplot(234);
        m1 = [mean(u1(i1)) mean(u1(i2)) mean(u1(i3))];
        bar(m1); hold on;
        errorbar( m1, [std(u1(i1)) std(u1(i2)) std(u1(i3))]./n, 'k.' );

        subplot(235);
        m2 = [mean(u2(i1)) mean(u2(i2)) mean(u2(i3))];
        bar(m2); hold on;
        errorbar( m2, [std(u2(i1)) std(u2(i2)) std(u2(i3))]./n, 'k.' );

        subplot(236);
        m3 = [mean(u3(i1)) mean(u3(i2)) mean(u3(i3))];
        bar(m3); hold on;
        errorbar( m3, [std(u3(i1)) std(u3(i2)) std(u3(i3))]./n, 'k.' );

        mx = max([m1 m2 m3]);
        if isnan(mx)
            continue;
        end
        subplot(234); ylim( [0 mx]*1.2 );
        subplot(235); ylim( [0 mx]*1.2 );
        subplot(236); ylim( [0 mx]*1.2 );

        %         figure( f2 );
        %         plot( CellTracesGroups(key) );
        %
        %         traces = fetchn( CellTraces(key), 'trace' );
        %         traces = [traces{:}];


        figure(f3);
        t = (-(L-1):(L-1))*max(fetchn(BrainStateCellCovars(key),'actual_binsize'));
        subplot(311); plot(t,mean(c1(:,i1),2),'r'); hold on;
        subplot(311); plot(t,mean(c1(:,i2),2),'g');  
        subplot(311); plot(t,mean(c1(:,i3),2));  grid on;  title( 'spont' );
        hold off;
        ylim([-0.2 1.2]*mx);

        subplot(312); plot(t,mean(c2(:,i1),2),'r'); hold on;
        subplot(312); plot(t,mean(c2(:,i2),2),'g');
        subplot(312); plot(t,mean(c2(:,i3),2));  grid on;  title( 'evoked' );
        hold off;
        ylim([-0.2 1.2]*mx);

        subplot(313); plot(t,mean(c3(:,i1),2),'r'); hold on;
        subplot(313); plot(t,mean(c3(:,i2),2),'g');
        subplot(313); plot(t,mean(c3(:,i3),2));  grid on;  title( 'noise' );
        hold off
        ylim([-0.2 1.2]*mx);
        drawnow;
    end
end

figure;
c = 'brg';
for r = 1:3
    m = [mean(Y1(:,r)), mean(Y2(:,r)), mean(Y3(:,r))];
    e = [std(Y1(:,r)), std(Y2(:,r)), std(Y3(:,r))]./sqrt([size(Y1,1),size(Y2,1),size(Y3,1)]);
    errorbar(m,e,c(r));
    hold on;
end
 set(gca,'XTick',1:3,'XTickLabel', {'d<=10', '10<d<=30', 'd>30'});
 legend('spont','evoked','noise');
 title( sprintf('Summary for %s state, binsize=%d', key.brain_state, key.binsize ) );
 ylim([0 0.75]); grid on; box off;
 legend boxoff

disp('done');

function d = oriDiff(a1,a2)
b1 = min(a1,a2);
b2 = max(a1,a2);
d = min(b2-b1,b1+180-b2);


function c=xcov2xcorr(c)
% convert cross-covariance to cross-correlation
N = sqrt(size(c,2));  % number of traces
L = (size(c,1)+1)/2;  % zero-lag element
v = sqrt(c(L,1:N+1:N*N));  % std deviations
k = 1;
for i=1:N
    for j=1:N
        c(:,k)=c(:,k)./(v(i)*v(j));
        k = k + 1;
    end
end