function [x2lags, Q, q] = analysis_3(t,dt,age,x1,x2,filtype,filtdeg,dif,max_dist);

%---- preamble ----%
max_distance = max_dist; % maximum distance to search between x1 and x2 peaks
pdist1 = 15; pdist2 = 15;% set minimum peak distance for x1 (1) and x2 (2) 
pheight1 = 0.075; pheight2 = 0.01; % set min peak height for x1 and x2: reccomended 0.075 for irdmar and 0.01 for bd13c

%---- Step 1: linear interpolation ----%
tl = fliplr(t.*-1)';
age1 = age(~isnan(x1));
x1 = x1(~isnan(x1));
xl1 = interp1(-age1,x1,tl,'linear','extrap');
age2 = age(~isnan(x2));
x2 = x2(~isnan(x2));
xl2 = interp1(-age2,x2,tl,'linear','extrap');

xl1 = normalize(xl1); % zscore normalization
xl2 = normalize(xl2); % zscore normalization

%---- Step 2: Filter ----%
if filtype == 1
    SM = filtdeg;
    f1s = filtfilt (ones(1,round(SM/dt))/round(SM/dt),1,xl1); % Moving-average
    f2s = filtfilt (ones(1,round(SM/dt))/round(SM/dt),1,xl2); % Moving average
end
if filtype == 2
    sgolay_degree = filtdeg;
    sgolay_len = 31; % Savitsky-Golay Filter Length
    f1s = sgolayfilt(xl1,sgolay_degree,sgolay_len); % Savitsky-Golay method
	f2s = sgolayfilt(xl2,sgolay_degree,sgolay_len); % Savitsky-Golay method
end
%---- Step 3: First Difference ----%
if dif == 1
    tDiff = tl(1:length(tl)-1,:) + diff(tl)/2; % generic timescale for 1st differential
    for i = 1 : length(tl)-1
        x1diff(i) = (f1s(i+1) - f1s(i))./(tl(i+1)-tl(i));
        x2diff(i) = (f2s(i+1) - f2s(i))./(tl(i+1)-tl(i));
    end
end

%---- Step 4: findpeaks ----%
[pks1,lcs1]=findpeaks(x1diff,tDiff,'minpeakdistance',pdist1,'minpeakheight',pheight1); % find x1 peaks
[pks2,lcs2]=findpeaks(x2diff,tDiff,'minpeakdistance',pdist2,'minpeakheight',pheight2); % find x2 peaks

for i = 1:length(pks2)
    l = lcs2(i);
    l2 = find(lcs1 >= l-max_distance & lcs1 < l+max_distance,1,'first');
    if isempty(l2)
        bin1(i) = NaN; % this means there is not ird peak within the max distance of the 13c peak
    else
        bin1(i)=lcs2(i)-lcs1(l2);
    end
end

%---- Step 5: checl with event synchronization (Rehfeld 2013) ----%
for i = 0:5
[Q(i+1),q(i+1)] = eventsynchro(tDiff,x1diff,tDiff,x2diff,i,0.9);
end

display(['Events Found: ',num2str(length(bin1))])
display(['Average Lag: ',num2str(nanmean(-bin1))])

x2lags = -bin1;
end

    

