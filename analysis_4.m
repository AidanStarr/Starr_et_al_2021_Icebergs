function [d18o_cum, ird_cum, time_cum, time_cum_norm, f3s, tl] = analysis_4(t,dt,age,d18o,x1,plt);

%---- Step 1: linear interpolation of bd18O ----%
tl = fliplr(t.*-1)';
age1 = age(~isnan(d18o));
x3 = d18o(~isnan(d18o));
x3l = interp1(-age1,x3,tl,'linear','extrap');
age2 = age(~isnan(x1));
x1 = x1(~isnan(x1));
x1l = interp1(-age2,x1,tl,'linear','extrap');

%---- Step 2: Filter bd18O ----%
SM = 10; % moving average filter degree
f3s = filtfilt (ones(1,round(SM/dt))/round(SM/dt),1,x3l); % Moving average


%---- Step 3: Divide d18O into Glacials ----%
lcs = load('glac_load.mat');
pk = lcs.glac_load;

%---- Step 4: Integrate under d18O curve for each glacial ----%
for i = 1 : length(pk)
    n1 = find(tl>=pk(i),1,'first');
    if i == length(pk)
        n2 = find((tl<=0) & (f3s>=(3.75-0.43)),1,'last'); % run between peak interglacial and the next time 18O crosses a threshold
    else
    n2 = find((tl<=pk(i+1)) & (f3s>=(3.75-0.43)),1,'last'); % run between peak interglacial and the next time 18O crosses a threshold
    end
    if isempty(n1)
        cp{i} = [];
        ircp{i} = [];
        d18c{i} = [];
    end
    if isempty(n2)
        cp{i} = [];
        ircp{i} = [];
        d18c{i} = [];
    end
    if (~isempty(n1)) && (~isempty(n2))
        xx = normalize(x3l(n1:n2),'range',[0 1]);
        cp{i} = cumtrapz(tl(n1:n2),xx); %cumtrapz integrates and records the results in cp
        ircp{i} = cumtrapz(tl(n1:n2),x1l(n1:n2));
        d18c{i} = x3l(n1:n2);
    end
    tp{i} = tl(n1:n2);
end

%---- Step 5: 'normalise' the results to percentages of total accumulation ----%
for i = 1:length(cp)
    cpp{i} = (cp{i}./(max(cp{i}))).*100;
    ircpp{i} = (ircp{i}./(max(ircp{i}))).*100;
end
% do the same for time
for i = 1:length(tp)
    tpp{i} = (tp{i}-min(tp{i}));
    tpp{i} = (tpp{i}./(max(tpp{i}))).*100;
end

%---- Step 6: Prepare outputs ----%
d18o_cum = cpp;
ird_cum = ircpp;
time_cum = tp;
time_cum_norm = tpp;

end