%--- Master script for Starr et al. in review, 'Shifting Antarctic Iceberg Melt' ---%
% June 2020
%- Citation: 
%              Aidan Starr*, Ian R. Hall, Stephen Barker, Sidney R.
%              Hemming, Jeroen van der Lubbe, Melissa A. Berke, Grant R.
%              Bigg, Alejandra Cartagena, Francisco J. Jiménez-Espejo, Jens
%              Gruetzner, Gregor Knorr, Nambiyathodi Lathika, Leah J.
%              LeVay, Rebecca Robinson, Martin Ziegler, Xu Zhang, Exp.
%              361 Science Party. 
%              "Shifting Antarctic Iceberg Melt Leads Glacial Deep Water-Mass Reorganizations". 
%              In Review, Nature
%
%- Corresponding Author*: Aidan Starr, StarrA1@Cardiff.ac.uk / aidan.m.starr@gmail.com  
%
%
%- Dependencies: nexcf.m; ar1sur.m;
%- Calls: analysis_1.m; analysis_2.m; analysis_3.m; analysis_4.m;
%        load_data.m;  create_analyseries_files.m
%
% Note: each part of the analysis is designed as a function. Caution is
% advised when running some sections as Monte Carlo simulations may be CPU
% intensive, to save time for users, a 'here's one I made earlier'
% simulation .mat file is available on my github
%
%--- Master script for Starr et al. 'Southern Ocean Lead' ---%
% Citation: 
%
%- Dependencies: nexcf.m; ar1sur.m; ETP.mat; glac_load.mat;
%- Calls: analysis_1.m; analysis_2.m; analysis_3.m; analysis_4.m;
%        load_data.m;  create_analyseries_files.m
%
% Note: each part of the analysis is designed as a function, figures are then
% generated seperately at the end of the script. Caution is advised when
% running some section as Monte Carlo simulations may be CPU intensive, to
% save some power, pre-run simulation .mat files are available on my github

%%%%--------- Preamble ---------%%%%%
% load the data available from pangaea.de (DOI xxx)
cd('/home/pablo/Insync/StarrA1@cardiff.ac.uk/OneDrive Biz/Paper 1/repo/Southern_Ocean_Lead');
addpath(genpath(pwd))
%path = 'https://github.com/AidanStarr/Southern_Ocean_Lead/raw';
data = readtable('raw/starr_submitted.csv');
bd13c = data.C_WuellerstorfiD13C__Vs_PDB_;
bd18o = data.C_WuellerstorfiD18O__Vs_PDB_;
irdmar = data.IRDApparentMAR_grains_cm2_kyr_;
age = data.Age_kyr_;

dt = 1.5; % set dt for linearly rescaled data (used later)
lt = [0 : dt : 1645]; % create linear time scale in kyrs

%%%%--------- Analysis 1: Cross-Correlation ---------%%%%%
[Cxy1, lag, Q1] = analysis_1(age,bd13c,log10(irdmar+1)); % --> analysis_1(time vector, variable 1, variable 2)
[Cxy2, ~, ~] = analysis_1(age,bd18o,log10(irdmar+1)); % --> analysis_1(time vector, variable 1, variable 2)

%%%%--------- Analysis 2: Frequency Analysis ---------%%%%%
% Get lags from Blackman-Tukey calculated in Analyseries
[ob_lag, ec_lag, pr_lag] = analysis_2('irdmar','bd13c'); % --> analysis_2('variable 1','variable 2');

%%%%--------- Analysis 3: Peak-Lag Algorithm ---------%%%%%
% USAGE: analysis_3(linear timescale, spacing, original time scale, variable 1, variable 2)
filtype = 1; % 1 means use Moving Average filter, 2 means use Savitsky-Golay Filter
filtdeg = 6; % if 1, this is moving average degree in dt
dif = 1; % use first difference?
max_dist = 10;
[bd13c_lags, Q, q] = analysis_3(lt,dt, age,log10(irdmar+1),-bd13c,filtype,filtdeg,dif,max_dist);


% Test sensitivity to filter degree
filtype = 1; % 1 means use Moving Average filter, 2 means use Savitsky-Golay Filter
for i = 1 : 11
    filtdeg = i; % if 1, this is moving average in kyrs, if 2, this is SGolay Degrees
    dif = 1; % Use First Difference of x1 and x2? 1 for yes, 0 for no
    temp = analysis_3(lt,dt, age,log10(irdmar+1),-bd13c,filtype,filtdeg,dif,max_dist);
    b_lag(i)=nanmean(temp);
end


%%%%--------- Analysis 3b: Red-Noise Peak-Lag Algorithm ---------%%%%%
% Create red-noise surrogates and run them through the peak-lag algorithm
 NoSur=1000; % number of surrogates for MC analysis. For most reliable results, set to 1000.
 XSur=ar1sur(age(~isnan(bd13c)),bd13c(~isnan(bd13c)),NoSur); % create AR1 Surrogates for bd13c
 filtype = 1; % 1 means use Moving Average filter, 2 means use Savitsky-Golay Filter
filtdeg = 6; % if 1, this is moving average degree in dt
dif = 1; % use first difference?
    for k=1:NoSur-1
        temp = analysis_3(lt,dt,age(~isnan(bd13c)),XSur(:,1),XSur(:,k+1),filtype,filtdeg,dif,max_dist);
        red_lags{k}=temp;
        clear temp
    end
for i = 1 : size(red_lags,2)
    red_lag(i) = nanmean(red_lags{i});
end
% use a 2-sided t test to determine whether the ird vs 13C data are
% significantly different to the 'random' red noise surrogates
[h,p,ci] = ttest2(bd13c_lags,red_lag);

%%%%--------- Analysis 4: Glacial Accumulation ---------%%%%%
plt = 0; % plot the process? 1 for yes, 0 for no
[d18o_cum, ird_cum, time_cum, time_cum_norm, f3s, tl] = analysis_4(lt,dt,age,bd18o,irdmar,plt); % Divide d18O into glacials  and intergrate under curve
% split into intervals
enwx = 1:1:100';
for i = 1:23
    newird_cum(:,i) = interp1(time_cum_norm{i},ird_cum{i},enwx,'linear');
end

for i = 1:23
    newtime_cum(:,i) = interp1(time_cum_norm{i},time_cum_norm{i},enwx,'linear');
end

for i = 1:23
    newd18o_cum(:,i) = interp1(time_cum_norm{i},d18o_cum{i},enwx,'linear');
end
% mean accumulation for IRD 
A = newird_cum(:,1:12);
B = newird_cum(:,13:18);
C = newird_cum(:,19:23);
pre_mpt = mean(A,2);
mpt = mean(B,2);
post_mpt = mean(C,2);
% calculate confidence intervals for IRD 
SEM_A = std(A,1,2)/sqrt(40); % get standard error of the mean
SEM_B = std(B,1,2)/sqrt(40);
SEM_C = std(C,1,2)/sqrt(40);
CI95_A = bsxfun(@plus, mean(A,2), bsxfun(@times, [-1  1]*1.96, SEM_A));   % 95% Bootstrap Confidence Intervals
CI95_B = bsxfun(@plus, mean(B,2), bsxfun(@times, [-1  1]*1.96, SEM_B));   
CI95_C = bsxfun(@plus, mean(C,2), bsxfun(@times, [-1  1]*1.96, SEM_C));   % 95% Confidence Intervals
% mean accumulation for d18O
A2 = newd18o_cum(:,1:11);
B2 = newd18o_cum(:,2:18);
C2 = newd18o_cum(:,19:23);
pre_mptd = mean(A2,2);
mptd = mean(B2,2);
post_mptd = mean(C2,2);

% calculate confidence intervals for d18O 
SEM_A2 = std(newd18o_cum(:,1:11),1,2)/sqrt(40);
SEM_B2 = std(newd18o_cum(:,12:18),1,2)/sqrt(40);
SEM_C2 = std(newd18o_cum(:,19:23),1,2)/sqrt(40);
CI95_A2 = bsxfun(@plus, mean(A2,2), bsxfun(@times, [-1  1]*1.96, SEM_A2));   % 95% Confidence Intervals
CI95_B2 = bsxfun(@plus, mean(B2,2), bsxfun(@times, [-1  1]*1.96, SEM_B2));   % 95% Confidence Intervals
CI95_C2 = bsxfun(@plus, mean(C2,2), bsxfun(@times, [-1  1]*1.96, SEM_C2));   % 95% Confidence Intervals

% analysis_4 OUTPUT
output_cumtrap(:,1) = pre_mpt;
output_cumtrap(:,2) = mpt;
output_cumtrap(:,3) = post_mpt;

output_cumtrap(:,4) = CI95_A(:,1);
output_cumtrap(:,5) =  CI95_A(:,2);
output_cumtrap(:,6) =  CI95_B(:,1);
output_cumtrap(:,7) =  CI95_B(:,2);
output_cumtrap(:,8) =  CI95_C(:,1);
output_cumtrap(:,9) =  CI95_C(:,2);

output_cumtrap(:,10) = pre_mptd;
output_cumtrap(:,11) = mptd;
output_cumtrap(:,12) = post_mptd;

output_cumtrap(:,13) = CI95_A2(:,1);
output_cumtrap(:,14) =  CI95_A2(:,2);
output_cumtrap(:,15) =  CI95_B2(:,1);
output_cumtrap(:,16) =  CI95_B2(:,2);
output_cumtrap(:,17) =  CI95_C2(:,1);
output_cumtrap(:,18) =  CI95_C2(:,2);

output_cumtrap(:,19) = enwx';
output_4 = array2table(output_cumtrap,'variablenames',{'IRD_mean_pre','IRD_mean_mpt','IRD_mean_post','IRDu_pre','IRDl_pre','IRDu_mpt','IRDl_mpt','IRDu_post','IRDl_post','d18o_mean_pre','d18o_mean_mpt','d18o_mean_post','d18ou_pre','d18ol_pre','d18ou_mpt','d18ol_mpt','d18ou_post','d18ol_post','time'});
writetable(output_4,'outputs/cumulative_glacials.csv')

%%%%--------- Output Table ---------%%%%
% column 1 is the peak-lag algorithm result
output{1,1} = nanmean(bd13c_lags); % mean
output{2,1} = ci(1); % upper confidence
output{3,1} = ci(2); % lower confidence
% columns 2,3,4 are the blackman-tukey phase estimates
output{1,2} = -1*pr_lag(1);
output{2,2} = -1*pr_lag(3);
output{3,2} = -1*pr_lag(2);
output{1,3} = -1*ob_lag(1);
output{2,3} = -1*ob_lag(3);
output{3,3} = -1*ob_lag(2);
output{1,4} = -1*ec_lag(1);
output{2,4} = -1*ec_lag(3);
output{3,4} = -1*ec_lag(2);
%column 5 is the cross-correlaion estimate
output{1,5} = lag(Cxy1==min(Cxy1));
output{2,5} = lag(Cxy1==min(Cxy1))-mean(diff(lag));
output{3,5} = lag(Cxy1==min(Cxy1))+mean(diff(lag));

output_table = cell2table(output,'variablenames',{'peak_lag','precession_BT','obliquity_BT','eccentricity_BT','cross_correlation'});
writetable(output_table,'outputs/lag_results.csv')