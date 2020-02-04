%--- Master script for Starr et al. in review, 'Southern Ocean Lead' ---%
% Jan 2020
%- Citation: 
%              Aidan Starr*, Ian R. Hall, Stephen Barker, Sidney R.
%              Hemming, Jeroen van der Lubbe, Melissa A. Berke, Grant R.
%              Bigg, Alejandra Cartagena, Francisco J. Jiménez-Espejo, Jens
%              Gruetzner, Gregor Knorr, Nambiyathodi Lathika, Leah J.
%              LeVay, Rebecca Robinson, Martin Ziegler, Xu Zhang, Exp.
%              361 Science Party. 
%              "A Southern Ocean Lead Over Deep Water-Mass Reorganizations 
%              during Pleistocene Glacials". 
%              In Review, Nature
%
%- Corresponding Author*: Aidan Starr, StarrA1@Cardiff.ac.uk  
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
% bd13c = benthic d13C
% bd18o = benthic d18O
% irdmar = IRD apparent mass accumulation (grains/cm2/kyr)
% age = Age (kyrs)
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
filtdeg = 7; % if 1, this is moving average degree in dt
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
filtdeg = 7; % if 1, this is moving average degree in dt
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
