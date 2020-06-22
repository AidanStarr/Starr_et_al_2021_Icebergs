function [Cxy, lag, Q1] = analysis_1(t,x1,x2);
% Cross-correlation 
% This follows the work of Kira Rehfeld (2014). It utilises their NESTool
%toolbox http://tocsy.pik-potsdam.de/nest.php

tx = t(~isnan(x1)); % tx is age with only points sampled for d13C
x = x1(~isnan(x1)); %x is d13C without gaps
ty = t(~isnan(x2)); %ty is age with only points sampled for IRD_MAR
y = x2(~isnan(x2)); %y is IRD_MAR without gaps

h = 0.5; %this sets the kernel width in units of dt (Rehfeld 2011)
lag = [-50 : 1 : 50]; %this sets the lags at which we want to measure correlation 
[Cxy] = nexcf(tx,x,ty,y,lag,h); %the NEXCF function uses gaussian kernel approach

% Monte Carlo routine on surrogate red noise data, this is time consuming
if ~exist('CSur1000.mat') % check to see if there is already a datafile for the surrogate
    NoSur=1000; % number of surrogates for MC analysis. For most reliable results, set to 1000.
    XSur=ar1sur(tx,x,NoSur); % create AR1 Surrogates for X/Y
    YSur=ar1sur(ty,y,NoSur);

    wa=waitbar(0,'Surrogate tests');
    for k=1:NoSur
        waitbar(k/NoSur)
        %gXCF
        [temp1] = nexcf(tx,XSur(:,k),ty,YSur(:,k),lag,h);
        CSur1000(k,:,:)=temp1;
        clear temp1 temp2
    end
    save('CSur1000.mat')
    delete(wa)
else
    load('CSur1000.mat'); %load a pre-run Monte Carlo series
end
% - quantiles for each point from Csur to obtain critical values
Q1=quantile(CSur1000,[.05 .95],1); %find the 5th and 95th percentile 

end