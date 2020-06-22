function [ob, ec, pr] = analysis_2(x1,x2)

% Here I import phase and coherence estimates computed using the
% Blackman-Tukey implementation in AnalySeries. I use a Bartlett window
% with 90 % confidence for the phases

% load data
data=xlsread('blackman_tukey_phases.xlsx',[x1, ' vs ', x2]);
f=data(:,1);
ps.phase=data(:,5);
ps.pu=data(:,7);
ps.pl=data(:,6);
ps.coh=data(:,2);

% Get lags
obb=interp1(f,ps.phase,1/41); % interpolate to find phase at 1/41 kyr
ob = (obb*41)/(2*pi); % convert to kyrs
ecc=interp1(f,ps.phase,1/100);
ec = (ecc*100)/(2*pi);
prr=interp1(f,ps.phase,1/23);
pr = (prr*23)/(2*pi);

 % Get confidence levels
 obu=interp1(f,ps.pu,1/41); % interpolate to find phase at 1/41 kyr
ob_up = (obu*41)/(2*pi); % convert to kyrs
 obl=interp1(f,ps.pl,1/41); % interpolate to find phase at 1/41 kyr
ob_lo = (obl*41)/(2*pi); % convert to kyrs

 ecu=interp1(f,ps.pu,1/100); % interpolate to find phase at 1/41 kyr
ec_up = (ecu*100)/(2*pi); % convert to kyrs
 ecl=interp1(f,ps.pl,1/100); % interpolate to find phase at 1/41 kyr
ec_lo = (ecl*100)/(2*pi); % convert to kyrs

pru=interp1(f,ps.pu,1/23);
pr_up = (pru*23)/(2*pi);
prl=interp1(f,ps.pl,1/23);
pr_lo = (prl*23)/(2*pi);

% combine for output
ob = [ob, ob_lo, ob_up];
ec = [ec, ec_lo, ec_up];
pr = [pr, pr_lo, pr_up];

end