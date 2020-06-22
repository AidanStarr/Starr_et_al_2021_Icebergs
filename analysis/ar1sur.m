function [Xsur]=ar1sur(tx,x,Nsur)
% AR1SUR computes autocorrelated irregular surrogate time series.
% [Xsur]=AR1SUR(tx,x,Nsur) calculates a matrix of NSur surrogate time series, Xsur, that fits the
% time axis of the  Sx with time axis Tx of a given Vector X(tx). 
% The persistence time (or: autocorrelation) is estimated from a least
% squares fit of AR1-processes to the data (standard) or, if the
% Optimization toolbox is not available, using the autocorrelation function
% of the data (e.g. in Octave).
% 
%
%    REFERENCE:
%
%       ESF: K. Rehfeld and J. Kurths: Similarity measures for irregular
%       and age uncertain time series, submitted to Clim. Past. Discuss,
%       2013.
%
% see also: NEXCF,GMI,SIMILARITY,CAR,TAR1
%
% (c) Potsdam Institute for Climate Impact Research 2013
% (c) Alfred Wegener Institute for Polar and Marine Research 2014
% Kira Rehfeld
% ver: 1.1  last rev: 2014-01-25

acf=0;

if length(tx) ~= length(x); error('tx and x not same length');
end;
%dtau=median(diff(tx));
    % use AR1-coefficient from autocorrelation function
    acf=true;
    Tau=NaN;



if Tau<min(diff(tx))||~isreal(Tau)||acf==true,
    % use ACF coefficients to estimate persistence
    lag=[1:2]*mean(diff(tx));
    
     Cxx=similarity(tx,x,tx,x,lag,'gXCF');
     if any(isnan(Cxx))% choose smaller decorrelation scale
     	Cxx=similarity(tx,x,tx,x,1:2,'gXCF');
     	end
     	
     	
     EstTau=-[1,2].*mean(diff(lag))./log(Cxx([1,2]))';
     Tau=mean(real(EstTau));
     est_phi=exp(-mean(diff(tx))/(Tau));

end

Xsur=zeros(length(tx),Nsur);
%stream=RandStream('mlfg6331_64','RandnAlg','Polar');
for k=1:Nsur
    x=randn(size(tx));
    % adjust standard deviation of the noise such that the
    % standard deviation of the resulting time series is 1:
    Eps=randn(length(tx),1);
    Eps=(Eps-mean(Eps))./std(Eps);
    Xsuri=Eps;
    % possible correction
    %%KK=sqrt((1-exp(-2*diff(tx)./TauKern))./(2./TauKern));
    KK=1;
    for i=2:length(tx)
        Xsuri(i)=exp(-(tx(i)-tx(i-1))./Tau).*Xsuri(i-1) +KK.*Eps(i);
    end
    
    Xsuri=(Xsuri-mean(Xsuri))./std(Xsuri) *std(x);
    
    Xsur(:,k)=Xsuri;
    clear Xsuri
    
end
