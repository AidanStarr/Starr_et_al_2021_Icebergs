function varargout = nexcf(tx, x, ty, y,lag,h, verbose)
% NEXCF Cross-correlation estimates for non-equidistantly sampled time series.
%    C = NEXCF(TX, X, TY, Y) returns the cross-correlation function estimates
%    between the time series X and Y, for which the sampling times are given in
%    the column vectors TX and TY. The lags for which the Correlation functions
%    are estimated are given in the vector LAG=-T:DT:T. If this vector is not
%    provided, where T=RANGE([MAX(TX(1),TY(1)),MIN(TX(END),TY(END))])/2 and
%    DT=T/10. Unless specified, the default kernel width is 0.25*DT.
%
%    C = NEXCF(TX, X, TY, Y, LAG, H) computes the correlation using a kernel
%    width H, given in units of DT. H=0.25 amounts to a kernel width of 0.25*DT.
%
%    C = NEXCF(TX, X, TY, Y, LAG) returns the cross-correlation function at lags
%    given by vector LAG.
%
%    [C,LAGS] = NEXCF(...) returns the correlation function and the vector of
%    time lags (LAGS). Note that LAG refers to the absolute time
%    differences, not the vector index differences (e.g. give -5 0 5 not -1
%    0 1 if the average spacing of the time vector is 5).
%
%
%    EXAMPLE:
%       % We take two time series from stochastic processes X and Y, {tx,x} and
%       % {ty,y}, which are defined as follows:
%       % Sampling times
%       tx = (1:201)'; ty = (1:201)';
%       % Original signals
%       x = randn(201,1);
%       x(2:end) = 0.5*x(1:end-1)+randn(200,1);
%       y = randn(201,1);
%       y(2:end) = 0.7*x(1:end-1)+randn(200,1);
%
%       % We have, however, in reality, observed only 50% of the points in y,
%       % yielding a time series {ty_miss,y_miss}:
%       rem = randi(length(ty), floor(length(ty)*0.5),1);
%       ty_miss = ty; ty_miss(rem) = [];
%       y_miss = y; y_miss(rem) = [];
%
%       % We can now calculate the ACF for the time series {ty_miss,y_miss}:
%       % First we create a lag vector,
%       lag = -50:50;
%       % and we define the kernel width we want to employ:
%       H = 0.25;
%       % Now the ACF of {ty_miss;y_miss} is given by
%       Cy = nexcf(ty_miss, y_miss, ty_miss, y_miss, lag, H);
%       plot(lag,Cy), xlabel('Lag'), ylabel('Auto-correlation C_{y}')
%       % and the CCF, indicating in this example the lag and the
%       % magnitude of the coupling, is given by
%       Cxy = nexcf(tx, x, ty, y, lag, H);
%       plot(lag,Cxy), xlabel('Lag'), ylabel('Cross-correlation C_{xy}')
%
%    REFERENCE:
%       Ref.: K. Rehfeld, N. Marwan, J. Heitzig, J. Kurths: Comparison of
%       correlation analysis techniques for irregularly sampled time series,
%       Nonlinear Processes in Geophysics, 18(3), 389-404 (2011).
%
%    See also SIMILARITY,NEMI,ESF.

% (c) Potsdam Institute for Climate Impact Research 2011-2013
% Kira Rehfeld, Norbert Marwan
% ver: 1.4   last rev: 2013-09-02


%% output check
if nargout > 2
    error('Too many output arguments')
end

%% input check
error(nargchk(2, 7, nargin, 'struct'))

verb = 0;
useCuda = 1; %default
pd=0;

if nargin == 7
    verb = verbose;
end

% check if pdist2 is available
if ~exist('pdist2','file')
    % octave doesn't have pdist2, use modified functions at the end
    pd=1;
end

% check cuda capability
try
   if exist('gpuDeviceCount','file')
       gpu = gpuDeviceCount;
   else
       gpu = 0;
   end
catch
    gpu = 0;
end



if useCuda && gpu && verb
    disp('Wow! CUDA support available! You can enjoy accelerated calculation!')
end
if ~useCuda && gpu && verb
    disp('CUDA support available but not used.')
end

% length of the time series x and y
Nx = length(x); Ny = length(y);

% check the orientation of the vectors (vector multiplication in pdist...)
tx=tx(:);ty=ty(:);x=x(:);y=y(:);


% check for differences in length of time vector and observation vector
if (Nx ~= length(tx)) || (Ny ~=length(ty)),
    error('Vectors for sampling times and observations have to be of same length.'),
end

% check for overlap in time-axes
t_min = max(min(tx), min(ty));
t_max = min(max(tx), max(ty));

if (t_max-t_min)<eps
    error('The time-series do not overlap.')
end

% check if h is given
if nargin < 6
    h = 0.25;
end

% check if lag is given
% (if not, we use some default lag)
if nargin < 5
    T = range([t_min,t_max])/2;
    DT = T/10;
    %DT=10*mean(diff(T));
    lag(:,1) = -T:DT:T;
end

% check for ACF or CCF
if (length(tx)==length(ty)) && ~any(x-y)
    flag_acf = 1;
    showType = 'Calculate ACF';
else
    flag_acf = 0;
    showType = 'Calculate XCF';
end
if verb, disp(showType), end

%% data preprocessing
% make the time index dimensionless to make l and tx/ty comparable
if numel(lag) == 1
    dtlag = 1;
else
    %dtlag = abs(mean(diff(lag)));
    dtlag=max(mean(diff(tx)),mean(diff(ty)));
end
tx = tx/dtlag;
ty = ty/dtlag;
normlag = lag/dtlag;

%normalise the data
x = (x - mean(x)) / std(x);
y = (y - mean(y)) / std(y);




%% main part of the calculation

% Gauss kernel
b = @(t,h2,h2pi)(exp(-(t.^2)./h2) / h2pi);

if useCuda && gpu
    % initialise correlation vector
    C = gpuArray(zeros(length(lag),1)); % correlation vector
    
    % pairwise time-distances between all observation
    t_dist = gpuArray(pdist2(tx,ty, @(xi,yi)(xi-yi))); t_dist = t_dist(:);
    t_distx = gpuArray(pdist2(tx,tx, @(xi,yi)(xi-yi)));
    t_disty = gpuArray(pdist2(ty,ty, @(xi,yi)(xi-yi)));
    
    % prodcut matrix of data
    xy_dist = gpuArray(pdist2(x,y, @(xi,yi)(xi.*yi))); xy_dist = xy_dist(:);
    
else
    % initialise correlation vector
    C = zeros(length(lag),1); % correlation vector
    if ~pd
        % pairwise time-distances between all observation
        t_dist = pdist2(tx,ty, @(xi,yi)(xi-yi)); t_dist = t_dist(:);
        t_distx = pdist2(tx,tx, @(xi,yi)(xi-yi));
        t_disty = pdist2(ty,ty, @(xi,yi)(xi-yi));
        
        % prodcut matrix of data
        xy_dist = pdist2(x,y, @(xi,yi)(xi.*yi)); xy_dist = xy_dist(:);
    else
        % octave -> use my_pdist
        t_dist = pdist2Diff(tx,ty); t_dist = t_dist(:);
        t_distx = pdist2Diff(tx,tx);
        t_disty = pdist2Diff(ty,ty);
        
        % prodcut matrix of data
        xy_dist = pdist2Prod(x,y); xy_dist = xy_dist(:);
    end
end

% lower and upper limit of time differences
min_dist = min(t_dist(:));
max_dist = max(t_dist(:));

h2 = 2*h^2;
h2pi = sqrt(h2 * pi);

if verb, fprintf('lags... %03d %%',0), end
for i = 1:length(lag)
    if verb
        percent = round((i/length(lag))*100);
        fprintf('\b\b\b\b\b\b %03d %%',percent);
    end
    
    % check whether lag is outside possible time difference range
    if (normlag(i) <= min_dist) || (normlag(i) >= max_dist)
        C(i) = NaN;
        continue
    end
    
    % new distance matrix corresponding to considered lag
    t_dist_d = t_dist+normlag(i);
    
    % count the number of observations falling inside the kernel
    kernel_content = (t_dist_d <= 5*h) & (t_dist_d >= -5*h);
    % check whether this number is smaller than 5
    % (if yes, we do not believe in the correlation value)
    if sum(kernel_content) < 5
        C(i) = NaN;
        warning([sprintf('Lag %0.3f outside confident time region. Consider smaller range for lags.',lag(i))]);
        continue
    end
    
    % matrix of the weights (= apply kernel)
    WL = b(t_dist_d,h2,h2pi);
    Weight = sum(WL);
    % sum of the product of kernel with products x*y is the correlation
    C_ = sum(xy_dist.*WL) ./ Weight;
    C(i) = C_;
end


%% normalise the correlation values to a range [-1:1] using the variances of the time series
%matrix of the weights (= Gaussian kernel) for lag = 0
%WL = (exp(-((t_dist(:)).^2)./(2*h^2))) / sqrt(2*pi*h^2);
WLx = (exp(-((t_distx(:)).^2)./h2)) / h2pi;
WLy = (exp(-((t_disty(:)).^2)./h2)) / h2pi;

WeightX = sum(WLx);
WeightY = sum(WLy);

% product matrix of data x and y
if ~pd
    x_dist = pdist2(x,x, @(xi,yi)(xi.*yi));
    y_dist = pdist2(y,y, @(xi,yi)(xi.*yi));
else
    x_dist = pdist2Prod(x,x);
    y_dist = pdist2Prod(y,y);
end

% corrected variance as weighted mean, taking into account the irregularity
varX = sum( x_dist(:).*WLx) ./ WeightX;
varY = sum( y_dist(:).*WLy) ./ WeightY;

%% normalise the correlation estimates
C = C/(sqrt(varX) * sqrt(varY));

if any(C(:)>1),C=C/max(C(:));end
% numerically, it can happen that the variance is underestimated-> the
% correlation, as the ratio of covariance to variance overestimated. As a
% rule, we forbid the Correlation to be larger than 1.
%%
% transform gpuArray data to a standard Matlab array
if useCuda && gpu
    C = gather(C);
end

%% output of results
if nargout == 1
    varargout = {C};
elseif nargout == 2
    varargout(1) = {C};
    varargout(2) = {lag};
elseif nargout == 0
    C
end

%%%%%%% Define some functions that are not available in Octave (pdist2-
%%%%%%% replacements)
    function D = pdist2Prod(X,Y)
        n = length(X);
        m = length(Y);
        D = X * Y';
    end
    function D = pdist2Diff(X,Y)
        n = length(X);
        m = length(Y);
        D = ones(n,1) * Y' - X * ones(1,m) ;
    end
end