% Default parameters for FISTA, MFISTA, MFISTA_TV
%
% copyright 2013 Hatef Monajemi (monajemi@stanford.edu)
function defopts = astro_setparams(opts)

%-------------------------------------------
%   BEGIN DEFAULT VALUES, ALL LOWER CASE
%-------------------------------------------
% ASTROSOLVE General params 

defopts.remblack=0;
defopts.blurthresh=0.01;  %threshold for bluring mask size
defopts.bluralpha=1.449;  %alpha in psf
defopts.blurR0=2.2364;    %r0 in psf
defopts.wavelet=[1 8];	  %wavelet type
defopts.nitr=1500;	      %Number of iterations
defopts.ps=0;			  %ps=1 if there are point sources
defopts.psclean=0;        %clean connected components after point sources estimate
defopts.smoothf=0;        %smooth profile 
defopts.showplot=1;       %show plot at the end
defopts.prevfw=[];        %give inirial values
defopts.sshrink=2;        %Add king function (set to 0 if just wavelets)
defopts.beta0=[];         %intercept
defopts.itersol=0;
defopts.forcepos=2;       
defopts.prevPS=[];
defopts.parallel=0;
defopts.center=0;         %uncertainty of the center (in fraction)
defopts.bootwindow=0    %uncertainty=m/64;

%----------------------------------------
%        END DEFAULT OPTS
%----------------------------------------

% return if no user opts
if nargin == 0 || isempty(opts) 
	return
end

% list of valid field names
vfields = fieldnames( defopts );

% Grab valid fields from user's opts
for i = 1:length(vfields)
    field = vfields{i};
    if isfield( opts, field );
        defopts.(field) = opts.(field);
    end
end
