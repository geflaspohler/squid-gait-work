prhdir = '/Users/yorkpatty/Documents/MATLAB_SQUID/SquidTest107/prh'
caldir = '/Users/yorkpatty/Documents/MATLAB_SQUID/SquidTest107/cal'
recdir = '/Users/yorkpatty/Documents/MATLAB_SQUID/SquidTest107/raw'
settagpath('cal', caldir, 'prh', prhdir);

prefix = 'SquidTest107';
df = 1;
X = d3readswv(recdir, prefix, df);
[CAL,D] = d3deployment(recdir, prefix, 'SquidTest107');
CAL = d3findcal('d401') ;
%%
[p,CAL]=d3calpressure(X,CAL,'full');
[A,CAL,fs] = d3calacc(X,CAL,'full') ;
[M,CAL] = d3calmag(X,CAL,'full') ;
d3savecal('SquidTest107','CAL',CAL)

%% If want to continue to estimate the orientation
TH = 0.2;
surface = 0.1;
METHOD = 2;
[PRH, T] = prhpredictor(p,A,fs,TH, METHOD, [], surface);

%% If want to continue to save all the data with p,r,h estimation
OTAB = [0 0 0 180 0];
OTAB = ((OTAB)./360)*2*pi;
[Aw,Mw] = tag2whale(A,M,OTAB,fs);
% All_Plot(p, Aw, Mw, fs)
d3savecal(prefix,'OTAB',OTAB);
% d3makeprhfile(recdir,prefix,deploy_name,df)
d3makeprhfile(recdir,prefix, prefix,df);

