function [Vd,I,f_main,m_main] = f_k_spectrum(V,x1,x2,x3,t,dx,param,displ)

% f_k_spectrum calculates the joint wavenumber-frequency spectrum S(k,f)
% from two signals, given their separation dx and a sampling
% period dt, using continuous wavelet transforms.

%   V:      data from voltage probe
%	x1:     data from segment 1
%	x2:	    data from segment 2 
%   x3:     data from segment 3 (largest anode segment)
%	t:      time vector
%	dx:	    angular separation;  Default dx = 0.175 rad (10 deg)
%	param:  [fmin, fmax, # of bins in f, # of bins in k, sigma] ; 
%            Default [0 1e7 2000 2000 1]
%	displ:  set displ>0 (default) to get display

% If x1 and x2 are matrices, the wavelet transform is computed
% columnwise and the columns are concatenated

%	S:      scattering function S(k,f)
%	s:	    conditional scattering fct : s(k|f)=S(k,f)/S(f)
%	H:	    scattering function histogram H(k,f)  (= # of events)
%	coh:    cross-coherence between components 1 and 2 (same size as S)
%   growth: growth rate (imag part of k, same size as S)
%	kur:    kurtosis (same size as S)
%	freq:   frequency axis [Hz]
%	k:	    wavenumbers/2pi  [1/m]
%
%   ThDdW 11/94
%   Modified 07/01 NG
%   Updated 03/25 Ryan Przybocki

if nargin<6,	displ = 1;		end
if nargin<5,	param = [0 0 0 0 0];	end
if nargin<4,	dx = 1;			end
if nargin<3,	dt = 1;			end
if dx<=0,	dx = 1;			end

% If inputs are tables, convert to matrices:
if isa(V, 'table'), V=table2array(V); end
if isa(x1, 'table'), x1=table2array(x1); end
if isa(x2, 'table'), x2=table2array(x2); end
if isa(x3, 'table'), x3=table2array(x3); end
if isa(t, 'table'), t=table2array(t); end

[n,m] = size(x1);
[n2,m2] = size(x2);
if any([n2-n m2-m]), error('** x1 and x2 must have same size **'); end

% Voltage and Current averages:
Vd = mean(V); % Discharge voltage (V)
I_tot = x1 + x2 + x3;
I = mean(I_tot)*1000; % Discharge current (mA)

% Variables: 
dt = t(2)-t(1);
kNyq = 1/dx/2*(2*pi);	% Nyquist mode number (dx = 10*pi/180)
fNyq = 1/dt/2;			% Nyquist frequency
fupper = fNyq/2.001;	% maximum resolvable frequency
flower = 4/dt/n;		% minimum resolvable frequency
threshold = 100;		% no display if # of events < threshold (default 100)
nf = param(3);

np = length(param);
if np<5,	param = [param(:);zeros(4-np,1);1];		end
if param(5)==0,	sigma=1;	else	sigma   = param(5);	end
if param(4)==0,	nk = 41;	else	nk   = param(4);        end
if param(2)==0,	fmax = fupper;	else	fmax = param(2);	end
if param(1)==0,	fmin = flower;	else	fmin = param(1);	end

if fmax>fNyq/2
     disp(['** fmax cannot be larger than ',num2str(fupper,3),' **'])
     fmax = fupper;
end

if fmin<flower
     disp(['** fmin cannot be smaller than ',num2str(flower),' **']);
     fmin = flower;
end


if nf<=0
	nf = ceil(2*log(fmax/fmin)/log(2));
	if displ 
        disp(['nf = ',int2str(nf)]); 
    end
end


Fs = 1/dt;
k = linspace(-kNyq,kNyq,nk)';   % array of  k values
dk = 2*kNyq / nk;               % resolution in k

x1 = detrend(x1);
x2 = detrend(x2);

x1_mod = x1./max(abs(x1));
x2_mod = x2./max(abs(x2));

x1 = x1_mod;
x2 = x2_mod;

% Continuous wavelet transforms
% Default: analytic Morlet wavelet 'amor'

[wx1, f] = cwt(x1, 'amor', Fs, 'VoicesPerOctave', 48, 'ExtendSignal', false) ;         %'FrequencyLimits', [1e5 1e7]);  For lower freq data
[wx2, ~] = cwt(x2, 'amor', Fs, 'VoicesPerOctave', 48, 'ExtendSignal', false) ;        %'FrequencyLimits', [1e5 1e7]);

sz = size(wx1);
nf = sz(1); 

S = zeros(nf,nk);	
H = S;		
kur = S;
S11 = complex(S, 0);                
S12 = complex(S, 0);        
S22 = complex(S, 0);
fsp = zeros(nf,1);
neff = fsp;

scale = 1 ./(f*dt);

for i=1:nf
    
    ndec = round(min([1.5*scale(i) n/8]));	% # of pts in transients    
    
    wwx1 = wx1(i, ndec:n-ndec);
    wwx2 = wx2(i, ndec:n-ndec);
    
    w11 = wwx1.*conj(wwx1);
    w22 = wwx2.*conj(wwx2);
    w12 = wwx1.*conj(wwx2);

	kx = angle(w12)/dx;      % Divide by 2*pi to get wavenumber in cycles/m
	neff(i) = length(kx);			% eff # of data points <= n
	ampl = abs(w11+w22)/2/neff(i);		% avg power spectral density
	ampl2 = ampl.*ampl;			% 4th moment

	pos = floor(abs((kx+kNyq))/(dk+eps)+1);   % indices for k in S(k,f) 
	fsp(i) = sum(ampl)/neff(i);	% power spectrum in f
    
	for j=1:neff(i)
		in = pos(j);
		H(i,in) = H(i,in) + 1;
		S(i,in) = S(i,in) + ampl(j);
		kur(i,in) = kur(i,in) + ampl2(j);
        S11(i,in) = S11(i,in) + w11(j);
        S22(i,in) = S22(i,in) + w22(j);
        S12(i,in) = S12(i,in) + w12(j);
	end
end

% Identify main frequency-wavenumber feature:
[~, idx] = max(S(:));
[i_f, i_t] = ind2sub(size(S), idx);
f_main = f(i_f);
m_main = k(i_t);


kur = kur.*H./(S.*S+eps);                   % kurtosis
ksp = sum(S)'/nf;                           % wavenumber power spectrum
s = S./(sum(S')'*ones(1,nk)+eps);
coh = abs(S12)./sqrt(real(S11.*S22)+eps);	% cross-coherence
phi = angle(S12);                           % cross-phase
gr1 = (abs(S12)./(S11+eps)-1);      		% growth rate from 1 to 2
gr2 = (abs(S12)./(S22+eps)-1);              % growth rate from 2 to 1
growth = (gr1-gr2)/2;                       % avg growth rate


if displ
    
   % Normalization:
   %S = S./(sum(S,'all'));
   S = S./max(S,[],'all'); % Default: normalized to strongest feature
   
   % Plot spectrum
   pcolor(k,[0;f/1e6], log10([zeros(1,nk); S]));
   
   % Set colormap
   mycolormap = jet(256);
   mycolormap(1,:) = [1 1 1];
   colormap(mycolormap);
   hold on;
   
   h2=colorbar;
   h2.Label.String = 'Normalized Power Spectral Density';
   
   axis([-kNyq kNyq 0 10])
   
   plot([0 0],[0 10],'k--','LineWidth',1.5);
   
   shading('interp')	;
   xlabel('Mode Number', 'FontSize', 12); 	
   ylabel('Frequency (MHz)', 'FontSize', 12);
   box on;
   
end