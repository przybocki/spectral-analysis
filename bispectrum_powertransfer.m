function [ data ] = bispectrum_powertransfer( x , y, fs , delta , nfft , overlap , wind , pwr)
% Estimation of bispectrum products and using a direct fft-based approach.
% Nonlinear power transfer functions obtained with the method of Ritz et al. (1989)
%
% Inputs:
%   x       : signal 1;   y - signal 2;
%   fs      : sampling frequency
%   delta   : spatial separation between x and y (m)
%   nfft    : fft length 
%   overlap : percentage overlap [default = 50]
%   wind    : Type of window for tappering ('rectangular', 'hann' or 'kaiser')
%   pwr     : if true, compute coupling coeffs and power transfer
%
% Outputs: 
%   data    : data structure containing bispectra products
%
%
 
% Original routine:
% September 10, 2020
% Kévin Martins - kevin.martins@u-bordeaux.fr
%
%
% Updated:
% August 2024
% Ryan Przybocki - ryancp@stanford.edu
%
% Added second time series input (y), spatially separated downstream from x by distance delta
%     - computes y power spectrum E[Y Y*], cross power spectrum E[Y X*], and cross bispectrum [X1* X2* Y]
% Computes additional terms relevant for nonlinear power transfer using Ritz method
%     - Second order moments E[|X1 X2|^2]
%     - Linear L(f3) and quadratic Q(f1,f2) transfer functions, computed via Ritz method with Millionshchikov hypothesis
%     - Wave-wave coupling coefficient Lambda(f1,f2)
%     - Linear growth rate gamma(f3)
%     - Nonlinear power transfer function T(f1,f2)
%     - Cumulative nonlinear power transfer T_sum(f3), quantifies net quadratic energy transfer in/out of mode at f3
%
% References:
%       C. Ritz, E. Powers, and R. Bengston, Phys. Fluids B 1, 153 (1989).
%       J. Kim et al., Phys. Plasmas 3, 3998 (1996).
%       C. Ritz and E. Powers, Phys. D: Nonlinear Phenom. 20, 320 (1986).

% --------------------- Various parameters -------------------------
% Note: nfft is forced to be even, this will be useful to get frequencies centered around 0.

  lx = length(x); x = x(:);
  if (exist('nfft') ~= 1)            
      nfft = 128; 
  end
  if (exist('overlap') ~= 1)      
      overlap = 50;  
  end
  if (isempty(wind) == 1)            
      wind = 'rectangular'; 
  end
  
  overlap = min(99,max(overlap,0));

  nfft     = nfft - rem(nfft,2);
  eadvance = fix(nfft * overlap / 100);
  nadvance = nfft - eadvance;
  nblock   = fix((lx - eadvance) / nadvance) + 1; % +1 for not throwing any data out


  % ---------------------- Initialization ------------------------
    
  if (rem(nfft,2) == 0)
    freqs = [-nfft/2:nfft/2]'/nfft*fs;
  end

  data.info    = 'Bispectra products - the bicoherence is computed with Haubrich (1965) normalisation';
  data.f       = freqs;
  data.f_info  = 'Frequency [Hz] (two-sided)';
  data.df      = abs(freqs(2)-freqs(1));
  data.df_info = 'Frequency resolution [Hz]';
  data.P       = zeros(nfft+1,1);
  data.P_info  = 'X Power spectrum';
  data.B       = zeros(nfft+1,nfft+1);
  data.B_info  = 'X Power bispectrum';
  
  data.Py = zeros(nfft+1,1);
  data.Py_info = 'Y Power spectrum';
  
  data.C    = zeros(nfft+1,nfft+1);
  data.C_info = 'Cross power bispectrum X1*X2*Y';
  data.Pc   = zeros(nfft+1, 1);
  data.Pc_info = 'Cross power spectrum Y*conj(X)';
  
  data.B2 = zeros(nfft+1,nfft+1);
  data.B2_info = 'Second order moments |XiXj|^2';
  
  % Transfer function and energy transfer terms:
  data.L = zeros(nfft+1, 1);
  data.L_info = 'Linear transfer function';
  data.Q = zeros(nfft+1, nfft+1);
  data.Q_info = 'Quadratic transfer function';
  
  data.LambdaL = zeros(nfft+1,1);
  data.LambdaL_info = 'Linear coupling coefficient';
  data.LambdaQ = zeros(nfft+1, nfft+1);
  data.LambdaQ_info = 'Quadratic coupling coefficient';
  
  data.gamma = zeros(nfft+1, 1);
  data.gamma_info = 'Linear growth rate';
  
  data.T = zeros(nfft+1, nfft+1);
  data.T_info = 'Nonlinear Power transfer function';
  data.Tsum = zeros(nfft+1, 1);
  data.Tsum_info = 'Summed nonlinear power transfer';

    
  format long g
  
  % ---------------------- Compute FFT ----------------------------

  % Initialization
  A = zeros(nfft+1,nblock);     % Fourier coefficients on x , for each block
  Ay = zeros(nfft+1,nblock);    % Fourier coefficients on y, for each block
  nmid = (nfft)/2 + 1;          % Middle frequency (f = 0)
  locseg = [1:nfft]';           % Indices for first block
    
  % Computing FFT (loop over blocks)

  for kk = 1:nblock-1
    % Preparing block kk timeseries
    % Treatment varies with the window
    % For the rectangular window, we force a certain continuity between blocks
    xseg = x(locseg);
    xseg = detrend(xseg);          % Detrend
    xseg = (xseg(:) - mean(xseg)); % De-mean
    
    yseg = y(locseg);
    yseg = detrend(yseg);          % Detrend
    yseg = (yseg(:) - mean(yseg)); % De-mean
    
    switch wind
      case 'hann'
        ww = window(@hann,nfft); normFactor = mean(ww.^2);
        xseg = xseg.*ww / sqrt(normFactor);
        yseg = yseg.*ww / sqrt(normFactor);
      case 'kaiser'
        ww = window(@kaiser,nfft,3.5); normFactor = mean(ww.^2);
        xseg = xseg.*ww / sqrt(normFactor);
        yseg = yseg.*ww / sqrt(normFactor);
      case 'rectangular'
        % Trying to make it periodic
        count = 0;
        while abs(xseg(end)-xseg(1)) > 0.2*nanstd(xseg)
          % Updating locseg
          if kk == 1
            locseg = locseg + 1;
          else
            locseg = locseg - 1;
          end
          % Updating xseg
          xseg = x(locseg);
          count = count + 1;
        end
        if count > 1
          xseg = detrend(xseg);          % Detrend
          xseg = (xseg(:) - mean(xseg)); % De-mean
        end
        
        count = 0;
        while abs(yseg(end)-yseg(1)) > 0.2*nanstd(yseg)
          % Updating locseg
          if kk == 1
            locseg = locseg + 1;
          else
            locseg = locseg - 1;
          end
          % Updating yseg
          yseg = y(locseg);
          count = count + 1;
        end
        if count > 1
          yseg = detrend(yseg);          % Detrend
          yseg = (yseg(:) - mean(yseg)); % De-mean
        end
        
        % Final windowing
        ww = window(@rectwin,nfft); normFactor = mean(ww.^2);
        xseg = xseg.*ww / sqrt(normFactor);
        yseg = yseg.*ww / sqrt(normFactor);
    end
    
    % FFT
    A_loc   = fft( xseg , nfft )/nfft;
    Ay_loc  = fft( yseg , nfft )/nfft;
    
    A(:,kk) = [ A_loc(nmid:nfft,:) ; A_loc(1:nmid,:) ]; % FFTshift x
    A(nmid,kk) = 0;
    
    Ay(:,kk) = [ Ay_loc(nmid:nfft,:) ; Ay_loc(1:nmid,:) ]; % FFTshift y
    Ay(nmid,kk) = 0;
    
    % Indices for next block
    locseg = locseg + nadvance;
  end
  
  % Last block, we are not throwing any data out
  for kk = nblock:nblock
    % Preparing block kk timeseries
    % Treatment varies with the window
    % For the rectangular window, we force a certain continuity between blocks
    locseg = [length(x)-nfft+1:length(x)]';
    xseg = x(locseg);
    xseg = detrend(xseg);          % Detrend
    xseg = (xseg(:) - mean(xseg)); % De-mean
    
    yseg = y(locseg);
    yseg = detrend(yseg);          % Detrend
    yseg = (yseg(:) - mean(yseg)); % De-mean
    
    switch wind
      case 'hann'
        ww = window(@hann,nfft); normFactor = mean(ww.^2);
        xseg = xseg.*ww / sqrt(normFactor);
        yseg = yseg.*ww / sqrt(normFactor);
      case 'kaiser'
        ww = window(@kaiser,nfft,3.5); normFactor = mean(ww.^2);
        xseg = xseg.*ww / sqrt(normFactor);
        yseg = yseg.*ww / sqrt(normFactor);
      case 'rectangular'
        % Trying to make it periodic
        count = 0;
        while abs(xseg(end)-xseg(1)) > 0.2*nanstd(xseg)
          % Updating locseg
          locseg = locseg - 1;
          % Updating xseg
          xseg = x(locseg);
          count = count + 1;
        end
        if count > 1
          xseg = detrend(xseg);          % Detrend
          xseg = (xseg(:) - mean(xseg)); % De-mean
        end
        
        count = 0;
        while abs(yseg(end)-yseg(1)) > 0.2*nanstd(yseg)
          % Updating locseg
          locseg = locseg - 1;
          % Updating yseg
          yseg = y(locseg);
          count = count + 1;
        end
        if count > 1
          yseg = detrend(yseg);          % Detrend
          yseg = (yseg(:) - mean(yseg)); % De-mean
        end
        
        % Final windowing
        ww = window(@rectwin,nfft); normFactor = mean(ww.^2);
        xseg = xseg.*ww / sqrt(normFactor);
        yseg = yseg.*ww / sqrt(normFactor);
    end
    
    % FFT
    A_loc   = fft( xseg , nfft )/nfft;
    Ay_loc   = fft( yseg , nfft )/nfft;
    
    A(:,kk) = [ A_loc(nmid:nfft,:) ; A_loc(1:nmid,:) ]; % FFTshift x
    A(nmid,kk) = 0;
    
    Ay(:,kk) = [ Ay_loc(nmid:nfft,:) ; Ay_loc(1:nmid,:) ]; % FFTshift y
    Ay(nmid,kk) = 0;
    % Indices for next block
    locseg = locseg + nadvance;
  end
    
  
  % -------------- Computing bispectrum products -------------------
    
  
  % Dealing with f1 + f2 = f3 indices
  [ ifr1 , ifr2 ] = meshgrid( 1:nfft+1 , 1:nfft+1 );
  ifr3 = nmid + (ifr1-nmid) + (ifr2-nmid);
  ifm3val = ((ifr3 >= 1) & (ifr3 <= nfft+1)); ifr3(~ifm3val) = 1;

  data.ifr3 = ifr3;
  % freq(ifr1(1, a)) + freq(ifr2(b, 1)) = fa+fb = freqs(ifr3(a,b))
  
  % Accumulating triple products (loop over blocks)
  for kk = 1:nblock
    % Block kk FFT
    A_loc  = A(:,kk);
    CA_loc = conj(A(:,kk));
    
    Ay_loc = Ay(:,kk);
    CAy_loc = conj(Ay(:,kk));
    
    % Compute bispectrum and PSD
    data.B = data.B + A_loc(ifr1) .* A_loc(ifr2) .* (CA_loc(ifr3)); % Bispectrum
    data.P = data.P + abs(A_loc.^2); % input (X) power spectrum
    data.Py = data.Py + abs(Ay_loc.^2); % output (Y) power spectrum
    
    % Compute cross bispectrum, cross power spectrum, 2nd order moments:
    data.C = data.C + CA_loc(ifr1) .* CA_loc(ifr2) .* Ay_loc(ifr3);     % cross bispectrum: conj(X1)conj(X2)(Y)
    data.Pc = data.Pc + Ay_loc.*CA_loc;                                 % Cross power spectrum = Y*conj(X)
    data.B2 = data.B2 + abs(A_loc(ifr1).*A_loc(ifr2)).^2;               % Squared 2nd order moment |X1 X2|^2
  end

  % Expected values
  data.B = data.B / nblock; data.B(~ifm3val) = 0;
  data.P = data.P / nblock;
  data.Py = data.Py / nblock;
  
  data.C = data.C / nblock; data.C(~ifm3val) = 0;
  data.Pc = data.Pc / nblock;
  data.B2 = data.B2 / nblock;   data.B2(~ifm3val) = 0;
  
  
  % Computing bicoherence and biphase, cross bicoherence

   data.Bic = abs(data.B) ./ sqrt( data.P(ifr1) .* data.P(ifr2) .* data.P(ifr3) ); % Default Normalization 
   % data.Bic = abs(data.B) ./ sqrt((data.B2 .* data.P(ifr3))) ; % Cauchy-Schwarz normalization
   data.Bic_info  = 'Auto-Bicoherence [-]';    
    
   data.cBic = abs(data.C) ./ sqrt( data.P(ifr1) .* data.P(ifr2) .* data.Py(ifr3) );
   % data.cBic = abs(data.C) ./ sqrt((data.B2 .* data.Py(ifr3))) ; % Cauchy-Schwarz normalization
   data.cBic_info = 'Cross-Bicoherence';
    
   data.Bip = atan2(imag(data.B),real(data.B));
   data.Bip_info  = 'Biphase [-]';

  
  % Linear and Quadratic coupling coefficients:
  
 
  if pwr
      format long; 
      bsum = zeros(nfft+1, 1);
      asum = zeros(nfft+1, 1);
      
      bsum_term = ( data.B .* data.C ) ./ data.B2;  % numerator sum
      asum_term = ( abs(data.B).^2 ) ./ data.B2;    % denominator sum

      bsum_term(isnan(bsum_term)) = 0;  % Remove NaN's for summation
      asum_term(isnan(asum_term)) = 0;
      
      % For each p, sum over all p1, p2 indices satisfying f1 + f2 = f
      % Factor of 0.5 accounts for double counting p1, p2
      
      for p = 1:nfft+1
          
          % If p is even (frequency index odd), sum over p1>=p2 is just 0.5 * sum(all p1, p2)
          % If p is odd, sum over p1>=p2 is 0.5 * sum(all p1, p2) + 0.5 * (p1, p2) term
          
          bsum(p) = bsum(p) + 0.5 * sum(bsum_term(ifr3 == p));
          asum(p) = asum(p) + 0.5 * sum(asum_term(ifr3 == p));
          
          if mod(p,2)
              p0 = p/2 + nfft/4 + 1/2;
              bsum(p) = bsum(p) + 0.5 * bsum_term(p0,p0);
              asum(p) = asum(p) + 0.5 * asum_term(p0,p0);
          end
          
          
      end

      clear p;

      epsilon = 0; % Can either add in small epsilon to denominator or remove NaN's later
      
      Exp_phase = data.Pc(ifr3) ./ abs(data.Pc(ifr3));   % Exp(i*delta_theta) from Ritz et al. (1989) eq 12
      
      data.L = ( data.Pc - bsum ) ./ (data.P - asum + epsilon);   % Linear transfer function
      % data.gamma = (abs(data.L).^2 - 1) ./ delta;   % Alt linear growth rate given by Kim et al. (1996)
      data.L(isnan(data.L)) = 0; 

      
      data.Lpp = zeros(size(nfft+1, nfft+1));
      data.Lpp = data.L(ifr3);
      
      % Construct linear coupling coefficient LambdaL(f)
      for p = 1:nfft+1;
          exp_phase_p = data.Pc(p) ./ abs(data.Pc(p));
          phase_p = log(exp_phase_p);
          data.LambdaL(p) = 1/delta * (data.L(p) / exp_phase_p + phase_p - 1);
          data.gamma(p) = 1/delta .* real(data.L(p) / exp_phase_p - 1);    % Linear growth rate in spectral power equation
      end
      data.gamma(isnan(data.gamma)) = 0;
      
      % Equivalent definitions of phase angle:
      data.Theta = angle(data.Pc);
      e_phase = data.Pc ./ abs(data.Pc);
      data.Theta1 = -1i*log(e_phase);
      
      % Quadratic transfer function
      Bstar = conj(data.B);
      data.Q = ( data.C - Bstar .* data.Lpp ) ./ data.B2; 
      data.Q(isnan(data.Q)) = 0;
      
      % Quadratic coupling coefficient LambdaQ
      data.LambdaQ = data.Q ./ Exp_phase ./ delta;
      data.LambdaQ(isnan(data.LambdaQ)) = 0;  % Remove NaN's for summation

      % Quadratic power transfer function
      data.T = real( data.LambdaQ .* data.B );
      
      % Sum over all f1, f2 contributing to power transfer into f3
      for p = 1:nfft+1
          data.Tsum(p) = sum(data.T(ifr3 == p));
      end
      
      clear p;
      
      
  end
  

  