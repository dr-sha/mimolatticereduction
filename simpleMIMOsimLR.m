% -----------------------------------------------------
%
% -- MIMO simulator for lattice reduction 
% -- Author: Shahriar Shahabuddin
% -- Centre for Wireless Communications, Finland
% -- Email: shahriar.cwc@gmail.com
% -- website: https://sites.google.com/site/shahriarshahabuddin/matlab_simulator
% -- Original MIMO simulator by Christoph Studer
% -- http://www.csl.cornell.edu/~studer/software_mimo.html
%
% -----------------------------------------------------

function simpleMIMOsim(varargin)

  % -- set up default/custom parameters
  
  if isempty(varargin)
    
    disp('using default simulation settings and parameters...')
        
    % set default simulation parameters 
    par.simName = 'ERR_4x4_16QAM'; % simulation name (used for saving results)
    par.runId = 0; % simulation ID (used to reproduce results)
    par.MR = 4; % receive antennas 
    par.MT = 4; % transmit antennas (set not larger than MR!) 
    par.mod = '16QAM'; % modulation type: 'BPSK','QPSK','16QAM','64QAM'
    par.trials = 10000; % number of Monte-Carlo trials (transmissions)
    par.SNRdB_list = 0:4:40; % list of SNR [dB] values to be simulated
    par.detector = {'ZF','LLL_ZF','Fix_ZF'}; % define detector(s) to be simulated  
  else     
    disp('use custom simulation settings and parameters...')    
    par = varargin{1}; % only argument is par structure
    
  end

  % -- initialization
  
  % use runId random seed (enables reproducibility)
  rng(par.runId); 

  % set up Gray-mapped constellation alphabet (according to IEEE 802.11)
  switch (par.mod)
    case 'BPSK',
      par.symbols = [ -1 1 ];
    case 'QPSK', 
      par.symbols = [ -1-1i,-1+1i, ...
                      +1-1i,+1+1i ];
    case '16QAM',
      par.symbols = [ -3-3i,-3-1i,-3+3i,-3+1i, ...
                      -1-3i,-1-1i,-1+3i,-1+1i, ...
                      +3-3i,+3-1i,+3+3i,+3+1i, ...
                      +1-3i,+1-1i,+1+3i,+1+1i ];                 
    case '64QAM',
      par.symbols = [ -7-7i,-7-5i,-7-1i,-7-3i,-7+7i,-7+5i,-7+1i,-7+3i, ...
                      -5-7i,-5-5i,-5-1i,-5-3i,-5+7i,-5+5i,-5+1i,-5+3i, ...
                      -1-7i,-1-5i,-1-1i,-1-3i,-1+7i,-1+5i,-1+1i,-1+3i, ...
                      -3-7i,-3-5i,-3-1i,-3-3i,-3+7i,-3+5i,-3+1i,-3+3i, ...
                      +7-7i,+7-5i,+7-1i,+7-3i,+7+7i,+7+5i,+7+1i,+7+3i, ...
                      +5-7i,+5-5i,+5-1i,+5-3i,+5+7i,+5+5i,+5+1i,+5+3i, ...
                      +1-7i,+1-5i,+1-1i,+1-3i,+1+7i,+1+5i,+1+1i,+1+3i, ...
                      +3-7i,+3-5i,+3-1i,+3-3i,+3+7i,+3+5i,+3+1i,+3+3i ];
                         
  end

  % extract average symbol energy
  par.Es = mean(abs(par.symbols).^2); 
  
  % precompute bit labels
  par.Q = log2(length(par.symbols)); % number of bits per symbol
  par.bits = de2bi(0:length(par.symbols)-1,par.Q,'left-msb');

  % track simulation time
  time_elapsed = 0;
  
  % -- start simulation 
  
  % initialize result arrays (detector x SNR)
  res.VER = zeros(length(par.detector),length(par.SNRdB_list)); % vector error rate
  res.SER = zeros(length(par.detector),length(par.SNRdB_list)); % symbol error rate
  res.BER = zeros(length(par.detector),length(par.SNRdB_list)); % bit error rate

  % generate random bit stream (antenna x bit x trial)
  bits = randi([0 1],par.MT,par.Q,par.trials);

  % trials loop
  tic
  for t=1:par.trials
  
    % generate transmit symbol
    idx = bi2de(bits(:,:,t),'left-msb')+1;
    s = par.symbols(idx).';
    s_real = [real(s);imag(s)]; % converting s from complex to real.

    % generate iid Gaussian channel matrix & noise vector
    n = sqrt(0.5)*(randn(par.MR,1)+1i*randn(par.MR,1));
    n_real = [real(n);imag(n)]; % converting n from complex to real.

    
    H = sqrt(0.5)*(randn(par.MR,par.MT)+1i*randn(par.MR,par.MT));
    H_real = [real(H),-imag(H);imag(H),real(H)];  % converting H from complex to real.

    % Calculations needed for lattice reduction
    T = LLL(H_real,par.MT); 
    s_LR = inv(T)*s_real;
    H_LR = H_real*T;
       
    % Fixed LLL
    T_fixed = MLLL(H,par.MT);
    s_fixed = inv(T_fixed)*s;
    H_fixed = H*T_fixed;
     
    
    
    % transmit over noiseless channel (will be used later)
    x = H*s;
    x_real = [real(x);imag(x)];
    x_LR = H_LR*s_LR;
    x_fixed = H_fixed*s_fixed;

    
    
    
    % SNR loop
    for k=1:length(par.SNRdB_list)
      
      % compute noise variance (average SNR per receive antenna is: SNR=MT*Es/N0)
      N0 = par.MT*par.Es*10^(-par.SNRdB_list(k)/10);
      
      % transmit data over noisy channel
      y = x+sqrt(N0)*n;      
      y_real = x_real + sqrt(N0)*n_real;     
      y_LR = x_LR + sqrt(N0)*n_real;
      
      y_fixed = x_fixed + sqrt(N0)*n;

    
      % algorithm loop      
      for d=1:length(par.detector)
          
        switch (par.detector{d}) % select algorithms
          case 'ZF', % zero-forcing detection
            [idxhat,bithat] = ZF(par,H,y,x,n);
          case 'LLL_ZF', % zero-forcing detection
            [idxhat,bithat] = LLL_ZF(par,H_LR,y_LR,T); 
          case 'Fix_ZF', % zero-forcing detection  
            [idxhat,bithat] = Fix_ZF(par,H_fixed,y_fixed,T_fixed);
          otherwise,
            error('par.detector type not defined.')      
        end

        % -- compute error metrics
        err = (idx~=idxhat);
        res.VER(d,k) = res.VER(d,k) + any(err);
        res.SER(d,k) = res.SER(d,k) + sum(err)/par.MT;    
        res.BER(d,k) = res.BER(d,k) + sum(sum(bits(:,:,t)~=bithat))/(par.MT*par.Q);      
      
      end % algorithm loop
                 
    end % SNR loop    
    
    % keep track of simulation time    
    if toc>10
      time=toc;
      time_elapsed = time_elapsed + time;
      fprintf('estimated remaining simulation time: %3.0f min.\n',time_elapsed*(par.trials/t-1)/60);
      tic
    end      
  
  end % trials loop

  % normalize results
  res.VER = res.VER/par.trials;
  res.SER = res.SER/par.trials;
  res.BER = res.BER/par.trials;
  res.time_elapsed = time_elapsed;
  
  % -- save final results (par and res structure)
    
  save([ par.simName '_' num2str(par.runId) ],'par','res');    
    
  % -- show results (generates fairly nice Matlab plot) 
  
  marker_style = {'bo-','rs--','mv-.','kp:','g*-','c>--','yx:'};
  figure(1)
  for d=1:length(par.detector)
    if d==1
      semilogy(par.SNRdB_list,res.BER(d,:),marker_style{d},'LineWidth',2)
      hold on
    else
      semilogy(par.SNRdB_list,res.BER(d,:),marker_style{d},'LineWidth',2)
    end
  end
  hold off
  grid on
  xlabel('average SNR per receive antenna [dB]','FontSize',12)
  ylabel('bit error rate (BER)','FontSize',12)
  axis([min(par.SNRdB_list) max(par.SNRdB_list) 1e-4 1])
  legend(par.detector,'FontSize',12)
  set(gca,'FontSize',12)
  
end

% -- set of detector functions 
%% zero-forcing (ZF) detector
function [idxhat,bithat] = ZF(par,H,y,x,n)
  %xhat = H\y;    % inv(H)*y
  xhat = pinv(H)*y;
  [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
  bithat = par.bits(idxhat,:);
end


%% lattice reduced zero-forcing (LR_ZF) detector
function [idxhat,bithat] = LLL_ZF(par,H_LR,y_LR,T)
  xhat_real = H_LR\y_LR; % inv(H)*y
  xhat_quan = MyQuan_LR(xhat_real,2*par.MT,T,length(par.symbols)); 
  xhat =  xhat_quan(1:par.MT,:)+j*xhat_quan((par.MT+1):2*par.MT,:);
  [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
  bithat = par.bits(idxhat,:);
end


%% Fixed lattice reduced zero-forcing (Fix_ZF) detector, this is where MLLL is used
function [idxhat,bithat] = Fix_ZF(par,H_fixed,y_fixed,T_fixed)
  xhat_col = H_fixed\y_fixed; % inv(H)*y
  %Quantization for lattice reduction, quantization is performed according
  %to the constellation points of inv(T)*s
  xhat_quan = MyQuan_CLR(xhat_col,2*par.MT,T_fixed,length(par.symbols)); 
  %Converting from real to complex
  xhat =  xhat_quan(1:par.MT,:)+j*xhat_quan((par.MT+1):2*par.MT,:);
  %Quantization and generating the index of the bits
  [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
  %Finding the corresponding bits according to the indexing
  bithat = par.bits(idxhat,:);
end

