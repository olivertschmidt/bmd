function [L,P,f,idx] = cbmd(X,varargin)
% CBMD Cross-Bispectreal Mode Decomposition
%
% CBMD has the same inputs and outputs as BMD, but X must contain all three
% variables Q, R, and S to be cross-correlated. The second index (after
% time) must be the variable index, such that X(:,1,...) = Q,  X(:,2,...) =
% R, and  X(:,3,...) = S. If a long-time mean is provided via OPTS.mean,
% then the first index of OPTS.mean must also be the variable index.
%
%  References:
%  [1] Schmidt, O. T., Bispectral mode decomposition of nonlinear flows,
%      Nonlinear Dynamics, 2020
%
% O. T. Schmidt (oschmidt@ucsd.edu)
% Last revision: 08-Nov-2020


single_prec     = false;
if nargin==6
    opts = varargin{5};
    if ~isfield(opts,'regions')
        opts.regions    = [1 2];
    end
    if isfield(opts,'precision')
        if strncmpi(opts.precision,'single',6)
            single_prec = true;
        end
    end
else
    opts.regions = [1 2];
end

% get problem dimensions
dim     = size(X);
nt      = dim(1);
nVar    = dim(2);
if ~(nVar==3),  error('Second index of X must be 3 (3 variables must be provided for cross-correlation).'), end
nx      = prod(dim(3:end));

% get default spectral estimation parameters and options
[window,weight,nOvlp,dt,nDFT,nBlks] = parser(nt,nx,varargin{:});

% determine correction for FFT window gain
winWeight   = 1/mean(window);

% optimizers for x*Ax
if isfield(opts,'solver')
    switch opts.solver
        case {'HeWatson','simpleIteration','eig'}
        otherwise
            error('Unknown solver.')
    end
else
    opts.solver = 'HeWatson';
end

% standard tolerance
if ~isfield(opts,'tol')
    opts.tol    = 1e-8;
end

% use long-time mean if provided
if isfield(opts,'mean')
    x1_mean     = squeeze(opts.mean(1,:));
    x2_mean     = squeeze(opts.mean(2,:));
    x3_mean     = squeeze(opts.mean(3,:));
    mean_name   = 'provided long-time mean';    
else
    x1_mean     = mean(X(:,1,:),1);
    x2_mean     = mean(X(:,2,:),1);
    x3_mean     = mean(X(:,3,:),1);
    mean_name   = 'data mean';
end
x1_mean     = x1_mean(:);
x2_mean     = x2_mean(:);
x3_mean     = x3_mean(:);

disp(['Mean                      : ' mean_name]);

% obtain frequency axis
[f,nFreq,idx,f_idx,f1_idx,f2_idx,f3_idx] = faxes(nDFT,dt,opts);

nTriads = length(idx);

% loop over number of blocks and generate Fourier realizations
disp(' ')
disp('Calculating temporal DFT')
disp('------------------------------------')
Q1_hat = zeros(nFreq,nx,nBlks);
Q2_hat = zeros(nFreq,nx,nBlks);
Q3_hat = zeros(nFreq,nx,nBlks);
for iBlk    = 1:nBlks    
    % get time index for present block
    offset                  = min((iBlk-1)*(nDFT-nOvlp)+nDFT,nt)-nDFT;
    timeIdx                 = (1:nDFT) + offset;
    disp(['block ' num2str(iBlk) '/' num2str(nBlks) ' (' ...
        num2str(timeIdx(1)) ':' num2str(timeIdx(end)) ')'])
    % assemble block and subtract mean
    Q1_blk   = bsxfun(@minus,squeeze(X(timeIdx,1,:)),x1_mean.');
    Q2_blk   = bsxfun(@minus,squeeze(X(timeIdx,2,:)),x2_mean.');
    Q3_blk   = bsxfun(@minus,squeeze(X(timeIdx,3,:)),x3_mean.');
    
    % window and Fourier transform block
    Q1_blk                   = bsxfun(@times,Q1_blk,window);
    Q1_blk_hat               = winWeight/nDFT*fft(Q1_blk);
    Q1_blk_hat               = fftshift(Q1_blk_hat,1);
    Q2_blk                   = bsxfun(@times,Q2_blk,window);
    Q2_blk_hat               = winWeight/nDFT*fft(Q2_blk);
    Q2_blk_hat               = fftshift(Q2_blk_hat,1);
    Q3_blk                   = bsxfun(@times,Q3_blk,window);
    Q3_blk_hat               = winWeight/nDFT*fft(Q3_blk);
    Q3_blk_hat               = fftshift(Q3_blk_hat,1);
    
    Q1_hat(:,:,iBlk)         = Q1_blk_hat;
    Q2_hat(:,:,iBlk)         = Q2_blk_hat;
    Q3_hat(:,:,iBlk)         = Q3_blk_hat;
end
clear X x Q1_blk_hat Q2_b lk_hat Q3_blk_hat Q1_blk x1_mean Q2_blk x2_mean Q3_blk x3_mean

% loop over all triads and calculate CBMD
L      = nan(nFreq,nFreq);
disp(' ')
disp('Calculating BMD')
disp('------------------------------------')
P       = zeros(2,nTriads,nx);

if single_prec
    P        = single(P);
    Q1_hat   = single(Q1_hat);
    Q2_hat   = single(Q2_hat);
    Q3_hat   = single(Q3_hat);    
end

for i=1:nTriads
    disp(['(' num2str(f_idx(f1_idx(i))) ',' num2str(f_idx(f2_idx(i))) ',' num2str(f_idx(f3_idx(i))) ') (' num2str(i) '/' num2str(nTriads) ')'])
    Q_hat_f1            = squeeze(Q1_hat(f1_idx(i),:,:));
    Q_hat_f2            = squeeze(Q2_hat(f2_idx(i),:,:));
    Q_hat_f3            = squeeze(Q3_hat(f3_idx(i),:,:));
    
    Q_hat_f12           = Q_hat_f1.*Q_hat_f2;
    B                   = Q_hat_f3'*bsxfun(@times,Q_hat_f12 ,weight)/nBlks;    
    
    % optimizer for x*Ax
    switch opts.solver
        case {'HeWatson'}
            %  He & Watson's sophisticated iteration
            [r,a]  = HeWatson(B,opts.tol);
        case {'simpleit'}
            %  Watson's simple iteration
            a           = rand(nBlks,1) + 1i*rand(nBlks,1); % random initial guess
            [r,a]  = simpleIteration(B,a);
        otherwise
            error('Unknown solver.')
    end
    
    % i+j component
    Psi           = Q_hat_f3*a(:,1);
    Psi           = Psi/(Psi'*(Psi.*weight)); % normalize by inner product
    P(1,i,:)      = Psi;
    % i*j component
    Psi           = Q_hat_f12*a(:,1);
    Psi           = Psi/(Psi'*(Psi.*weight)); % normalize by inner product
    P(2,i,:)      = Psi;
    L(f1_idx(i),f2_idx(i)) ...
                  = r;
end
P   = reshape(P,[2 nTriads dim(3:end) 1]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [window,weight,nOvlp,dt,nDFT,nBlks] = parser(nt,nx,varargin)
% PARSER Parser for BMD parameters

% read input arguments from cell array
window = []; weight = []; nOvlp = []; dt = [];
nvarargin = length(varargin);

if nvarargin >= 1
    window = varargin{1};
    if nvarargin >= 2
        weight   = varargin{2};
        if nvarargin >= 3
            nOvlp   = varargin{3};
            if nvarargin >= 4
                dt      = varargin{4};
            end
        end
    end
end

window = window(:); weight = weight(:);

% check arguments and determine default spectral estimation parameters
% window size and type
if isempty(window)
    nDFT        = 2^floor(log2(nt/5)); if nDFT>256, nDFT=256; end
    window      = hammwin(nDFT);
    window_name = 'Hamming';
elseif length(window)==1
    nDFT        = window;
    window      = hammwin(window);
    window_name = 'Hamming';
elseif length(window) == 2^nextpow2(length(window))
    nDFT        = length(window);
    window_name = 'user specified';
else
    nDFT        = length(window);
    window_name = 'user specified';
end

% block overlap
if isempty(nOvlp)
    nOvlp = floor(nDFT/2);
elseif nOvlp > nDFT-1
    error('Overlap too large.')
end

% time step between consecutive snapshots
if isempty(dt)
    dt = 1/nDFT;
end

% inner product weight
if isempty(weight)
    weight      = ones(nx,1);
    weight_name = 'uniform';
elseif numel(weight) ~= nx
    error('Weights must have the same spatial dimensions as data.');
else
    weight_name = 'user specified';
end

% number of blocks
nBlks = floor((nt-nOvlp)/(nDFT-nOvlp));

% test feasibility
if nDFT < 4 || nBlks < 2
    error('Spectral estimation parameters not meaningful.');
end

% display parameter summary
disp(' ')
disp('BMD parameters')
disp('------------------------------------')
disp(['No. of snaphots per block : ' num2str(nDFT)])
disp(['Block overlap             : ' num2str(nOvlp)])
disp(['No. of blocks             : ' num2str(nBlks)])
disp(['Windowing fct. (time)     : ' window_name])
disp(['Weighting fct. (space)    : ' weight_name])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [window] = hammwin(N)
% HAMMWIN standard Hamming window of lenght N
window = 0.54-0.46*cos(2*pi*(0:N-1)/(N-1))';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w,z] = simpleIteration(A,z_0,tol)
% SIMPLEITERATION simple power iteration from He & Watson that is not
% guaranteed to find global optimum; small tolerences proposed in original
% paper tends to prevent convergence. 1e-8 works in most cases.
z       = z_0/sqrt(z_0'*z_0);
w       = Inf;
w_err   = Inf;
it_max  = 100;
it      = 0;
% tol   = 10*length(A)*eps*norm(A,1);
while w_err > tol
    
    w_old   = w;
    w       = z'*A*z;
    w_err   = abs(w-w_old);
    z       = w*A'*z + w'*A*z;
    z       = z/sqrt(z'*z);
    
    it      = it+1;
    if it>it_max, break, end
end
w       = z'*A*z;
% disp(['Watson's simple iteration required ' num2str(it) ' iterations.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w,z] = HeWatson(A,tol)
% HEWATSON 'An Algorithm' from He & Watson (1997) that is guranteed to find
% the global optimum upon convergence
N       = size(A,1);
normA   = norm(A,1);
z       = rand(N,1) + 1i*rand(N,1);
lb      = 0;
ub      = normA;
Z       = zeros(N);
I       = eye(N);
S       = [A Z; Z I];

it  = 0;
while (ub-lb)>tol || it==0
    it = it + 1;
    
    [w,z]   = simpleIteration(A,z,tol);
    lb      = max(lb,abs(w));
    alpha   = lb + tol;
    R       = [2*alpha*I -A'; I Z];
    [V,D]   = eig(R,S);
    D       = diag(D);
    
    ucirc   = abs(abs(D)-1) < (sqrt(eps)*normA);
    if sum(ucirc)==0
        break
    elseif mod(it,100)==0
        disp(['Optimizer did not converge in ' num2str(it) ' iterations! Trying new initial guess...']);
        z       = rand(N,1) + 1i*rand(N,1);
    elseif it>500
        disp(['Optimizer did not converge in ' num2str(it) ' iterations!']);
        break
    else
        idx = find(ucirc==1);
        z   = V(end-N+1:end,idx(1));
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f,nFreq,idx,f_idx,f1_idx,f2_idx,f3_idx] = faxes(nDFT,dt,opts)
% FAXES obtain frequency axes and indices
f_idx = (0:nDFT-1);
if mod(nDFT,2)==0
    f_idx(nDFT/2+1:end)     = f_idx(nDFT/2+1:end)-nDFT;
else
    f_idx((nDFT+1)/2+1:end) = f_idx((nDFT+1)/2+1:end)-nDFT;
end
f_idx   = fftshift(f_idx);
f       = f_idx/dt/nDFT;
fNyq_idx= -f_idx(1);
nFreq   = numel(f_idx);
if isfield(opts,'nfreq')
    f_idx_max   = opts.nfreq;
else
    f_idx_max   = fNyq_idx;
end

region  = nan(nFreq,nFreq);
idx     = nan(nFreq^2,1);
f1_idx  = idx;
f2_idx  = idx;
f3_idx  = idx;
count   = 0;
for i=1:nFreq
    for j=1:nFreq
        f1plus2  = f_idx(i)+f_idx(j);
        if abs(f1plus2)<fNyq_idx && abs(f_idx(i))<=f_idx_max && abs(f_idx(j))<=f_idx_max
            if sum(opts.regions==1)>0 && f_idx(i)>=0 && f_idx(j)>=0 && f_idx(i)>=f_idx(j)            % region #1
                region(i,j)     = 1;
            end
            if sum(opts.regions==2)>0 && f_idx(i)>=0 && f_idx(j)<=0 && f_idx(i)>=abs(f_idx(j))       % region #2
                region(i,j)     = 2;
            end
            if sum(opts.regions==3)>0 && f_idx(i)>=0 && f_idx(j)<=0 && f_idx(i)<=abs(f_idx(j))       % region #3
                region(i,j)     = 3;
            end
            if sum(opts.regions==4)>0 && f_idx(i)<=0 && f_idx(j)<=0 && f_idx(i)>=(f_idx(j))          % region #4
                region(i,j)     = 4;
            end
            if sum(opts.regions==5)>0 && f_idx(i)<=0 && f_idx(j)<=0 && f_idx(i)<=f_idx(j)            % region #5
                region(i,j)     = 5;
            end
            if sum(opts.regions==6)>0 && f_idx(i)<=0 && f_idx(j)>=0 && abs(f_idx(i))>=f_idx(j)       % region #6
                region(i,j)     = 6;
            end
            if sum(opts.regions==7)>0 && f_idx(i)<=0 && f_idx(j)>=0 && abs(f_idx(i))<=f_idx(j)       % region #7
                region(i,j)     = 7;
            end
            if sum(opts.regions==8)>0 && f_idx(i)>=0 && f_idx(j)>=0 && f_idx(i)<=f_idx(j)            % region #8
                region(i,j)     = 8;
            end
        end
        
        if ~isnan(region(i,j))
            count           = count + 1;
            idx(count)      = sub2ind([nFreq nFreq],i,j);
            f1_idx(count)   = i;
            f2_idx(count)   = j;
            f3_idx(count)   = find(f_idx==f1plus2);
        end
    end
end
f1_idx  = f1_idx(1:count);
f2_idx  = f2_idx(1:count);
f3_idx  = f3_idx(1:count);
idx     = idx(1:count);

end