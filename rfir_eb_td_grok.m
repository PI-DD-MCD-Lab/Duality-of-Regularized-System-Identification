function  [th, obj, hp] = rfir_eb_td_grok(kernel,data, pre_ind, configs)
% RFIR_EB_TD Estimates impulse response using empirical Bayes with a Gaussian process prior.
% Inputs:
%   kernel  - String specifying the kernel type (only 'DC' supported currently)
%   data    - iddata object containing input and output signals (single experiment)
%   pre_ind - Numeric vector/scalar of prediction indices (converted to column if row)
%   configs     - Cell array with configuration options (e.g., {'frf'}, {'eir'}, {'cut', n})
% Outputs:
%   th      - Estimated impulse response coefficients
%   obj     - Optimized objective value (negative log-likelihood)
%   hp      - Optimized hyperparameters
% @ Xiaozhu Fang 2025/11/17

narginchk(4,4);
validateattributes(data, {'iddata'}, {});
validateattributes(pre_ind, {'numeric'}, {'vector', 'nonnegative'});
if  isrow(pre_ind)
    pre_ind=pre_ind';
end
kernel = upper(kernel);
validatestring(kernel, {'DC'});
validateattributes(configs, {'cell'}, {});

% if ~realdata(data)
%     error('Only real-valued signals supported.');
% end

% main body
[hpini, Phi, Yt, lb, ub] = initialize_parameters(data, pre_ind, kernel, configs);
objFunc = @(x)objective_function(x, Phi, Yt, pre_ind, kernel, lb, ub, configs);
[hp, ~, ~] = fminsearch(objFunc, hpini);
% Evaluate best solution
[obj1, ~, th1] = objective_function(hpini, Phi, Yt, pre_ind, kernel, lb, ub, configs);
[obj2, ~, th2] = objective_function(hp, Phi, Yt, pre_ind, kernel, lb, ub, configs);

if obj1 > obj2
    th = th2;
    obj = obj2;
else
    th = th1;
    obj = obj1;
end

function [hpini, Phi,Yt, lb, ub] = initialize_parameters(z, pre_ind, kernel, configs)
% Inputs:
    %   z       - iddata object
    %   pre_ind - Prediction indices
    %   kernel  - Kernel type ('DC')
    %   configs    - Configuration cell array
    % Outputs:
    %   hpini   - Initial hyperparameters
    %   Phi     - Input signal (time domain)
    %   Yt      - Output signal (time domain)
    %   lb, ub  - Hyperparameter bounds
y = pvget(z,'OutputData');
u = pvget(z,'InputData');
y = y{1}; % Single experiment
u = u{1};
switch configs{1}
    case 'frf' 
%         Phi = circulant(flip(u), 1);
        Phi = circulant([flip(u)]',1);
        Yt = y;
    case "eir"
        N = length(y);
        Phi = flipud(hankel(flip([0; u(1:end-1)]), zeros(length(u), 1)))';
        Phi = Phi(1:N, :)';
        Yt = y;
    case "eir-no-t"
        N = length(y);
        Phi = flipud(hankel(flip([0; u(1:end-1)]), zeros(length(u), 1)))';
        Phi = Phi(1:N, :)';
        Yt = y;
    case "cut"
        n = configs{2};
        Yt = y((n+1):end);
        N = length(Yt);
        Phi_full = flipud(hankel(flip([0; u(1:end-1)]), zeros(length(u), 1)))';
        Phi = Phi_full(1:n, (n+1):(n+N))';
    otherwise
        error(['Unsupported tmp{1} option: ' configs{1}]);
end


% Define bounds for DC kernel
if sum(strcmp(kernel, {'DC'}))
    lb = -inf(5,1);
    ub = inf(5,1);
    lbga = -10 * ones(5,1);
    ubga = 10 * ones(5,1);
end

% Optimize initial hyperparameters using MultiStart
ratio = 5;
if sum(strcmp(kernel, {'TC','DC'})) %      c_sign is useless for these kernels
    objFunc =@(x) objective_function(x, Phi, Yt, pre_ind, kernel, lb, ub, configs);
    ms =  MultiStart('FunctionTolerance',1e-10,'UseParallel',false);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon', 'x0', rand(size(lbga)), 'objective', objFunc, 'lb', lbga, 'ub', ubga);
    [hpini ,~,~] =  run(ms, problem, ratio*length(lbga));
end

function [obj, condNb, xest] =  objective_function(hyper, Phi, Yt, pre_ind, kernel, lb, ub, configs)
% NGLGLKLHD_SIMP Computes negative log-likelihood and impulse response
    % Inputs:
    %   hyper   - Hyperparameters [log(alpha), log(beta), log(c1), log(c2), log(sigma2)]
    %   Ut      - Input signal (time domain)
    %   Yt      - Output signal (time domain)
    %   pre_ind - Prediction indices
    %   kernel  - Kernel type 
    %   lb, ub  - Hyperparameter bounds
    %   configs     - Configuration cell array
    % Outputs:
    %   obj     - Negative log-likelihood
    %   condNb  - Conditional number of covariance matrix
    %   xest    - Estimated impulse response
warning('off', 'MATLAB:nearlySingularMatrix');
chol_flag = 1;
cond_flag = 1;
ini_flag = 1;
N= length(Yt);

ini_flag_di = sum(hyper<=ub) + sum(hyper>=lb);
if ini_flag_di ~= 2*length(hyper)
    ini_flag = 0;
end
c1=exp(hyper(end-2));
c2=exp(hyper(end-1));
sigma2= exp(hyper(end));
switch kernel
    case 'DC' 
        vec=1:N;
        alpha=exp(hyper(1));
        beta=exp(hyper(2));
        switch configs{1}
            case  "frf"
                Pi= kdc_gp(vec', vec, beta, alpha, N);

                if nargout >1
                    Pipre=kdc_ggp( pre_ind,vec, beta, alpha, N) ;
                    if sum(sum(isnan(Pipre)))
                        Pipre = zeros(length(pre_ind), N);
                        for n= 0 :N:(100*N)
                            Pipre= Pipre+ kdc(pre_ind, vec + n, beta, alpha);
                        end
                    end
                end
                Wf = (0:(N-1))*2*pi/N;
                Ti_temp = 1./(1-exp(-beta-1j*(Wf'-Wf))).* (1-exp(-2*alpha - beta-1j*(Wf'-Wf)))./(1-exp(-alpha- beta/2-1j*(Wf')))./(1-exp(-alpha- beta/2+1j*(Wf)));
                Ti = fft(ifft(Ti_temp,[],1),[],2);
                Ti = real(Ti);
                % yvec=1:N;
                % Ti= exp(-alpha*abs(yvec'-yvec)) .* exp(-beta*(yvec'+yvec)/2) ;
            case "eir"
                Pi=exp(-beta/2*(vec'+vec)-alpha*abs(vec'-vec));
                if nargout >1
                    Pipre = exp(-beta/2*(pre_ind+vec)-alpha*abs(pre_ind-vec));
                end
                yvec=1:N;
                Ti= exp(-alpha*abs(yvec'-yvec)) .* exp(-beta*(yvec'+yvec)/2) ;

                if strcmp(configs{2},"off")
                else
                    Pi(configs{2}:end,:)=0;
                    Pi(:, configs{2}:end)=0;
                    Ti(configs{2}:end,:)=0;
                    Ti(:, configs{2}:end)=0;
                end

            case"eir-no-t"
                Pi=exp(-beta/2*(vec'+vec)-alpha*abs(vec'-vec));
                c2=0;
                if nargout >1
                    Pipre = exp(-beta/2*(pre_ind+vec)-alpha*abs(pre_ind-vec));
                end
                yvec=1:N;
                Ti= exp(-alpha*abs(yvec'-yvec)) .* exp(-beta*(yvec'+yvec)/2) ;
            case "cut"
                vec=1:configs{2};
                c2=0;
                Pi=exp(-beta/2*(vec'+vec)-alpha*abs(vec'-vec));
                if nargout >1
                    Pipre = exp(-beta/2*(pre_ind+vec)-alpha*abs(pre_ind-vec));
                end
                yvec=1:N;
                Ti= exp(-alpha*abs(yvec'-yvec)) .* exp(-beta*(yvec'+yvec)/2) ;
        end

        
end

if ini_flag == 0
    obj = inf;
    return;
end
if  cond_flag ~= 1 || chol_flag ~=1
    obj = 1/eps;
else
    M = c1 * Phi * Pi * Phi' + c2 * Ti + sigma2 * eye(N);
    try
        cholFactor = chol(M, 'lower');
    catch
        try
            cholFactor = chol((M + M') / 2, 'lower');
        catch
            obj = 1 / eps;
            xest = zeros(length(pre_ind), 1);
            condNb = nan;
            return;
        end
    end
    a = cholFactor' \ (cholFactor \ Yt);
    obj = real( Yt' * a + 2 * sum(log(diag(cholFactor))));
end
% Compute optional outputs
if nargout > 1
    xest = c1 * Pipre * Phi' * a;
    condNb = cond(M);
end

function C = circulant(vec,direction)
if (nargin<1) || (nargin > 2)
  error('circulant takes only one or two input arguments')
end
if (nargin < 2) || isempty(direction)
  direction = 1;
elseif ~ismember(direction,[1,-1])
  error('direction must be either +1 or -1 if it is supplied')
end
if ~isvector(vec)
  error('vec must be a vector')
elseif length(vec) == 1
  C = vec;
  return
end
n = length(vec);
n1 = n-1;
if direction == -1
  C = repmat(0:n1,n,1);
  C = vec(mod(C+C',n)+1);
else
  if size(vec,1) == 1
    rind = 1:n;
    cind =  n + 2 - rind' ;
    cind(cind == (n+1)) = 1;
  else
    cind = (1:n)';
    rind =  n + 2 - cind';
    rind(rind == (n+1)) = 1;
  end
  C = vec(toeplitz(cind,rind));
end

function Pi= kdc_ggp(tau, taup, beta, alpha, NTs)
Pi=( exp(-beta/2*(tau+taup)).*exp(-alpha*(taup-tau))  ) ./(exp(NTs*(beta/2+alpha))-1)+ exp(-beta/2*(tau+taup)-alpha*abs(tau-taup)) ;

function Pi= kdc_gp(tau, taup, beta, alpha, NTs)
Pi= exp(-beta/2*(tau+taup)).*(exp(-alpha*abs(tau-taup))  +  exp(alpha*abs(tau-taup)-NTs*(beta/2+alpha) ) )/(1-exp(-beta*NTs))/(1-exp(-NTs*(beta/2+alpha)));

function Pi= kdc(tau, taup , beta, alpha)
Pi = exp(-beta/2*(tau+taup)).*exp(-alpha*abs(tau-taup));



