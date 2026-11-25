function  [th, obj, hp] = rfir_eb_fd_grok(kernel,data, pre_ind, configs)
% RFIR_EB_FD_GROK Estimates impulse response in frequency domain using empirical Bayes
% Inputs:
%   kernel  - String specifying the kernel type (only 'DC' supported currently)
%   data    - iddata object containing input and output signals (single experiment)
%   pre_ind - Numeric vector/scalar of prediction indices (converted to column if row)
%   configs     - Cell array with configuration options (e.g., {'frf'}, {'eir'}, {'cut', n})
% Outputs:
%   th      - Estimated transfer function coefficients
%   obj     - Optimized negative log-likelihood
%   hp      - Optimized hyperparameters
% @ Xiaozhu Fang 2025/11/17

narginchk(4,4); % data, [nb nk], kernel
validateattributes(data, {'iddata'},{});
validateattributes(pre_ind, {'numeric'}, {'vector', 'nonnegative'});
if  isrow(pre_ind)
    pre_ind=pre_ind';
end
kernel = upper(kernel);
validatestring(kernel, {'DC'});
validateattributes(configs, {'cell'}, {});

% main body
[hpini, Ut, Yt,lb,ub] = initialize_parameters(data, pre_ind, kernel, configs);
objFunc = @(x)objective_function(x, Ut, Yt, pre_ind, kernel, lb, ub, configs);
[hp, ~, ~] = fminsearch(objFunc,hpini);
[obj1, ~, th1] = objective_function(hpini, Ut, Yt, pre_ind, kernel, lb, ub, configs);
[obj2, ~, th2] = objective_function(hp, Ut, Yt, pre_ind, kernel, lb, ub, configs);
if obj1>obj2

    th=th2;
    obj=obj2;
else
    th=th1;
    obj=obj1;
end

function [hpini, Ut, Yt, lb, ub] = initialize_parameters(z,pre_ind,kernel,configs)
% Inputs:
%   z       - iddata object
%   pre_ind - Prediction indices
%   kernel  - Kernel type ('DC')
%   configs    - Configuration cell array
% Outputs:
%   hpini   - Initial hyperparameters
%   Ut      - Input signal (frequency domain)
%   Yt      - Output signal (frequency domain)
%   lb, ub  - Hyperparameter bounds
y = pvget(z,'OutputData');
u = pvget(z,'InputData');
Ut = u{1};
Yt = y{1};

% Define bounds for DC kernel
if sum(strcmp(kernel, {'DC'}))
    lb = -inf(5,1);
    ub = inf(5,1);
    lbga = -10 * ones(5,1);
    ubga = 10 * ones(5,1);
end
% Optimize initial hyperparameters using MultiStart
ratio=5;
if sum(strcmp(kernel, {'DC' }))
    objFunc =@(x) objective_function(x, Ut, Yt, pre_ind, kernel, lb, ub, configs);
    ms = MultiStart('FunctionTolerance', 1e-6, 'UseParallel', false);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon', 'x0', rand(length(lbga),1), 'objective', objFunc, 'lb', lbga, 'ub', ubga);
     [hpini,~,~]  = run(ms, problem, ratio*length(lbga));
end

function [obj, condNb, xest] = objective_function(hyper, Ut, Yt, pre_ind, kernel, lb, ub, configs)
% NGLGLKLHD_SIMP Computes negative log-likelihood and impulse response
    % Inputs:
    %   hyper   - Hyperparameters [log(alpha), log(beta), log(c1), log(c2), log(sigma2)]
    %   Ut      - Input signal (frequency domain)
    %   Yt      - Output signal (frequency domain)
    %   pre_ind - Prediction indices
    %   kernel  - Kernel type
    %   lb, ub  - Hyperparameter bounds
    %   configs     - Configuration cell array
    % Outputs:
    %   obj     - Negative log-likelihood
    %   condNb  - Conditional number of covariance matrix (optional)
    %   xest    - Estimated Transfer Function (optional)
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
        alpha=exp(hyper(1));
        beta=exp(hyper(2));
        switch configs{1}
            case "frf"
                Wf = (0:(N-1))*2*pi/N;
                UtUtT=Ut.*(Ut');
                PiY = UtUtT.*( exp(-beta)*exp(-1j*(Wf'-Wf))./(1-exp(-beta-1j*(Wf'-Wf))).* (1-exp(-2*alpha - beta-1j*(Wf'-Wf)))./(1-exp(-alpha- beta/2-1j*(Wf')))./(1-exp(-alpha- beta/2+1j*(Wf))));
                Ti = 1./(1-exp(-beta-1j*(Wf'-Wf))).* (1-exp(-2*alpha - beta-1j*(Wf'-Wf)))./(1-exp(-alpha- beta/2-1j*(Wf')))./(1-exp(-alpha- beta/2+1j*(Wf)));
            case "eir"
                u = ifft(Ut)*sqrt(N);
                Phizero=flipud(hankel(flip([0; u(1:end-1)]), zeros(length(u),1)))';
                Phizero=[Phizero(1:N,:)]';
                vec=1:(N);
                K_DC_td= kdc(vec',vec,beta,alpha);
                PiY = ifft(fft(Phizero*K_DC_td*(Phizero'),[],1),[],2);
                yvec= 1: N;
                Ti_temp = exp(-alpha*abs(yvec'-yvec)) .* exp(-beta*(yvec'+yvec)/2);
                Ti = ifft(fft(Ti_temp,[],1),[],2);
        end

        M = c1*PiY+ c2*Ti+ sigma2*eye(N) ;
end

if ini_flag == 0
    obj = inf;
    return;
end

if  cond_flag ~= 1 || chol_flag ~=1
    obj = 1/eps;
else

    try
        cholFactor = chol(M,'lower');
    catch
        try
            cholFactor = chol((M+M')/2,'lower');
        catch
            obj = 1/eps;
            xest=zeros(length(pre_ind),1);
            condNb = nan;
            return
        end
    end
    a=cholFactor'\(cholFactor\Yt);
    obj = Yt'*a+ 2*sum(log(diag(cholFactor)));
    obj=real(obj);
end

if nargout >1
    if configs{1}== "frf"
        Pipre = exp(-beta)*exp(-1j*(pre_ind-Wf))./(1-exp(-beta-1j*(pre_ind-Wf))).* (1-exp(-2*alpha - beta-1j*(pre_ind-Wf)))./(1-exp(-alpha- beta/2-1j*(pre_ind)))./(1-exp(-alpha- beta/2+1j*(Wf)));
        xest=c1 * Pipre * diag(Ut') * a;
    elseif configs{1} == "eir"
        xest= c1 * ifft(kdc_Gg( beta,alpha, pre_ind, vec) * (Phizero'),[],2) * sqrt(N) * a;
    end
    condNb = cond(M);
end

function Pi= kdc(tau, taup , beta, alpha)
Pi = exp(-beta/2*(tau+taup)).*exp(-alpha*abs(tau-taup));

% function Pipre = kdc_Gg( beta,alpha, Wpre, tvec)
%     % First term
%     num1 = exp(alpha + alpha * (-tvec) + beta * tvec / 2);  % e^{alpha + alpha (-tp) + beta tp / 2}
%     den1 = -exp(alpha) + exp(beta/2 + 1j * Wpre);           % -e^{alpha} + e^{beta/2 + i w}
%     term1 = num1 ./ den1;
% 
%     % Second term
%     num2 = exp(-1j * tvec .* Wpre);                            % e^{-i tp w}
%     den2 = -1 + exp(alpha + beta/2 + 1j * Wpre);            % -1 + e^{alpha + beta/2 + i w}
%     term2 = num2 ./ den2;
% 
%     % Third term
%     num3 = exp(alpha - 1j * tvec .* Wpre);                     % e^{alpha - i tp w}
%     den3 = exp(alpha) - exp(beta/2 + 1j * Wpre);            % e^{alpha} - e^{beta/2 + i w}
%     term3 = num3 ./ den3;
% 
%     % Combine all terms
%     Pipre =  exp(beta * (-tvec)) .* (term1 + term2 + term3);

% function Pipre = kdc_Gg(beta, alpha, Wpre, tvec)
%     % Precompute common exponentials
%     exp_alpha = exp(alpha);
%     exp_beta2_iW = exp(beta/2 + 1j * Wpre);
%     
%     % Compute term1 + term3
%     den1 = -(1 - exp(-alpha) * exp_beta2_iW);  % -(1 - e^{-alpha} * e^{beta/2 + i Wpre})
%     num13 = exp(-alpha * tvec + (beta * tvec) / 2) - exp(-1j * tvec * Wpre);
%     term13 = num13 ./ den1;
%     
%     % Compute term2
%     den2 = -1 + exp_alpha * exp_beta2_iW;  % -1 + e^{alpha + beta/2 + i Wpre}
%     num2 = exp(-1j * tvec * Wpre);
%     term2 = num2 ./ den2;
%     
%     % Combine with outer factor
%     Pipre = exp(-beta * tvec) .* (term13 + term2);


% 
function Pipre = kdc_Gg(beta, alpha, omega, t_i)
    % Compute gamma = beta/2 + i*Wpre
    gamma = beta/2 + 1j * omega;
    
    % Compute delta1 = gamma - alpha and delta2 = alpha + gamma
    delta1 = gamma - alpha;
    delta2 = alpha + gamma;
    
    % Precompute common terms
    exp_neg_i_tw = exp(-1j * t_i .* omega);
    
    % Handle term1 and term3 using delta1
    % To avoid division by near zero, use expm1 for denominator
    tvec_beta_alpha = t_i * (-beta/2 - alpha);
    term1 = exp(tvec_beta_alpha) ./ expm1(delta1);
    term3 = -exp_neg_i_tw ./ expm1(delta1);
    
    % Handle term2 using delta2
    term2 = exp_neg_i_tw ./ expm1(delta2);
    
    % Combine all terms and multiply by exp(-beta*tvec)
    Pipre = term1+ exp(-beta * t_i) .* ( term2 + term3);

