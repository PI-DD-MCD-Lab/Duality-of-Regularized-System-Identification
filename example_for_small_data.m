% kernel-based regularization method for both DT impulse response and
% frequency response function
% 
% The key feature: both the system (impusle response or FRF) and the transient are
% modelled and estiamted in a nonparametric way. It is signficantly better
% than ML/PEM when the data is short (i.e., transient is severa)
%
% More details can refer to a preprint paper submmited to TAC
%   Xiaozhu Fang, Yu Xu, Tianshi Chen. Duality of Time and Frequency Domain
%   Kernel-Based Regularized System Identification
%
% @ Xiaozhu Fang 2025/11/17

% More efficienct codes with matrix approxiamtion can refer to Xuyu's codes

%% Initialization
clear; clc; close all;  % Simplified clear commands
echo off
configs{3} = 5; % MultiStart runs per hyperparameter (5 gives 25 starts).


% 
% c = parcluster('local');
% c.NumWorkers =70;
% parpool(70);



%% Main processing
databank = {'s1d1star_sysfrom1to10'};
%	S1D1$^*$: fast systems,  $N = 500$, SNR = 10;
%   S2D1$^*$: slow systems,  $N = 500$, SNR = 10;
%	S1D2$^*$: fast systems,  $N = 375$, SNR = 1;
%	S2D2$^*$: slow systems, $N = 375$, SNR = 1.
num_iters = 10;
n = 125;

for ibank = 1:numel(databank)
    textfile = databank{ibank};
    load([textfile '.mat']);
    if textfile(4)=='1'
        N = 500;
    elseif textfile(4)=='2'
        N = 375;
    end
    % Preallocate results

      gFIT = cell(num_iters,1); % fit of impulse response (first 125)
      yFIT =  cell(num_iters,1); % fit of predicted output
      GINF =  cell(num_iters,1); % normalized Hinifty norm of $G-\hat{G}$
      OBJ = cell(num_iters,1); % cost function of hyperparameter estimation
      HP = cell(num_iters,1); % hyperparameters
      algorithm = {
          {'DC', {"frf","off",configs{3}}, 'fd'}
          {'DC', {"eir","off",configs{3}}, 'fd'}
          {'DC', {"frf","off",configs{3}}, 'td'}
          {'DC', {"eir","off",configs{3}}, 'td'}
          {'DC', {"eir-no-t","off",configs{3}}, 'td'}
          {'DC', {"eir",125,configs{3}}, 'td'}
          {'DC', {"cut",125,configs{3}}, 'td'}
          {'DC', {"cut",125}, 'pem'}
          };
    for ii = 1:num_iters
        fprintf('Processing iteration %d/%d\n', ii, num_iters);
        upre = randn(50*n,1); % input for prediction
        % Data preparation
        [mt, z] = prepare_data(sys{ii}, dat{ii});
        z = z((end-N+1):end);
        [Y, U, w, we] = compute_spectral_components(z, n); 
        % w: training frequency bins
        % we: predicted frequency bins
        
        % Model estimation
        [th, obj, hp, ypre, TFpre] = estimate_models(z, Y, U, we, n, upre, algorithm);
        
        %% Performance metrics
        [gfit, yfit, Ginf] = calculate_metrics(mt, th, ypre, TFpre, n, we, upre );
        %%
        % Store results
        gFIT{ii} = gfit;
        yFIT{ii} = yfit;
        GINF{ii} = Ginf;
        OBJ{ii} = obj;
        HP{ii} = hp;
    end
    
    save(fullfile(['Dataset_TAC2025_', textfile '.mat']),'yFIT','gFIT','GINF','OBJ','HP','algorithm');
end

%% Helper functions
function [mt, z] = prepare_data(sys_entry, data_entry)
mt = sys_entry;
mt.Ts = 1;
z = data_entry(1:end);
end

function [Y, U, w, we] = compute_spectral_components(z, n_points)
N = length(z.u);
w = 2*pi*(0:N-1)'/N;        % w: training frequency bins
we = 2*pi*(0:(50*n_points-1))'/(50*n_points);        % we: predicted frequency bins

U = fft(z.u)/sqrt(N);  % here we normalize the U and Y to avoid the numerical issue, this will not affect the estimate
Y = fft(z.y)/sqrt(N);  
end

function [th_all, obj_all, hp_all, ypre_all, TFpre_all] = estimate_models(z, Y, U, we, n, upre, algorithm)


num_configs = size(algorithm, 1);
[th_all, obj_all, hp_all] = deal(cell(num_configs, 1));
[ypre_all, TFpre_all] = deal(zeros(50*n, num_configs));

        ratio_yu = std(z.y)/std(z.u);
        z.y= z.y/ratio_yu;
        Y= Y/ratio_yu;
% ratio_yu = std(z.y) / std(z.u);
% z.u = z.u / std(z.u);
% z.y = z.y / std(z.y);
% U = U / std(z.u);  
% Y = Y / std(z.y);  

% we normalize the data, and denormalize the estiamte. The normalization is essential to find proper range of the scaling hyperparameters. 
% This operation can be regard as a part of estimation algorithms and put into each subfunctions.   
for c_idx = 1:num_configs
    kernel = algorithm{c_idx,:}{1};
    configs = algorithm{c_idx,:}{2};
    domain = algorithm{c_idx,:}{3};
    % Model estimation
    switch domain
        case 'fd'
            [fh, obj, hp] = rfir_eb_fd(kernel, iddata(Y, U, 1), we, configs);
            th = ifft(fh, 'symmetric'); % impulse response g(0) to g(n-1)
        case 'td'
            [th, obj, hp] = rfir_eb_td(kernel, z, 1:(n-1), configs);
            th = [0; th]; % impulse response  g(0) to g(n-1), where g(0)=0
        case 'pem'
            sysest=  tfest(z, n4sid(z,'best',  n4sidOptions('Focus','simulation'))) ;
            th = lsim(sysest, [1,zeros(1,n-1)]);
            obj=nan;
            hp=nan;
        case 'pw'
            [th, obj, hp] = rfir_eb_td(kernel, [iddata(zeros(configs{2}-1,1),zeros(configs{2}-1,1)); z], 1:(n-1), configs);
            th = [0; th]; % impulse response  g(0) to g(n-1), where g(0)=0
    end
    
    % Store results
    th_all{c_idx} = th(1:n) * ratio_yu;
    obj_all{c_idx} = obj;
    hp_all{c_idx} = hp;
    
    % Simulation
    mpre = idpoly(1, th' * ratio_yu);
    ypre_all(:,c_idx) = lsim(mpre, upre);
    
    % Frequency response
    switch domain
        case 'fd'
            TFpre_all(:,c_idx) = fh*ratio_yu;
        case 'td'
            TFpre_all(:,c_idx) = squeeze(freqresp(mpre, we));
    end
end
end

function [gfit, yfit, Ginf] = calculate_metrics(mt, th, ypre, TFpre, n, we, upre)
% Time domain metrics
th0 = impulse(mt, 0:(n-1));
th= cell2mat(th');
gfit = 100 - 100*sqrt(sum((th0(1:125) - th(1:125,:)).^2) ./ sum((th0(1:125) - mean(th0(1:125))).^2));

% Simulation metrics
ytrue = sim(mt, upre);
ind_test = (49*n+1):50*n;
yfit = 100 - 100*sqrt(sum((ytrue(ind_test) - ypre(ind_test,:)).^2) ./ sum((ytrue(ind_test) - mean(ytrue(ind_test))).^2));

% Frequency domain metrics
TF0 = squeeze(freqresp(mt, we));
%Ginf = max(abs(TFpre - TF0)) ./ max(abs(TF0)) ;
Ginf = max(abs(TFpre - TF0)) ./ sqrt(sum(abs(TF0).^2)/ length(TF0)*2*pi) ;
end
