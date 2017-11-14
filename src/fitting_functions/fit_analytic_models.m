function [fit, params, resnorm, residual, output] = fit_analytic_models (tac,scanTime,IFparams,ppCp, model)

warning off;

if max(scanTime(:))>180
    scanTime = scanTime./60; % time has to be in minutes
end

% set reasonable start values for fitting parameters
switch(model)
    case '2TCr_analytic'
        %         vB  sa   K1   k2    k3    k4
        par0  = [0.01  1  0.1  0.1  0.01  0.01];
        fixed = [0     0    0    0     0     0];
    case '2TCr_analytic_nodelay'
        %         vB  sa   K1   k2    k3    k4
        par0  = [0.01  1  0.1  0.1  0.01  0.01];
        fixed = [0     0    0    0     0     0];
end

% create options
options = optimset('Display', 'none');
options.TolFun=1e-6;
options.TolX=1e-6;
options.MaxFunEvals=1e6;
options.MaxIter=1e6;
lb    = [0.  0.  0.   0.   0.   0.]; % all parameters should be at least positive
ub    = [1.  1.  8.   8.   8.   8.]; % and at most 1

% fit = zeros(size(tac))+eps;
% params = zeros(6,size(tac,2))+eps;

% func = @(k,t) analyticFDG(k, t, ppCp, sa, tol);
func = @(par,t) analytic_models(model, par, t, IFparams, ppCp, 1, 0);
[params,resnorm,residual,~,output] = fit_nl(func, par0, scanTime, tac, fixed, lb, ub, options);

Ki = params(3)*params(5)/(params(4)+params(5)); %k1*k3/k2+k3
params = [params Ki]; % vB sa k1 k2 k3 k4 Ki

fit = analytic_models(model, params, scanTime, IFparams, ppCp, 1, 0);

% figure,
% plot(mean(scanTime,2),tac,'*'), hold on
% plot(mean(scanTime,2),fit,'-')
% plot(mean(scanTime,2),ppCp)
% title([{['v_B=',num2str(params(1)),'  ||  sa=',num2str(params(2)),...
%     '  ||  K_i=',num2str(params(7))]};...
%     {['k=',num2str(params(3:6))]}])
% legend('Measured TAC','Fitted TAC','Input function')
