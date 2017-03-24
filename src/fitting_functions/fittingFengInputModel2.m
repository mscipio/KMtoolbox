function [Cp, ppCp, ifParams] = fittingFengInputModel2(IF, time, model)

% With time(r,1) containing the time points in sec
%      IF(r,1)   containing the radioactivity or cps
%      (e.g. r=1:1000);
% ifParams: [delay, A1, A2, A3, -lambda1, -lambda2, -lambda3]
% Use togheter with fengInput.m

% fit data to model 2 using

ord = double(floor(single(log(max(IF))./log(10))));
if ord<0; ord=-0.5; end

if model == 2 || model == 3
    x0 = double([0 10*10^ord 0.5*10^ord 2*10^ord -5 -0.05 -0.005]);
    lb = double([0 0 0 0 -Inf -Inf -Inf]);
    ub = double([Inf Inf Inf Inf 0 0 0]);
elseif model == 4 || model == 6
    x0=double([0 10*10^ord 0.5*10^ord -5 -0.05]);
    lb = double([0 0 0 -Inf -Inf]);
    ub = double([Inf Inf Inf 0 0]);
elseif model == 5
    x0=double([0 10*10^ord -5 -0.05]);
    lb = double([0 0 -Inf -Inf]);
    ub = double([Inf Inf 0 0]);
end

options=optimset('Display','none','MaxFunEvals',1e10, 'MaxIter',1e10);
func = @(par,t) fengInput([model par],t);
p = lsqcurvefit(func,x0,time/60,IF,lb,ub,options);

Cp = fengInput([model p],time/60);
ppCp = spline(time/60,fengInput([model p],time/60));
ifParams = p';
