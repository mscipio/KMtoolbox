<<<<<<< HEAD
<<<<<<< HEAD
function [Cp, ppCp, ifParams] = fittinFengInputModel2(IF, time)

% With time(r,1) containing the time points in sec
%      IF(r,1)   containing the radioactivity or cps
%      (e.g. r=1:1000);
% ifParams: [delay, A1, A2, A3, -lambda1, -lambda2, -lambda3]
% Use togheter with IF_models.m

% fit data to model 2 using
x0=[0 100000 5000 20000 -5 -0.05 -0.005];
options=optimset('Display','none','MaxFunEvals',20000, 'MaxIter',2000);
[p,resnorm,residual,exitflag,output]= ...
    lsqcurvefit(inline('IF_models([2 x],xdata)','x','xdata'),...
    x0,time/60,IF,[],[],options);

Cp = IF_models([2 p],time/60);
ppCp = spline(time,IF_models([2 p],time));
ifParams = p';
=======
=======
>>>>>>> 5fa9534adcf332ff010dadcf5e9217186ab53ed8
function [Cp, ppCp, ifParams] = fittinFengInputModel2(IF, time)

% With time(r,1) containing the time points in sec
%      IF(r,1)   containing the radioactivity or cps
%      (e.g. r=1:1000);
% ifParams: [delay, A1, A2, A3, -lambda1, -lambda2, -lambda3]
% Use togheter with IF_models.m

% fit data to model 2 using
x0=[0 100000 5000 20000 -5 -0.05 -0.005];
options=optimset('Display','none','MaxFunEvals',20000, 'MaxIter',2000);
[p,resnorm,residual,exitflag,output]= ...
    lsqcurvefit(inline('IF_models([2 x],xdata)','x','xdata'),...
    x0,time/60,IF,[],[],options);

Cp = IF_models([2 p],time/60);
ppCp = spline(time,IF_models([2 p],time));
ifParams = p';
<<<<<<< HEAD
>>>>>>> 5fa9534adcf332ff010dadcf5e9217186ab53ed8
=======
>>>>>>> 5fa9534adcf332ff010dadcf5e9217186ab53ed8
