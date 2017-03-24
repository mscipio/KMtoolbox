function y = fengInput(parm, t)
% Evaluate Feng input function (FDG)
% y = FengInput(parm, t)
% parm =  [modelNumber delay a(1:end) lambda(1:end)];
%      or [modelNumber delay];
%      or [modelNumber];
%      or [];
% in the latter cases, default values will be used
% t = time in minutes
% y = input function
%
% Example usage:
% t=0:60;
% plot(t,FengInput(2,t))
%
% D Feng, SC Huang, X Wang
% Models for Computer simulation Studies of Input Functions for 
% Tracer Kinetic Modeling with Positron emission Tomography
% International Journal of Biomedical Computing 32(2):95-110, 1993 (March).
% (Int J of Biomed Comput)
%
% Fitting example:
% With a(r,1) containing the time points in sec
%      a(r,2) containing the radioactivity or cps
%      (e.g. r=1:1000);
% fit data to model 2 using
% x0=[.5 16000 2188 2081 -7.5 -10 -0.01043];
% options=optimset('Display','iter','MaxFunEvals',20000, 'MaxIter',2000);
% [p,resnorm,residual,exitflag,output]=lsqcurvefit(inline('fengInput([2 x],xdata)','x','xdata'),x0,a(:,1)/60,a(:,2),[],[],options);
% plot(a(:,1),a(:,2),'.',a(:,1),fengInput([2 p],a(:,1)/60)); 
% current version: $Revision: 1.1 $ $Date: 2007/10/08 19:46:34 $

if length(parm) < 1,
   modelNumber = 2;
else
   modelNumber = parm(1);
end

switch modelNumber
    case 2
        if length(parm) < 2,
            delay = 0.735;
        else
            delay = parm(2);
        end
        
        if length(parm) < 3,
            a = [851.1 21.88 20.81];
            lambda = [-4.134 -0.1191 -0.01043];
        else
            a=parm(3:5);
            lambda=parm(6:8);
        end
        
        y = zeros(size(t));
        idx = find(t > delay);
        y(idx) = (a(1)*(t(idx)-delay)-a(2)-a(3)).*exp(lambda(1)*(t(idx)-delay)) + ...
                  a(2)*exp(lambda(2)*(t(idx)-delay)) + ...
                  a(3)*exp(lambda(3)*(t(idx)-delay));
    case 3
        delay = 0;

        if length(parm) < 3,
            a = [851.1 21.88 20.81];
            lambda = [-4.134 -0.1191 -0.01043];
        else
            a=parm(3:5);
            lambda=parm(6:8);
        end
        
        y = zeros(size(t));
        idx = find(t > delay);
        y(idx) = (a(1)*(t(idx)-delay)-a(2)-a(3)).*exp(lambda(1)*(t(idx)-delay)) + ...
            a(2)*exp(lambda(2)*(t(idx)-delay)) + ...
            a(3)*exp(lambda(3)*(t(idx)-delay));
        
        
    case 4
        if length(parm) < 2,
            delay = 0.731;
        else
            delay = parm(2);
        end
        
        if length(parm) < 3,
            a = [892.5 36.8];
            lambda = [-3.8862 -0.0262];
        else
            a = parm(3:4);
            lambda = parm(5:6);
        end
        
        y = zeros(size(t));
        idx = find(t > delay);
        y(idx) = (a(1)*(t(idx)-delay)-a(2)).*exp(lambda(1)*(t(idx)-delay)) + ...
                                        a(2)*exp(lambda(2)*(t(idx)-delay));
        
    case 5
        if length(parm) < 2,
            delay = 0.781;
        else
            delay = parm(2);
        end
        
        if length(parm) < 3,
            a = [90.64];
            lambda = [-77.53 -0.341];
        else
            a = parm(3);
            lambda = parm(4:5);
        end
        
        y=zeros(size(t));
        idx=find(t > delay);
        y(idx) = a(1)*exp(lambda(1)*(t(idx)-delay)) + ...
            a(1)*exp(lambda(2)*(t(idx)-delay));  % note the second term is + here
                                                 % and - in eq. 5 of Feng 1993
                                                 % this must be a typo in Feng 1993
        
    case 6

        delay = parm(2);
        a = parm(3:4);
        lambda = parm(5:6);
        
        y=zeros(size(t));
        idx=find(t > delay);
        y(idx) = a(1)*exp(lambda(1)*(t(idx)-delay)) + ...
                 a(2)*exp(lambda(2)*(t(idx)-delay));
    
    otherwise
        error(sprintf('FengInput: model %d not implemented', modelNumber))
        y = [];
end

return



