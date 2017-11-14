function OUT = analytic_models(model, params,scanTime,IFparams, ppCp, fitting, plotfig)

% numeric impelmentation of analytic solution of FDG model
%
% [PET, Ct] = analyicFDG(ppCp, k, scanTime)
%
% Outputs
%   OUT.tsample - resampled time vector
%   OUT.IF      - interpolated input function Cp (pmol/ml) in respect to 
%                 OUT.tsample
%   OUT.TAC     - instantaneous tissue concentration (pmol/ml) = Ce + Cm
%                 multiplied by exponentially decaying specific activity (uCi/ml)
%
% Inputs
%   model         - label identifying the chosen kinetic model
%   k             - vector of rate constants of kinetic model [k1,k2,k3,k4]
%   scanTime(:,1) - frame start times
%   scanTime(:,2) - frame end times
%   ppCp          - input function Cp (pmol/ml)
%   sa            - specific activity measured at t=0 (uCi/pmol)
%   vB            - blood fraction in tissue (% value [0,1])
%   plot          - [boolean] auto plot both OUT.IF and OUT.TAC
%
% SVN Version Information__________________________________________________
% $Rev: 1 $:
% $Author: Michele Scipioni $:
% $Date: 2016-11-11 $:
% $Id: compartmentalModels.m $:
% _________________________________________________________________________

    vB = params(1);
    sa = params(2);
    k  = params(3:end);
    
    dk = log(2)/109.8;  % radioactive decay constant for F-18
    if max(scanTime(:))>180
        scanTime = scanTime./60; % time has to be in minutes
    end 
    time = mean(scanTime,2);%[scanTime(1,1); scanTime(:,2)];

    sol = FDGImpulseResponse([k, vB], time, IFparams, model);
    sol = sol * sa .* exp(-dk * time);
    Ct = (1-vB)*sol + vB* ppCp;    
    Ct(Ct<=0) = eps;
    
    if nargin > 4
        if fitting
            OUT = Ct;
        else
            OUT.tsample = time;
            OUT.TAC = Ct;
        end
    else       
        OUT.tsample = time;
        OUT.TAC = Ct;
    end

    if nargin > 5
        if plotfig
            figure,plot(OUT.tsample,ppCp,'r',OUT.tsample,OUT.TAC),
        end
    end
    
    
    
 function h = FDGImpulseResponse(k,t,IFpar,model)   
     switch(model)
         case '2TCr_analytic'
             h = reversible_2c_model(k, t, IFpar);
         case '2TCr_analytic_nodelay'
             h = reversible_2c_model_nodelay(k, t, IFpar);
     end


 function h = reversible_2c_model(k, time, IFparams)
     
% auxiliar parameters
d    = abs(sqrt((k(2)+k(3)+k(4)).^2 - 4*k(2)*k(4)));
p(2) = (k(2) + k(3) + k(4) + d) ./ 2;   %L1
p(4) = (k(2) + k(3) + k(4) - d) ./ 2;   %L2
p(1) = (k(1)*(p(2) - k(3) - k(4)))./ d; %B1
p(3) = (k(1)*(-p(4) + k(3) + k(4)))./ d; %B2
p(5) = k(5); % vB = k5

delay = IFparams(1);
idx = find(time > delay);
t=time(idx)-delay;

A = IFparams(2:4); % [A1 A2 A3]
lambda = -IFparams(5:7); % [lambda1 lambda2 lambda3] has to be > 0

h = zeros(size(time,1),1);

for i = 1:2:4 % 2 compartimenti
    sumTerm = 0;
    Bi = p(i);
    Li = p(i+1);
    Ahat = [-A(2)-A(3)-(A(1)/(Li -lambda(1)))   A(2)  A(3)];
    
    for j = 1:3 % IF modeled as Feng
        Ahat_j = Ahat(j);
        lambda_j = lambda(j);
        delta = 1./(Li-lambda_j);
        
        sumTerm = sumTerm + Ahat_j*delta*(exp(-lambda_j*t)-exp(-Li*t));
    end
    sumTerm = sumTerm .*Bi + ((A(1)*Bi*t) / (Li-lambda(1))).*exp(-lambda(1)*t);
    h(idx) = h(idx) + sumTerm;
end

function h = reversible_2c_model_nodelay(k, time, IFparams)
     
% auxiliar parameters
d    = abs(sqrt((k(2)+k(3)+k(4)).^2 - 4*k(2)*k(4)));
p(2) = (k(2) + k(3) + k(4) + d) ./ 2;   %L1
p(4) = (k(2) + k(3) + k(4) - d) ./ 2;   %L2
p(1) = (k(1)*(p(2) - k(3) - k(4)))./ d; %B1
p(3) = (k(1)*(-p(4) + k(3) + k(4)))./ d; %B2
p(5) = k(5); % vB = k5

A = IFparams(2:4); % [A1 A2 A3]
lambda = -IFparams(5:7); % [lambda1 lambda2 lambda3] has to be > 0

h = zeros(size(time,1),1);

for i = 1:2:4 % 2 compartimenti
    sumTerm = 0;
    Bi = p(i);
    Li = p(i+1);
    Ahat = [-A(2)-A(3)-(A(1)/(Li -lambda(1)))   A(2)  A(3)];
    
    for j = 1:3 % IF modeled as Feng
        Ahat_j = Ahat(j);
        lambda_j = lambda(j);
        delta = 1./(Li-lambda_j);
        
        sumTerm = sumTerm + Ahat_j*delta*(exp(-lambda_j*time)-exp(-Li*time));
    end
    sumTerm = sumTerm .*Bi + ((A(1)*Bi*time) / (Li-lambda(1))).*exp(-lambda(1)*time);
    h = h + sumTerm;
end

