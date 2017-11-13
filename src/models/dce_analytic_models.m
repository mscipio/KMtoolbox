<<<<<<< HEAD
<<<<<<< HEAD

function OUT = dce_analytic_models(model, params,scanTime,IFparams, ppCp, fitting, plotfig)

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
    k  = params(3:end);
    
    if max(scanTime(:))>180
        scanTime = scanTime./60; % time has to be in minutes
    end 
    time = mean(scanTime,2);%[scanTime(1,1); scanTime(:,2)];

    sol = FDGImpulseResponse([k, vB], time, IFparams, model);
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
         case 'DCE-MRI an_2TCr'
             h = reversible_2c_model(k, t, IFpar);
         case 'DCE-MRI an_2TCi'
             h = irreversible_2c_model(k, t, IFpar);
     end


function h = reversible_2c_model(k, time, IFparams)
     
% auxiliar parameters
d    = abs(sqrt((k(2)+k(3)+k(4)).^2 - 4*k(2)*k(4)));
p(2) = (k(2) + k(3) + k(4) + d) ./ 2;    %B1
p(4) = (k(2) + k(3) + k(4) - d) ./ 2;    %B2
p(1) = (k(1)*(p(2) - k(3) - k(4)))./ d;  %a1
p(3) = (k(1)*(-p(4) + k(3) + k(4)))./ d; %a2

delay = IFparams(1);
idx = find(time > delay);
t=time(idx)-delay;

A = IFparams(2:3); % [A1 A2]
lambda = -IFparams(4:5); % [lambda1 lambda2 ] has to be > 0

h = zeros(size(time,1),1);

for i = 1:2:4 % 2 compartimenti
    sumTerm = 0;
    Bi = p(i);
    Li = p(i+1);
    Ahat = [-A(2)-(A(1)/(Li -lambda(1)))   A(2)];
    
    for j = 1:2 % IF modeled as Feng model 4
        Ahat_j = Ahat(j);
        lambda_j = lambda(j);
        delta = 1./(Li-lambda_j);
        
        sumTerm = sumTerm + Ahat_j*delta*(exp(-lambda_j*t)-exp(-Li*t));
    end
    sumTerm = sumTerm .*Bi + ((A(1)*Bi*t) / (Li-lambda(1))).*exp(-lambda(1)*t);
    h(idx) = h(idx) + sumTerm;
end

 function h = irreversible_2c_model(k, time, IFparams)
     
% auxiliar parameters
p(2) = k(2) + k(3);         %B1
p(4) = 0;                   %B2
p(1) = (k(1)* k(2))./ p(2); %a1
p(3) = (k(1)* k(3))./ p(2); %a2

delay = IFparams(1);
idx = find(time > delay);
t=time(idx)-delay;

A = IFparams(2:3); % [A1 A2]
lambda = -IFparams(4:5); % [lambda1 lambda2 ] has to be > 0

h = zeros(size(time,1),1);

for i = 1:2:4 % 2 compartimenti
    sumTerm = 0;
    Bi = p(i);
    Li = p(i+1);
    Ahat = [-A(2)-(A(1)/(Li -lambda(1)))   A(2)];
    
    for j = 1:2 % IF modeled as Feng model 4
        Ahat_j = Ahat(j);
        lambda_j = lambda(j);
        delta = 1./(Li-lambda_j);
        
        sumTerm = sumTerm + Ahat_j*delta*(exp(-lambda_j*t)-exp(-Li*t));
    end
    sumTerm = sumTerm .*Bi + ((A(1)*Bi*t) / (Li-lambda(1))).*exp(-lambda(1)*t);
    h(idx) = h(idx) + sumTerm;
end























=======
=======
>>>>>>> 5fa9534adcf332ff010dadcf5e9217186ab53ed8

function OUT = dce_analytic_models(model, params,scanTime,IFparams, ppCp, fitting, plotfig)

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
    k  = params(3:end);
    
    if max(scanTime(:))>180
        scanTime = scanTime./60; % time has to be in minutes
    end 
    time = mean(scanTime,2);%[scanTime(1,1); scanTime(:,2)];

    sol = FDGImpulseResponse([k, vB], time, IFparams, model);
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
         case 'DCE-MRI an_2TCr'
             h = reversible_2c_model(k, t, IFpar);
         case 'DCE-MRI an_2TCi'
             h = irreversible_2c_model(k, t, IFpar);
     end


function h = reversible_2c_model(k, time, IFparams)
     
% auxiliar parameters
d    = abs(sqrt((k(2)+k(3)+k(4)).^2 - 4*k(2)*k(4)));
p(2) = (k(2) + k(3) + k(4) + d) ./ 2;    %B1
p(4) = (k(2) + k(3) + k(4) - d) ./ 2;    %B2
p(1) = (k(1)*(p(2) - k(3) - k(4)))./ d;  %a1
p(3) = (k(1)*(-p(4) + k(3) + k(4)))./ d; %a2

delay = IFparams(1);
idx = find(time > delay);
t=time(idx)-delay;

A = IFparams(2:3); % [A1 A2]
lambda = -IFparams(4:5); % [lambda1 lambda2 ] has to be > 0

h = zeros(size(time,1),1);

for i = 1:2:4 % 2 compartimenti
    sumTerm = 0;
    Bi = p(i);
    Li = p(i+1);
    Ahat = [-A(2)-(A(1)/(Li -lambda(1)))   A(2)];
    
    for j = 1:2 % IF modeled as Feng model 4
        Ahat_j = Ahat(j);
        lambda_j = lambda(j);
        delta = 1./(Li-lambda_j);
        
        sumTerm = sumTerm + Ahat_j*delta*(exp(-lambda_j*t)-exp(-Li*t));
    end
    sumTerm = sumTerm .*Bi + ((A(1)*Bi*t) / (Li-lambda(1))).*exp(-lambda(1)*t);
    h(idx) = h(idx) + sumTerm;
end

 function h = irreversible_2c_model(k, time, IFparams)
     
% auxiliar parameters
p(2) = k(2) + k(3);         %B1
p(4) = 0;                   %B2
p(1) = (k(1)* k(2))./ p(2); %a1
p(3) = (k(1)* k(3))./ p(2); %a2

delay = IFparams(1);
idx = find(time > delay);
t=time(idx)-delay;

A = IFparams(2:3); % [A1 A2]
lambda = -IFparams(4:5); % [lambda1 lambda2 ] has to be > 0

h = zeros(size(time,1),1);

for i = 1:2:4 % 2 compartimenti
    sumTerm = 0;
    Bi = p(i);
    Li = p(i+1);
    Ahat = [-A(2)-(A(1)/(Li -lambda(1)))   A(2)];
    
    for j = 1:2 % IF modeled as Feng model 4
        Ahat_j = Ahat(j);
        lambda_j = lambda(j);
        delta = 1./(Li-lambda_j);
        
        sumTerm = sumTerm + Ahat_j*delta*(exp(-lambda_j*t)-exp(-Li*t));
    end
    sumTerm = sumTerm .*Bi + ((A(1)*Bi*t) / (Li-lambda(1))).*exp(-lambda(1)*t);
    h(idx) = h(idx) + sumTerm;
end























<<<<<<< HEAD
>>>>>>> 5fa9534adcf332ff010dadcf5e9217186ab53ed8
=======
>>>>>>> 5fa9534adcf332ff010dadcf5e9217186ab53ed8
