function OUT = dce_models(model, params, scanTime, ppCp, fitting, plotfig)
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


    k  = params(3:end);
    vB = params(1);
    
    if max(scanTime(:))>180
        scanTime = scanTime./60; % time has to be in minutes
    end
    
    % use fine time sampling and implement the convolution
    % using discrete mathematics.  
    dt = .01;  % time step, in minutes    
    tm0 = mean(scanTime,2); %[scanTime(1,1); scanTime(:,2)];
    t = (scanTime(1,1):dt:scanTime(end,2))';
    Cp = max(0,interp1(tm0,ppCp,t,'linear',0) ); 

    sol = conv(FDGImpulseResponse(k, t, model), Cp) * dt;
    tsol = (0:dt:2*max(t))';
    % Taking into account blood fraction in tissue
    sol = max(0,interp1(tsol,sol,tm0,'linear',0) ); 
    Ct = (1-vB)*sol + vB* ppCp;
    Ct(Ct<=0) = eps;
    
    if nargin > 4
        if fitting
            OUT = Ct;
        else
            OUT.tsample = scanTime(:,2);
            OUT.TAC = Ct;
        end
    else       
        OUT.tsample = scanTime(:,2);
        OUT.TAC = Ct;
    end

    if nargin > 5
        if plotfig
            figure,plot(OUT.tsample,ppCp,'r',OUT.tsample,OUT.TAC),
        end
    end

return


function h = FDGImpulseResponse(k,t,model)
    switch(model)
        case 'DCE_MRI 1TC'
            h = reversible_1c_model(k, t);
        case 'DCE_MRI 2TCr'
            h = reversible_2c_model(k, t);
        case 'DCE_MRI 2TCi'
            h = irreversible_2c_model(k, t);
    end

function h = reversible_1c_model(k, t)

    if length(k)<2
        warning('You need to provide at least the value for K1 and k2')
    elseif length(k)>2
        warning('The 1TC model will use only K1 and k2 values')
    end
    
    h = (t >= 0) .* (k(1)*exp(-(k(2))*t));
return

function h = irreversible_2c_model(k, t)

    if length(k)<3
        warning('You need to provide at least the value for K1, k2 and k3')
    elseif length(k)>3
        warning('The 1TC model will use only K1, k2 and k3 values')
    end
    
    h = (t >= 0) .* (k(1)*...
   ((k(2)/(k(2)+k(3)))*exp(-(k(2)+k(3))*t) + (k(3)/(k(2)+k(3)))));
return

function h = reversible_2c_model(k, t)

    if length(k)<4
        warning('You need to provide at least the value for K1, k2, k3 and k4')
    elseif length(k)>4
        warning('The 1TC model will use only K1, k2, k3 and k4 values')
    end
    
    beta = [ ...
        (k(2)+k(3)+k(4)-sqrt((k(2)+k(3)+k(4))^2-4*k(2)*k(4)))/2;
        (k(2)+k(3)+k(4)+sqrt((k(2)+k(3)+k(4))^2-4*k(2)*k(4)))/2];

    h = (t >= 0) .* (k(1)/(beta(2)-beta(1))*...
   ((k(3)+k(4)-beta(1))*exp(-beta(1)*t) + (beta(2)-k(3)-k(4))*exp(-beta(2)*t)));
return


