function [x,resnorm,residual,exitflag,output] = fit_nl(model, x0, xdata, ydata, fixed, lb, ub, opt)
% Extends lsqcurvefit to incorporate built-in functions and fixed parameters.
%
% This function mainly calls the function lsqcurvefit from the Optimization
% toolbox. For supporting automatic operations all Warnings or Display
% outputs of lsqcurvefit are suppressed. Additionally parameters can
% be set as fixed.
%
% Syntax
%   function [x, resnorm] = fit_nl(model, x0, xdata, ydata, fixed, lb, ub, opt)
%
% Input parameters
%   model   Defines the function which is used for the fit. This is either
%           a keyword for built-in functions (see fit_func.m) or a function
%           pointer to a function equivalent to fun used in lsqcurvefit,
%           i.d. function F = myfun(x, xdata)
%   x0      The initial guess for the model parameter to be estimated.
%           For the nubmer of meaning of parameter for the built-in
%           functions, see fit_func.m
%   xdata   Input data for the model function/ indepedent variables.
%           For the xdata which is required by the built-in functions see
%           fit_func.m
%   ydata   Output data to be matched by the model function
%
%   [optional parameters]
%   fixed   Array of 0 or 1 determining which parameters should stay fixed
%           in the fitting process.
%           default: set to zeros(size(x0) - no fixed parameter
%           Even if some parameter are set fixed the model function still
%           should take all parameter (including the fixed).
%   lb      Lower bound for x
%           default: set to -Infinity
%   ub      Upper bound for x
%           default: set to +Infinity
%   opt     the options which can be forwarded to lsqcurvefit
%           default: all warnings and printing will be suppressed
%   Simply pass [], if you want any of these options to be set to default
%   values.
%
% Output parameters
%   x       Resulting paramters that minimize the distance between model
%           and ydata
%   resnorm Quadratic norm of (model(x, xdata) - ydata)
%
% Comment
%   This is a least squares fitting. So the error should be distributed
%   according to a Gaussian distribution around the true curves. Linear
%   parameters can also be fitted, but the algorithm will than be slower
%   than standard regression. N-D function can be given in the model and
%   will be linearized internally.

% parameter checking and setting to default values

% need at least 4 parameter
if nargin < 4
    error('Not enough arguments given!');
end

x0    = double(x0); % lsqcurvefit needs double input
ydata = double(ydata);

s = warning('query', 'all'); % save state of warnings

% the last 4 parameters can be set to default values
% if values are set they are check for equal size
if nargin < 8, opt = [];
    if nargin < 7, ub = [];
        if nargin < 6, lb = [];
            if nargin < 5, fixed = [];
            end, end, end, end
if isempty(opt)
    % no options given, standard is to set Display to off and shutdown
    % warnings so that automatic fitting is facilitated
    opt = optimset('Display','off');
    warning('off', 'all');
end
if isempty(ub)
    ub = zeros(size(x0)) + Inf;
end
if isempty(lb)
    lb = zeros(size(x0)) - Inf;
end
if isempty(fixed)
    fixed = zeros(size(x0));
end
if size(fixed) ~= size(x0) | size(lb) ~= size(x0) | size(ub) ~= size(x0)
    error('Parameter x0, fixed, lb and ub have not same size!');
end

% assemble 'global' struct accessible from subfunction F
p = [];
p.fix = fixed;  % fixed parameters
p.x   = x0;     % initial parameters (needed for knowing which are fixed)
p.fun = model;  % the model function (either keyword or function pointer)

% reduce x0 to parameters not fixed (also lb and ub)
x0 = x0(fixed == 0);
lb = lb(fixed == 0);
ub = ub(fixed == 0);

% the call for the internal Matlab function
[x,resnorm,residual,exitflag,output] = lsqcurvefit(@F, x0, xdata, ydata, lb, ub, opt);

% mix fixed parameters in again
p.x(fixed == 0) = x;
x = p.x;

% restore warning settings
warning(s);

% functions ---------------------------------------------------------------

    function y = F(x, xdata)
        % Computes the output of a possible n-D function for given
        % parameters x and xdataset xdata while fixed parameters which are
        % defined in the global variable p are included. The function is
        % defined by p.m in connection with function fit_func().

        % first mix fixed and variable parameters
        h = x;
        x = p.x;
        x(p.fix == 0) = h;

        % calculate the function (external/built-in)
        y = fit_func(p.fun, x, xdata);
    end
end
