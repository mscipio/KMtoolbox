function y = fit_func(fun, x, xdata)
%Either calls a function pointer or uses a built-in function.
%
% fit_func is called by fit_nl and fit_nl_ex and either simply forward
% everything and call y=fun(x, xdata) or calls a built-in function
% depending on the keyword given in fun.
% Built-in are 1D-3D Gaussian/Lorentzian models with a variable number of
% Gaussian/Lorentzian peaks.
%
% Syntax
%   function y = fit_func(fun, x, xdata)
%
% Input parameters
%   fun     Either a keyword or a function pointer for a function of type
%           y = myfun(x, xdata)
%           Possible keywords are 'xd_gaussian' or 'xd_lorentzian'
%           with x=1,2,3
%   x       Parameter vector. Number of parameters for built-in functions
%           are given below.
%   xdata   The input data for the model function. For built-in functions
%           depending on the dimension this will be a struct with fields x,
%           y and z which are arrays created by ndgrid(), i.d. they contain
%           the x, y or z values of an (x,y,z) position
%
% Output parameters
%   y       Output value of function fun. Can be multi-dimensional.
%
%
% Built-in functions
% ------------------ (Read this if You want to use them!)
%
% Gaussian/Exponential function and sums of them:
%
% The basic form of one(!) 3D Gaussian peak is:
%   p1 * 2^(-(x-p2)^2/(p5/2)^2 -(y-p3)^2/(p6/2)^2 -(z-p4)^2/(p8/2)^2)
% where p1 is the amplitude (p2, p3, p4) defines the center of the peak and
% (p5, p6, p7) are the FWHMs in direction x, y and z. 7 different
% parameters are needed.
%
% For one 2D Gaussian 5 parameters are needed: p1 - amplitude, (p2, p3)
% center position and (p4, p5) FWHMs.
%   [p1 * 2^(-(x-p2)^2/(p4/2)^2 -(y-p3)^2/(p5/2)^2)]
%
% For one 1D Gaussian 3 parameters are enough: p1 - amplitude, p2 - center,
% p3 - FWHM. [p1 * 2^(-(x-p2)^2/(p3/2)^2)]
%
% The general model of the built-in functions is that we have a constant
% background plus the sum of m different gaussian peaks. Therefore for 1D
% depending on the number of Gaussians we will need 1+3m parameters,
% for 2D 1+5m and for 3D 1+7m parameters.
%
% That means that when the dimension of the gaussian is given by the
% keyword ('1/2/3d_gaussian') and the number of parameters is given by
% length(x) than the built-in function knows how many peaks are wanted.
%
% The order of the parameter is first constant background and then the
% parameters for the first gaussian peak, than next gaussian peaks and so
% on...
%
% Example: Two 1D Gaussian peaks
%   keyword is '1d_gaussian', length(x) must be 7
%   function to be computed will be (pi=x(i)=parameter):
%   p1 + p2*2^(-(x-p3)^2/(p4/2)^2)] + p5*2^(-(x-p6)^2/(p7/2)^2)]
%
%
% Lorentzian functions and sums of them:
%
% This is completely equivalent to the system above except the function
% definitions are different.
% The basic definition of a 3D Lorentzian is:
%   p1 / ((x-p2)^2/(p5/2)^2 + (y-p3)^2/(p6/2)^2 + (z-p4)^2/(p7/2)^2 + 1)
%   with amplitude p1, center (p2, p3, p4) and FWHMS (p5, p6, p7)
% And for a 2D Lorentzian:
%   p1 / ((x-p2)^2/(p4/2)^2 + (y-p3)^2/(p5/2)^2 + 1)
% And for a 1D Lorentzian:
%   p1 / ((x-p2)^2/(p3/2)^2 + 1)
%
% Example: Two 1D Lorentzian peaks
%   keyword is '1d_lorentzian', length(x) must be 7
%   function to be computed will be (pi=x(i)=parameter):
%   p1 + p2 / ((x-p3)^2/(p4/2)^2 + 1) + p5 / ((x-p6)^2/(p7/2)^2 + 1)

if nargin < 3
    error('Not enough arguments given!');
end

% check fun for keywords and select function pointers for built-in
% functions
m = [];
if ischar(fun)
    switch fun
        case '1d_gaussian'
            m.fun = @f_gaussian;
            m.dim  = 1;
        case '2d_gaussian'
            m.fun = @f_gaussian;
            m.dim  = 2;
        case '3d_gaussian'
            m.fun = @f_gaussian;
            m.dim  = 3;
        case '1d_lorentzian'
            m.fun = @f_lorentzian;
            m.dim  = 1;
        case '2d_lorentzian'
            m.fun = @f_lorentzian;
            m.dim  = 2;
        case '3d_lorentzian'
            m.fun = @f_lorentzian;
            m.dim  = 3;
        otherwise
            error('Unknown built-in function!');
    end
    % check if the number of parameters is consistent with chosen built-in
    % function
    if mod(length(x) - 1 , 2 * m.dim + 1) ~= 0
        % not all numbers of parameters are allowed
        error('Size of x is not 1+(2*d+1)m with d-dimensionality and m-number of peaks!');
    end
    % check if x_data is a struct with fields x, y, z (depending on
    % dimensionality)
    if ~isstruct(xdata)
        error('Argument xdata must be struct with fields x(,y,z) for built-in functions!');
    end
    if ~isfield(xdata, 'x')
        error('Field x missing in struct xdata!');
    end
    if m.dim > 1 & ~isfield(xdata, 'y')
        error('Field y missing in struct xdata!');
    end
    if m.dim == 3 & ~isfield(xdata, 'z');
        error('Field z missing in struct xdata!');
    end
    % check if size of x,y and z arrays is the same
    if m.dim == 2 & ~isequal(size(xdata.x), size(xdata.y))
        error('Field x and y in xdata must have same size!');
    end
    if m.dim == 3 & ~isequal(size(xdata.x), size(xdata.y), size(xdata.z))
        error('Field x, y and z in xdata must have same size!');
    end
else
    % not a built-in function just forward
    m.fun = fun;
end

% the rest is simple, just call the function
y = m.fun(x, xdata);

% builtin-functions -------------------------------------------------------

    function y = f_gaussian(x, xdata)
        % we can assume that the number of parameters in x and the
        % structure in xdata is correct because this was checked above

        % the number of different gaussians in the sum is n
        n = (length(x) - 1) / (2 * m.dim + 1);

        % set background
        y = zeros(size(xdata.x)) + x(1);

        % for all peaks
        for ki = 1 : n
            % get the right parameters for this peak out of x
            si = 2 + (ki - 1) * (2 * m.dim + 1);
            p = x(si:si+2*m.dim);
            % depending on dimensionality add a peak
            switch m.dim
                case 1
                    y = y + p(1) * power(2., -(xdata.x - p(2)).^2 / (p(3) / 2.)^2);
                case 2
                    y = y + p(1) * power(2., -(xdata.x - p(2)).^2 / (p(4) / 2.)^2 -(xdata.y - p(3)).^2 / (p(5) / 2.)^2);
                case 3
                    y = y + p(1) * power(2., -(xdata.x - p(2)).^2 / (p(5) / 2.)^2 -(xdata.y - p(3)).^2 / (p(6) / 2.)^2 -(xdata.z - p(4)).^2 / (p(7) / 2.)^2);
            end % no error handling necessary m.dim was set by us
        end
    end

    function y = f_lorentzian(x, xdata)
        % we can assume that the number of parameters in x and the
        % structure in xdata is correct because this was checked above

        % the number of different lorentzians in the sum is n
        n = (length(x) - 1) / (2 * m.dim + 1);

        % set background
        y = zeros(size(xdata.x)) + x(1);

        % for all peaks
        for ki = 1 : n
            % get the right parameters for this peak out of x
            si = 2 + (ki - 1) * (2 * m.dim + 1);
            p = x(si:si+2*m.dim);
            % depending on dimensionality add a peak
            switch m.dim
                case 1
                    y = y + p(1) ./ ((xdata.x - p(2)).^2 / (p(3) / 2.)^2 + 1);
                case 2
                    y = y + p(1) ./ ((xdata.x - p(2)).^2 / (p(4) / 2.)^2 + (xdata.y - p(3)).^2 / (p(5) / 2.)^2 + 1);
                case 3
                    y = y + p(1) ./ ((xdata.x - p(2)).^2 / (p(5) / 2.)^2 + (xdata.y - p(3)).^2 / (p(6) / 2.)^2 + (xdata.z - p(4)).^2 / (p(7) / 2.)^2 + 1);
            end % no error handling necessary m.dim was set by us
        end
    end

end