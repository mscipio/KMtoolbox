<<<<<<< HEAD
<<<<<<< HEAD
%' Patlak plot implementation
%' 
%' Implements a Patlak plot linear analysis, another graphical method for 
%' analysis of tracers that can be modeled after an irreversible two-tissue 
%' compartment model.
%' 
%' @param input.function Input function TAC.
%' @param tissue Tissue TAC.
%' @param time.start Initial acquisition time for each frame.
%' @param time.end Final acquisition time for each frame.
%' @param plot Should the Patlak plot be displayed? Defaults to \code{TRUE}.
%' @param maxerror \code{maxerror} parameter to be passed to
%'  \code{processlinear}.
%' @param ... Additional parameters passed to the \code{plot} function.
%' 
%' @details Please refer to the documentation of \code{\link{processlinear}}
%'   for implementation details. 
%'   
%' @return Returns a list with three fields: \code{kparms}, the computed kinetic
%'   parameters; \code{stderrors}, the standard errors for each parameter as a 
%'   percentage; \code{fit}, the actual fitted object.
%'   
%' @references C. Patlak, R. Blasberg, and J. Fenstermacher, "Graphical
%'   evaluation of blood-to-brain transfer constants from multiple-time uptake
%'   data," J Cereb Blood Flow Metab, 1983.
%' 
%' @seealso \code{\link{logan.plot}}.

function res = patlak_plot (IF, tac, scant, plotfig, maxerror)

if max(scant(:))>180
    scant = scant./60; % time has to be in minutes
end

if nargin<4
    maxerror=0.03;
    plotfig = 0;
end
if nargin<5
    maxerror=0.03;
end

% Compute data on both axis
% IF = IF(2:end);
% tac = tac(2:end);
termy = tac ./ IF;
dt = scant(:,2)-scant(:,1);
termx = cumsum(IF .* dt) ./ IF;

% The fitting and plotting is done by this helper function
res = processlinear(termx, termy, maxerror, plotfig);
end

                     
%% ' Helper function for linear models
%' 
%' This function helps process linear models, such as {patlak_plot}
%' or {logan_plot}. It finds the optimal fitting point by calling
%' {findbestfit} and proceeeds to compute the kinetic parameters,
%' that will mostly consist of the intercept and slope of the linear regression,
%' along with their standard errors.
%'
%' If {plotfig} is set to {TRUE}, the regression will be shown on screen.
%' The data points used to compute the regression are displayed with a filled
%' style, while the ones that have not been used are displayed using open
%' circles.
%'
%' @param termx,termy - The x and y values of the graph_
%' @param plotfig     - Should the result be shown?
%' @param maxerror    - {maxerror} parameter to be passed to {findbestfit}
%'
%' @out - Returns a list with four fields: 
%'     {kparms}     - the computed kinetic parameters; 
%'     {stderrors}  - the standard errors for each parameter;
%'     {stderrorsp} - the standard errors for each parameter as a percentage;
%'     {fit}        - the actual fitted object_

function out = processlinear (termx, termy, maxerror, plotfig)

fit = findbestfit(termx, termy, maxerror);
t0 = fit.initial_time_point;
mdl = fit.fit;

res = mdl.Coefficients.Estimate; % [intercepts, slope]
% stdabs = coef(summary(fit.fit))[, 2];
% stderr = data_frame(intercept = stdabs[1], slope = stdabs[2]);
% stderrors = (coef(summary(fit.fit))[, 2] / abs(res)) * 100;
% stderrp = data_frame(intercept = stderrors[1], slope = stderrors[2]);

if plotfig
    figure,
    plot(termx, termy,'LineStyle','none',...
        'Marker','o','MarkerEdgeColor','black'); hold on;
    plot(termx(t0:end), termy(t0:end),'LineStyle','none',...
        'Marker','o','MarkerEdgeColor','black',...
        'MarkerFaceColor','black')
    plot(termx, res(1) + res(2)*termx,'r-')
    xlim([min(termx) max(termx)])
    ylim([min(termy) max(termy)])
    title(['Patlak plot -->  K_i=',num2str(res(2))])
    xlabel('Int(Cplasma) / Cplasma [min]')
    ylabel('Ctissue / Cplasma [unitless]')
end

out.kparms = res;
% out.stderrors = stderr;
% out.stderrorsp = stderrp;
out.fit = mdl;
out.t0 = t0;
out.termx = termx;
out.termy = termy;
end

%% ' Finds optimal fitting point
%'
%' Finds the number of points that have to be used for an optimal fit 
%' between the {x} and {y} variables
%'
%' @param x,y       - Regression data
%' @param minpoints - Minimum data points to be used for the regression
%' @param maxerror  - The maximum error rate allowed
%'
%' @details This function finds the earliest point that allows to have a
%'   regression with less than 10% error between the chosen points (this
%'   parameter is controlled by the {maxerror} variable). If no point
%'   allows to have this error rate, the point that yields the minimum
%'   regression error is used.
%'   
%' @out  -  A list containing 
%'    {initial_time_point} - the time point chosen
%'    {fit}                - the actual fitted object

function out = findbestfit (x, y, maxerror, minpoints)

if nargin<4
    minpoints=3;
end

% It is possible that x and y have NAs from the previous computations.
% The easiest way of dealing with this is to remove them
navalues = isnan(x) | isnan(y) | isinf(x) | isinf(y);
x = x(~navalues);
y = y(~navalues);

fitdata = [x, y];
n = size(fitdata,1);
limit = n - minpoints + 1;
respoint = 1;

for i = 1:limit
    fd = fitdata(i:n,:);
    mdl = fitlm(fd(:,1),fd(:,2),'linear','RobustOpts','on');
%     lm1 = lm(y ~ x, data = fd)
    maxresid = max(abs(mdl.Residuals.Raw / fd(:,2)));
    if (maxresid < maxerror)
        respoint = i;
        break
    end
end


lmres = mdl;
out.initial_time_point = respoint;
out.fit = lmres;

end


=======
=======
>>>>>>> 5fa9534adcf332ff010dadcf5e9217186ab53ed8
%' Patlak plot implementation
%' 
%' Implements a Patlak plot linear analysis, another graphical method for 
%' analysis of tracers that can be modeled after an irreversible two-tissue 
%' compartment model.
%' 
%' @param input.function Input function TAC.
%' @param tissue Tissue TAC.
%' @param time.start Initial acquisition time for each frame.
%' @param time.end Final acquisition time for each frame.
%' @param plot Should the Patlak plot be displayed? Defaults to \code{TRUE}.
%' @param maxerror \code{maxerror} parameter to be passed to
%'  \code{processlinear}.
%' @param ... Additional parameters passed to the \code{plot} function.
%' 
%' @details Please refer to the documentation of \code{\link{processlinear}}
%'   for implementation details. 
%'   
%' @return Returns a list with three fields: \code{kparms}, the computed kinetic
%'   parameters; \code{stderrors}, the standard errors for each parameter as a 
%'   percentage; \code{fit}, the actual fitted object.
%'   
%' @references C. Patlak, R. Blasberg, and J. Fenstermacher, "Graphical
%'   evaluation of blood-to-brain transfer constants from multiple-time uptake
%'   data," J Cereb Blood Flow Metab, 1983.
%' 
%' @seealso \code{\link{logan.plot}}.

function res = patlak_plot (IF, tac, scant, plotfig, maxerror)

if max(scant(:))>180
    scant = scant./60; % time has to be in minutes
end

if nargin<4
    maxerror=0.03;
    plotfig = 0;
end
if nargin<5
    maxerror=0.03;
end

% Compute data on both axis
% IF = IF(2:end);
% tac = tac(2:end);
termy = tac ./ IF;
dt = scant(:,2)-scant(:,1);
termx = cumsum(IF .* dt) ./ IF;

% The fitting and plotting is done by this helper function
res = processlinear(termx, termy, maxerror, plotfig);
end

                     
%% ' Helper function for linear models
%' 
%' This function helps process linear models, such as {patlak_plot}
%' or {logan_plot}. It finds the optimal fitting point by calling
%' {findbestfit} and proceeeds to compute the kinetic parameters,
%' that will mostly consist of the intercept and slope of the linear regression,
%' along with their standard errors.
%'
%' If {plotfig} is set to {TRUE}, the regression will be shown on screen.
%' The data points used to compute the regression are displayed with a filled
%' style, while the ones that have not been used are displayed using open
%' circles.
%'
%' @param termx,termy - The x and y values of the graph_
%' @param plotfig     - Should the result be shown?
%' @param maxerror    - {maxerror} parameter to be passed to {findbestfit}
%'
%' @out - Returns a list with four fields: 
%'     {kparms}     - the computed kinetic parameters; 
%'     {stderrors}  - the standard errors for each parameter;
%'     {stderrorsp} - the standard errors for each parameter as a percentage;
%'     {fit}        - the actual fitted object_

function out = processlinear (termx, termy, maxerror, plotfig)

fit = findbestfit(termx, termy, maxerror);
t0 = fit.initial_time_point;
mdl = fit.fit;

res = mdl.Coefficients.Estimate; % [intercepts, slope]
% stdabs = coef(summary(fit.fit))[, 2];
% stderr = data_frame(intercept = stdabs[1], slope = stdabs[2]);
% stderrors = (coef(summary(fit.fit))[, 2] / abs(res)) * 100;
% stderrp = data_frame(intercept = stderrors[1], slope = stderrors[2]);

if plotfig
    figure,
    plot(termx, termy,'LineStyle','none',...
        'Marker','o','MarkerEdgeColor','black'); hold on;
    plot(termx(t0:end), termy(t0:end),'LineStyle','none',...
        'Marker','o','MarkerEdgeColor','black',...
        'MarkerFaceColor','black')
    plot(termx, res(1) + res(2)*termx,'r-')
    xlim([min(termx) max(termx)])
    ylim([min(termy) max(termy)])
    title(['Patlak plot -->  K_i=',num2str(res(2))])
    xlabel('Int(Cplasma) / Cplasma [min]')
    ylabel('Ctissue / Cplasma [unitless]')
end

out.kparms = res;
% out.stderrors = stderr;
% out.stderrorsp = stderrp;
out.fit = mdl;
out.t0 = t0;
out.termx = termx;
out.termy = termy;
end

%% ' Finds optimal fitting point
%'
%' Finds the number of points that have to be used for an optimal fit 
%' between the {x} and {y} variables
%'
%' @param x,y       - Regression data
%' @param minpoints - Minimum data points to be used for the regression
%' @param maxerror  - The maximum error rate allowed
%'
%' @details This function finds the earliest point that allows to have a
%'   regression with less than 10% error between the chosen points (this
%'   parameter is controlled by the {maxerror} variable). If no point
%'   allows to have this error rate, the point that yields the minimum
%'   regression error is used.
%'   
%' @out  -  A list containing 
%'    {initial_time_point} - the time point chosen
%'    {fit}                - the actual fitted object

function out = findbestfit (x, y, maxerror, minpoints)

if nargin<4
    minpoints=3;
end

% It is possible that x and y have NAs from the previous computations.
% The easiest way of dealing with this is to remove them
navalues = isnan(x) | isnan(y) | isinf(x) | isinf(y);
x = x(~navalues);
y = y(~navalues);

fitdata = [x, y];
n = size(fitdata,1);
limit = n - minpoints + 1;
respoint = 1;

for i = 1:limit
    fd = fitdata(i:n,:);
    mdl = fitlm(fd(:,1),fd(:,2),'linear','RobustOpts','on');
%     lm1 = lm(y ~ x, data = fd)
    maxresid = max(abs(mdl.Residuals.Raw / fd(:,2)));
    if (maxresid < maxerror)
        respoint = i;
        break
    end
end


lmres = mdl;
out.initial_time_point = respoint;
out.fit = lmres;

end


<<<<<<< HEAD
>>>>>>> 5fa9534adcf332ff010dadcf5e9217186ab53ed8
=======
>>>>>>> 5fa9534adcf332ff010dadcf5e9217186ab53ed8
