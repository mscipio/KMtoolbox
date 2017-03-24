%% ' Logan plot implementation
%' 
%' Implements a Logan plot, or Logan graphical analysis, a 'graphical method
%' of analysis applicable to ligands that bind reversibly  to receptors or
%' enzymes requiring the simultaneous measurement of plasma and tissue
%' radioactivities for multiple times after the injection of a radiolabeled 
%' tracer is presented'
%' 
%' @param input_function Input function TAC
%' @param tissue Tissue TAC
%' @param time_start Initial acquisition time for each frame
%' @param time_end Final acquisition time for each frame
%' @param plot Should the Logan plot be displayed? Defaults to \code{TRUE}
%' @param ___ Additional parameters passed to the {plot} function
%' 
%' @details Please refer to the documentation of \code{\link{processlinear}}
%'   for implementation details_ 
%' 
%' @return Returns a list with three fields: \code{kparms}, the computed kinetic
%'   parameters; \code{stderrors}, the standard errors for each parameter as a 
%'   percentage; \code{fit}, the actual fitted object_
%'   
%' @references J_ Logan, J_ S_ Fowler, N_ D_ Volkow, A_ P_ Wolf, S_ L_
%' Dewey, D_ J_ Schlyer, R_ R_ MacGregor, R_ Hitzemann, B_ Bendriem, and S_ J_
%' Gatley, 'Graphical analysis of reversible radioligand binding from
%' time-activity measurements applied to [N-11C-methyl]-(-)-cocaine PET studies
%' in human subjects_,' Journal of cerebral blood flow and metabolism : official
%' journal of the International Society of Cerebral Blood Flow and Metabolism,
%' vol_ 10, no_ 5, pp_ 740-7, Sep_ 1990_
%'
%' @seealso {patlak_plot}

function res = logan_plot(IF, tac, scant, plotfig, maxerror)

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
dt = scant(:,2)-scant(:,1);
termy = cumsum(tac .* dt) ./ tac;
termx = cumsum(IF .* dt) ./ tac;

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
    title(['Logan plot -->  K_i=',num2str(res(2))])
    xlabel('Int(Cplasma) / Ctissue [min]')
    ylabel('Int(Ctissue) / Ctissue [min]')
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


