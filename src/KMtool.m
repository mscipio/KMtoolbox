<<<<<<< HEAD
function varargout = KMtool(varargin)
% KMTOOL MATLAB code for KMtool.fig
%      KMTOOL, by itself, creates a new KMTOOL or raises the existing
%      singleton*.
%
%      H = KMTOOL returns the handle to a new KMTOOL or the handle to
%      the existing singleton*.
%
%      KMTOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KMTOOL.M with the given input arguments.
%
%      KMTOOL('Property','Value',...) creates a new KMTOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before KMtool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to KMtool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help KMtool

% Last Modified by GUIDE v2.5 08-Mar-2017 10:06:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @KMtool_OpeningFcn, ...
    'gui_OutputFcn',  @KMtool_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);

if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT


% --- Executes just before KMtool is made visible.
function KMtool_OpeningFcn(hObject, eventdata, handles, varargin)

if isfield(handles,'FIT'); handles = rmfield(handles,'FIT'); end
if isfield(handles,'S'); handles = rmfield(handles,'S'); end
axes(handles.TAC_preview) ; cla
axes(handles.FittingPlot) ; cla
axes(handles.ResidualPlot); cla
set(handles.slice_num,'String',' ');
set(handles.pixel_num,'String',' ');
set(handles.sd_mean  ,'String',' ');
set(handles.ki_value ,'String',' ');
set(handles.k1_value ,'String',' ');
set(handles.k2_value ,'String',' ');
set(handles.k3_value ,'String',' ');
set(handles.k4_value ,'String',' ');
set(handles.vb_value ,'String',' ');
set(handles.sa_value ,'String',' ');
set(handles.Ki_tag   ,'String','Ki')
set(handles.FittingErrorBox,'String',' ')
set(handles.text23   ,'String',' ')

if length(varargin) == 1
    
    handles.immagine = varargin{1};
    handles.clims = [min(handles.immagine(:)) max(handles.immagine(:))];
    
    handles.slices = size(handles.immagine,3);
    handles.frames = size(handles.immagine,4);
    
    set(handles.slice_slider,'max',handles.slices);
    set(handles.slice_slider,'SliderStep',[1/(handles.slices-1) 1/(handles.slices-1)]);
    set(handles.time_slider,'max',handles.frames);
    set(handles.time_slider,'SliderStep',[1/(handles.frames-1) 1/(handles.frames-1)]);
    
    handles.idx_slice=round(handles.slices/2);
    handles.idx_time=round(handles.frames/2);
    set(handles.title,'String',...
        ['Time frame: ',num2str(handles.idx_time),'  |  Slice: ',num2str(handles.idx_slice)]);
    set(handles.time_slider,'Value',handles.idx_time);
    set(handles.slice_slider,'Value',handles.idx_slice);
    set(handles.clim_slider,'min',handles.clims(1)+eps);
    set(handles.clim_slider,'max',handles.clims(2));
    %     set(handles.clim_slider,'SliderStep',[0.1 0.1]);
    set(handles.clim_slider,'Value',handles.clims(2));
    plotMainDataset(handles);
    
elseif length(varargin) > 1
    h = warndlg('!!! ERROR !!! Input data format error');
else
    
end

% Choose default command line output for KMtool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



% UIWAIT makes KMtool wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = KMtool_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;


% ------------------------------ SLIDERS ----------------------------------
function slice_slider_Callback(hObject, eventdata, handles)
handles.idx_slice=round(get(hObject,'Value'));
set(handles.title,'String',['Time frame: ',num2str(handles.idx_time),'  |  Slice: ',num2str(handles.idx_slice)]);
plotMainDataset(handles);
guidata(hObject, handles);
function slice_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject,'min',1);
set(hObject,'max',10);
set(hObject,'SliderStep',[1/(10-1) 1/(10-1)]);
set(hObject,'Value',1);
handles.idx_slice=round(get(hObject,'Value'));
guidata(hObject, handles);


function time_slider_Callback(hObject, eventdata, handles)
handles.idx_time=round(get(hObject,'Value'));
set(handles.title,'String',['Time frame: ',num2str(handles.idx_time),'  |  Slice: ',num2str(handles.idx_slice)]);
plotMainDataset(handles);
guidata(hObject, handles);
function time_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject,'min',1);
set(hObject,'max',10);
set(hObject,'SliderStep',[1/(10-1) 1/(10-1)]);
set(hObject,'Value',1);
handles.idx_time=round(get(hObject,'Value'));
guidata(hObject, handles);

function clim_slider_Callback(hObject, eventdata, handles)
handles.clims(2) = get(hObject,'Value');
plotMainDataset(handles);
guidata(hObject, handles);
function clim_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject,'min',1);
set(hObject,'max',10);
set(hObject,'Value',1);
guidata(hObject, handles);


% ------------------------------- LOAD FILE -------------------------------
function filetype_popup_Callback(hObject, eventdata, handles)
type = get(hObject,'Value');
options = get(hObject,'String');
ft = options(type);
handles.filetype = ft{1};
guidata(hObject, handles);
function filetype_popup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'Value',1);
options = get(hObject,'String');
ft = options(1);
handles.filetype = ft{1};
guidata(hObject, handles);

function load_file_Callback(hObject, eventdata, handles)
[filename, pathname] = uigetfile(['*',handles.filetype], ['Load 4D ',handles.filetype,' image']);

if isfield(handles,'FIT'); handles = rmfield(handles,'FIT'); end
if isfield(handles,'S'); handles = rmfield(handles,'S'); end
axes(handles.TAC_preview) ; cla
axes(handles.FittingPlot) ; cla
axes(handles.ResidualPlot); cla
set(handles.slice_num,'String',' ');
set(handles.pixel_num,'String',' ');
set(handles.sd_mean  ,'String',' ');
set(handles.ki_value ,'String',' ');
set(handles.k1_value ,'String',' ');
set(handles.k2_value ,'String',' ');
set(handles.k3_value ,'String',' ');
set(handles.k4_value ,'String',' ');
set(handles.vb_value ,'String',' ');
set(handles.sa_value ,'String',' ');
set(handles.Ki_tag   ,'String','Ki')
set(handles.FittingErrorBox,'String',' ')
set(handles.text23   ,'String',' ')

if strcmp(handles.filetype,'.mat') 
    tmp = load([pathname,filename]);
    whichVariables = fieldnames(tmp); 
    if numel(whichVariables) == 1
        handles.immagine = tmp.(whichVariables{1});
        handles.clims = [min(handles.immagine(:)) max(handles.immagine(:))];
        handles.slices = size(handles.immagine,3);
        handles.frames = size(handles.immagine,4);

        set(handles.slice_slider,'max',handles.slices);
        set(handles.slice_slider,'SliderStep',[1/(handles.slices-1) 1/(handles.slices-1)]);
        set(handles.time_slider,'max',handles.frames);
        set(handles.time_slider,'SliderStep',[1/(handles.frames-1) 1/(handles.frames-1)]);


        handles.idx_slice=round(handles.slices/2);
        handles.idx_time=round(handles.frames/2);
        set(handles.title,'String',...
            ['Time frame: ',num2str(handles.idx_time),'  |  Slice: ',num2str(handles.idx_slice)]);
        set(handles.time_slider,'Value',handles.idx_time);
        set(handles.slice_slider,'Value',handles.idx_slice);

        set(handles.clim_slider,'min',handles.clims(1)+eps);
        set(handles.clim_slider,'max',handles.clims(2));
        set(handles.clim_slider,'Value',handles.clims(2));

        plotMainDataset(handles);

        h = warndlg('Image data loaded!');
    else
        h = warndlg('!!! ERROR !!! Input file format error');
    end
    
elseif strcmp(handles.filetype,'.img') || strcmp(handles.filetype,'.dcm')
    [PETData, PETInfo] = read4Ddicom(pathname);
    handles.immagine = PETData;
    handles.S.time_vec= [PETInfo.time(1,1); PETInfo.time(:,2)]; % mean(Info.time,2);
    handles.S.scant = PETInfo.time;
    handles.clims = [min(handles.immagine(:)) max(handles.immagine(:))];
    handles.slices = size(handles.immagine,3);
    handles.frames = size(handles.immagine,4);

    set(handles.slice_slider,'max',handles.slices);
    set(handles.slice_slider,'SliderStep',[1/(handles.slices-1) 1/(handles.slices-1)]);
    set(handles.time_slider,'max',handles.frames);
    set(handles.time_slider,'SliderStep',[1/(handles.frames-1) 1/(handles.frames-1)]);


    handles.idx_slice=round(handles.slices/2);
    handles.idx_time=round(handles.frames/2);
    set(handles.title,'String',...
        ['Time frame: ',num2str(handles.idx_time),'  |  Slice: ',num2str(handles.idx_slice)]);
    set(handles.time_slider,'Value',handles.idx_time);
    set(handles.slice_slider,'Value',handles.idx_slice);

    set(handles.clim_slider,'min',handles.clims(1)+eps);
    set(handles.clim_slider,'max',handles.clims(2));
    set(handles.clim_slider,'Value',handles.clims(2));

    plotMainDataset(handles);

    h = warndlg('Image data loaded!');
end

guidata(hObject, handles);


% --- Executes on button press in transpose_btn.
function transpose_btn_Callback(hObject, eventdata, handles)
handles.immagine = permute(handles.immagine,[2,1,3,4]);
plotMainDataset(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function title_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','Time frame: --  |  Slice: --');
guidata(hObject, handles);

% ---------------------------- COLORMAP -----------------------------------
function colormap_Callback(hObject, eventdata, handles)
values = cellstr(get(hObject,'String'));
handles.cmap=values{get(hObject,'Value')};
plotMainDataset(handles);
guidata(hObject, handles);

function colormap_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
values = cellstr(get(hObject,'String'));
set(hObject,'Value',7);
handles.cmap=values{get(hObject,'Value')};
guidata(hObject, handles);


% ------------------------- LOAD TIME VECTOR ------------------------------
function LoadTime_btn_Callback(hObject, eventdata, handles)

[filename, pathname] = uigetfile('*.mat', 'Load PETInfo file');

tmp = load([pathname,filename]);
whichVariables = fieldnames(tmp);
if numel(whichVariables) == 1
    Info = tmp.(whichVariables{1});
    handles.S.scant = Info.time;
    handles.S.time_vec= mean(Info.time,2); % [Info.time(1,1); Info.time(:,2)];
    
    h = warndlg('Time vector loaded!') ;
    
    if isfield(handles.S,'TAC')
        axes(handles.TAC_preview)
        cla
        plot(handles.S.time_vec,handles.S.TAC,'*-');
        xlabel('Time(sec)');ylabel('Counts')
        xlim([min(handles.S.time_vec) max(handles.S.time_vec)] )
    end
    
else
    h = warndlg('!!! ERROR !!! Input file format error');
end
guidata(hObject, handles);

% ------------------------- ROI SELECTION ---------------------------------
function roi_btn_Callback(hObject, eventdata, handles)

%plotMainDataset(handles);
handles = externalFigureRoiSelection(hObject,handles);

h = imfreehand; %position = wait(h);
ROI = h.createMask();

k = 0;
for i = 1:size(ROI,1)
    for j = 1:size(ROI,2)
        if ROI(i,j)==0
            continue
        else
            k=k+1;
            TAC(k,:) = squeeze(handles.immagine(i,j,handles.idx_slice,:));
        end
    end
end

set(handles.slice_num,'String',num2str(handles.idx_slice));
set(handles.pixel_num,'String',num2str(size(TAC,1)));
%set(handles.sd_mean,'String',num2str(mean(std(TAC)./mean(TAC),'omitnan')));

handles.S.TAC   = mean(TAC)'; %[0; mean(TAC)'];
handles.S.slice = handles.idx_slice;
handles.S.ROI   = ROI;

axes(handles.TAC_preview)
cla
if ~isfield(handles.S,'time_vec')
    plot(handles.S.TAC,'*-');
    xlabel('Frame number');ylabel('Counts')
    xlim([0 length(handles.S.TAC)])
else
    plot(handles.S.time_vec,handles.S.TAC,'*-');
    xlabel('Time(sec)');ylabel('Counts')
    xlim([min(handles.S.time_vec) max(handles.S.time_vec)] )
end

close(handles.H)
guidata(hObject, handles);

% ---------------------- SEND TACs TO MODELING TOOL -----------------------

function setIF_btn_Callback(hObject, eventdata, handles)

if ~isempty(handles.S)
    if isfield(handles.S,'time_vec')
        handles.FIT.IF     = double(handles.S.TAC);
        handles.FIT.time   = handles.S.time_vec;
        handles.FIT.scant  = handles.S.scant;
        
        if isfield(handles.FIT,'IF_fit')
            handles.FIT = rmfield(handles.FIT,'IF_fit');
            handles.FIT = rmfield(handles.FIT,'IF_params');
            set(handles.IFfit_chk,'Value',0.0);
        end
        set(handles.IFmeas_chk,'Value',1.0);        
        plotFIT(hObject, handles);
    else
        h = warndlg('You are currently using a uniform time vector!') ;
        handles.FIT.IF     = double(handles.S.TAC);
        handles.FIT.scant  = [linspace(0,length(handles.S.TAC)-1,length(handles.S.TAC)); ...
                                linspace(1,length(handles.S.TAC),length(handles.S.TAC))]';
        handles.FIT.time   = mean(handles.FIT.scant,2);
        if isfield(handles.FIT,'IF_fit')
            handles.FIT = rmfield(handles.FIT,'IF_fit');
            handles.FIT = rmfield(handles.FIT,'IF_params');
            set(handles.IFfit_chk,'Value',0.0);
        end
        set(handles.IFmeas_chk,'Value',1.0);        
        plotFIT(hObject, handles);
    end
else
    h = warndlg('You need to draw a ROI first!') ;
end
guidata(hObject, handles);

% --- Executes on button press in setTAC_btn.
function setTAC_btn_Callback(hObject, eventdata, handles)

if ~isempty(handles.S)
    if isfield(handles.S,'time_vec')
        handles.FIT.tissue = double(handles.S.TAC);
        handles.FIT.time   = handles.S.time_vec;
        handles.FIT.scant  = handles.S.scant;
        disp(max(handles.S.scant(:)))
        
        if isfield(handles.FIT,'tissue_fit')
            handles.FIT = rmfield(handles.FIT,'tissue_fit');
            handles.FIT = rmfield(handles.FIT,'tissue_params');
            handles.FIT = rmfield(handles.FIT,'residual');
            set(handles.TACfit_chk,'Value',0.0);
            axes(handles.ResidualPlot)
            cla
        end
        set(handles.TACmeas_chk,'Value',1.0);
        plotFIT(hObject, handles);
    else
        h = warndlg('You are currently using a uniform time vector!') ;
        handles.FIT.tissue = double(handles.S.TAC);
        handles.FIT.scant  = [linspace(0,length(handles.S.TAC)-1,length(handles.S.TAC)); ...
                                linspace(1,length(handles.S.TAC),length(handles.S.TAC))]';
        handles.FIT.time   = mean(handles.FIT.scant,2);
        if isfield(handles.FIT,'tissue_fit')
            handles.FIT = rmfield(handles.FIT,'tissue_fit');
            handles.FIT = rmfield(handles.FIT,'tissue_params');
            handles.FIT = rmfield(handles.FIT,'residual');
            set(handles.TACfit_chk,'Value',0.0);
            axes(handles.ResidualPlot)
            cla
        end
        set(handles.TACmeas_chk,'Value',1.0);
        plotFIT(hObject, handles);
    end
else
    h = warndlg('You need to draw a ROI first!') ;
end
guidata(hObject, handles);




% ------------------------ EXPORT DATA FROM GUI ---------------------------
function saveROI_btn_Callback(hObject, eventdata, handles)
if isfield(handles,'S')
    ROI = handles.S.ROI;
    slice_num = handles.S.slice;
    TAC = handles.S.TAC;
    uisave({'ROI','slice_num','TAC'},'ROI.mat')
end
guidata(hObject, handles);

function saveFIT_btn_Callback(hObject, eventdata, handles)
if isfield(handles,'FIT')
    FIT = handles.FIT;
    uisave({'FIT'},'FIT.mat')
end
guidata(hObject, handles);

% ------------------------------ CHECKBUTTONS -----------------------------
function IFmeas_chk_Callback(hObject, eventdata, handles)
plotFIT(hObject, handles);guidata(hObject, handles);

function IFfit_chk_Callback(hObject, eventdata, handles)
plotFIT(hObject, handles);guidata(hObject, handles);

function TACmeas_chk_Callback(hObject, eventdata, handles)
plotFIT(hObject, handles);guidata(hObject, handles);

function TACfit_chk_Callback(hObject, eventdata, handles)
plotFIT(hObject, handles);guidata(hObject, handles);

% ---------------------------TAC MODEL SELECTION --------------------------
function TACmodels_switch_Callback(hObject, eventdata, handles)

values = cellstr(get(hObject,'String'));
handles.TACmodel=values{get(hObject,'Value')};

if isfield(handles.FIT,'IF_fit')
    IF = handles.FIT.IF_fit;
else
    IF = handles.FIT.IF;
end

switch handles.TACmodel
    
    case {'DCE_MRI 1TC','DCE_MRI 2TCr','DCE_MRI 2TCi'}
        [fit, params, resnorm, residual, output] = fit_dce_models (handles.FIT.tissue,handles.FIT.scant,IF,handles.TACmodel);
        handles.FIT.tissue_fit = fit;
        handles.FIT.tissue_params = params;
        handles.FIT.residual = residual;

        msg1 = struct2cell(output);
        msg2 = fieldnames(output);
        set(handles.FittingErrorBox,'String',msg1(1:6))
        set(handles.text23,'String',msg2(1:6))
        set(handles.TACfit_chk,'Value',1.0);
                    
        plotFIT(hObject, handles);
        plotResidual(hObject,handles)
        
        set(handles.ki_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(7),3)));
        set(handles.k1_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(3),3)));
        set(handles.k2_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(4),3)));
        set(handles.k3_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(5),3)));
        set(handles.k4_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(6),3)));
        set(handles.vb_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(1),3)));
        set(handles.sa_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(2),3)));
        set(handles.Ki_tag,'String','Ki')
    
    case 'Logan'
        
        if isfield(handles.FIT,'residual')
            handles.FIT = rmfield(handles.FIT,'tissue_fit');
            handles.FIT = rmfield(handles.FIT,'tissue_params');
            handles.FIT = rmfield(handles.FIT,'residual');
            axes(handles.ResidualPlot)
            cla
        end
        
        res = logan_plot(IF, handles.FIT.tissue, handles.FIT.scant);
        
        axes(handles.FittingPlot)
        cla
        plot(res.termx, res.termy,'LineStyle','none',...
            'Marker','o','MarkerEdgeColor','black'); hold on;
        plot(res.termx(res.t0:end), res.termy(res.t0:end),'LineStyle','none',...
            'Marker','o','MarkerEdgeColor','black',...
            'MarkerFaceColor','black')
        plot(res.termx, res.kparms(1) + res.kparms(2)*res.termx,'r-')
        xlim([min(res.termx) max(res.termx)])
        ylim([min(res.termy) max(res.termy)])
        xlabel('Int(Cplasma) / Ctissue [min]')
        ylabel('Int(Ctissue) / Ctissue [min]')
        
        handles.FIT.tissue_fit = res.fit;
        handles.FIT.tissue_params = res.kparms;
        
        set(handles.ki_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(2),3)));
        set(handles.k1_value,'String',' ');
        set(handles.k2_value,'String',' ');
        set(handles.k3_value,'String',' ');
        set(handles.k4_value,'String',' ');
        set(handles.vb_value,'String',' ');
        set(handles.sa_value,'String',' ');
        set(handles.Ki_tag,'String','DV')
        
        
    case 'Patlak'
                
        if isfield(handles.FIT,'residual')
            handles.FIT = rmfield(handles.FIT,'tissue_fit');
            handles.FIT = rmfield(handles.FIT,'tissue_params');
            handles.FIT = rmfield(handles.FIT,'residual');
            axes(handles.ResidualPlot)
            cla
        end
        
        res = patlak_plot(IF, handles.FIT.tissue, handles.FIT.scant);
        
        axes(handles.FittingPlot)
        cla
        plot(res.termx, res.termy,'LineStyle','none',...
            'Marker','o','MarkerEdgeColor','black'); hold on;
        plot(res.termx(res.t0:end), res.termy(res.t0:end),'LineStyle','none',...
            'Marker','o','MarkerEdgeColor','black',...
            'MarkerFaceColor','black')
        plot(res.termx, res.kparms(1) + res.kparms(2)*res.termx,'r-')
        xlim([min(res.termx) max(res.termx)])
        ylim([min(res.termy) max(res.termy)])
        xlabel('Int(Cplasma) / Ctissue [min]')
        ylabel('Int(Ctissue) / Ctissue [min]')
        
        handles.FIT.tissue_fit = res.fit;
        handles.FIT.tissue_params = res.kparms;
        
        set(handles.ki_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(2),3)));
        set(handles.k1_value,'String',' ');
        set(handles.k2_value,'String',' ');
        set(handles.k3_value,'String',' ');
        set(handles.k4_value,'String',' ');
        set(handles.vb_value,'String',' ');
        set(handles.sa_value,'String',' ');
        set(handles.Ki_tag,'String','Ki')
        
    case 'Tissue modeling'
        return;
        
    case {'2TCr_analytic','2TCr_analytic_nodelay'}
        [fit, params, resnorm, residual, output] = fit_analytic_models (handles.FIT.tissue,handles.FIT.scant,handles.FIT.IF_params,handles.FIT.IF_fit, handles.TACmodel);
        handles.FIT.tissue_fit = fit;
        handles.FIT.tissue_params = params;
        handles.FIT.residual = residual;

        msg1 = struct2cell(output);
        msg2 = fieldnames(output);
        set(handles.FittingErrorBox,'String',msg1(1:6))
        set(handles.text23,'String',msg2(1:6))
        set(handles.TACfit_chk,'Value',1.0);
                    
        plotFIT(hObject, handles);
        plotResidual(hObject,handles)
        
        set(handles.ki_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(7),3)));
        set(handles.k1_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(3),3)));
        set(handles.k2_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(4),3)));
        set(handles.k3_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(5),3)));
        set(handles.k4_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(6),3)));
        set(handles.vb_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(1),3)));
        set(handles.sa_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(2),3)));
        set(handles.Ki_tag,'String','Ki')
        
    case {'DCE-MRI an_2TCr','DCE-MRI an_2TCi'}
        disp(max(handles.FIT.scant(:)))
        [fit, params, resnorm, residual, output] = fit_dce_analytic_models (handles.FIT.tissue,handles.FIT.scant,handles.FIT.IF_params,handles.FIT.IF_fit, handles.TACmodel);
        handles.FIT.tissue_fit = fit;
        handles.FIT.tissue_params = params;
        handles.FIT.residual = residual;

        msg1 = struct2cell(output);
        msg2 = fieldnames(output);
        set(handles.FittingErrorBox,'String',msg1(1:6))
        set(handles.text23,'String',msg2(1:6))
        set(handles.TACfit_chk,'Value',1.0);
                    
        plotFIT(hObject, handles);
        plotResidual(hObject,handles)
        
        set(handles.ki_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(7),3)));
        set(handles.k1_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(3),3)));
        set(handles.k2_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(4),3)));
        set(handles.k3_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(5),3)));
        set(handles.k4_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(6),3)));
        set(handles.vb_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(1),3)));
        set(handles.sa_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(2),3)));
        set(handles.Ki_tag,'String','Ki')
        
    otherwise
        [fit, params, resnorm, residual, output] = fit_compartmental_models (handles.FIT.tissue,handles.FIT.scant,IF,handles.TACmodel);
        handles.FIT.tissue_fit = fit;
        handles.FIT.tissue_params = params;
        handles.FIT.residual = residual;

        msg1 = struct2cell(output);
        msg2 = fieldnames(output);
        set(handles.FittingErrorBox,'String',msg1(1:6))
        set(handles.text23,'String',msg2(1:6))
        set(handles.TACfit_chk,'Value',1.0);
                    
        plotFIT(hObject, handles);
        plotResidual(hObject,handles)

        set(handles.ki_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(7),3)));
        set(handles.k1_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(3),3)));
        set(handles.k2_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(4),3)));
        set(handles.k3_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(5),3)));
        set(handles.k4_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(6),3)));
        set(handles.vb_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(1),3)));
        set(handles.sa_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(2),3)));
        set(handles.Ki_tag,'String','Ki')
end

guidata(hObject, handles);


function TACmodels_switch_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'Value',1);


% --------------------------- IF MODEL SELECTION --------------------------
function IFmodels_switch_Callback(hObject, eventdata, handles)
handles.IFmodel=get(hObject,'Value');

switch handles.IFmodel
    case 1
        return;
    otherwise
        [Cp, ppCp, ifParams] = fittingFengInputModel2(handles.FIT.IF, handles.FIT.time, handles.IFmodel);
        handles.FIT.IF_fit = Cp;
        handles.FIT.IF_params = ifParams;
end
isTissue = isfield(handles.FIT,'tissue');
isIF = isfield(handles.FIT,'IF');
isTissue_fit = isfield(handles.FIT,'tissue_fit');
isIF_fit = isfield(handles.FIT,'IF_fit');
checks = [isTissue, isIF, isTissue_fit, isIF_fit];
if isTissue;     set(handles.TACmeas_chk,'Value',1);end
if isIF;         set(handles.IFmeas_chk,'Value',1);end
if isTissue_fit; set(handles.TACfit_chk,'Value',1);end
if isIF_fit;     set(handles.IFfit_chk,'Value',1);end
plotFIT(hObject, handles);

guidata(hObject, handles);

function IFmodels_switch_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'Value',1);


% ---------------------------------------------------------------------------- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function k1_value_Callback(hObject, eventdata, handles)
function k1_value_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function vb_value_Callback(hObject, eventdata, handles)
function vb_value_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function k4_value_Callback(hObject, eventdata, handles)
function k4_value_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function k3_value_Callback(hObject, eventdata, handles)
function k3_value_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function k2_value_Callback(hObject, eventdata, handles)
function k2_value_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sa_value_Callback(hObject, eventdata, handles)
function sa_value_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ki_value_Callback(hObject, eventdata, handles)
function ki_value_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function slice_num_Callback(hObject, eventdata, handles)
function slice_num_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pixel_num_Callback(hObject, eventdata, handles)
function pixel_num_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sd_mean_Callback(hObject, eventdata, handles)
function sd_mean_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% -------------------------- UTILITY FUNCTION -----------------------------

function ret = optimizeParamVisualiztion(n,dec)
f = 10.^dec;
ret = round(f*n)/f;

function plotMainDataset(handles)
axes(handles.plot_axes)
imagesc(squeeze(handles.immagine(:,:,handles.idx_slice,handles.idx_time)),handles.clims),colormap(handles.cmap);

function handles = externalFigureRoiSelection(hObject, handles)
handles.H = figure(1);
imagesc(squeeze(handles.immagine(:,:,handles.idx_slice,handles.idx_time)),handles.clims),colormap(handles.cmap),axis equal tight;
guidata(hObject, handles);


function plotFIT (hObject, handles)

isTissue = isfield(handles.FIT,'tissue') && get(handles.TACmeas_chk,'Value');
isIF = isfield(handles.FIT,'IF') && get(handles.IFmeas_chk,'Value');
isTissue_fit = isfield(handles.FIT,'tissue_fit') && get(handles.TACfit_chk,'Value');
isIF_fit = isfield(handles.FIT,'IF_fit') && get(handles.IFfit_chk,'Value');


if isTissue;     set(handles.TACmeas_chk,'Value',1);end
if isIF;         set(handles.IFmeas_chk,'Value',1);end
if isTissue_fit; set(handles.TACfit_chk,'Value',1);end
if isIF_fit;     set(handles.IFfit_chk,'Value',1);end

if isTissue; TAC = handles.FIT.tissue;
else TAC = zeros(length(handles.FIT.time))*NaN;end;

if isTissue_fit; TAC_fit = handles.FIT.tissue_fit;
else TAC_fit = zeros(length(handles.FIT.time))*NaN;end;

if isIF; IF = handles.FIT.IF;
else IF = zeros(length(handles.FIT.time))*NaN;end;

if isIF_fit; IF_fit = handles.FIT.IF_fit;
else IF_fit = zeros(length(handles.FIT.time))*NaN;end;

axes(handles.FittingPlot)
cla
plot(handles.FIT.time,TAC,'*-c'); hold on
plot(handles.FIT.time,IF,'.-m');
plot(handles.FIT.time,TAC_fit,'-b');
plot(handles.FIT.time,IF_fit,'-r'); hold off
xlim([min(handles.FIT.time) max(handles.FIT.time)] )
xlabel('Time(sec)');ylabel('Counts')

guidata(hObject, handles);

function plotResidual(hObject,handles)
axes(handles.ResidualPlot)
cla
plot(handles.FIT.time,zeros(length(handles.FIT.time)),'-r'); hold on
plot(handles.FIT.time,handles.FIT.residual,'*-k');hold off
xlim([min(handles.FIT.time) max(handles.FIT.time)] )
xlabel('Time(sec)');ylabel('Residual')

guidata(hObject, handles);


function [PETData, PETInfo] = read4Ddicom(pathname)

home = pwd;
cd(pathname);
% COPY AND RENAME FILES
directory = dir(pathname); % list all file in directory
old_filename = {directory.name}; % Get the filenames
old_filename = old_filename(3:end);

if old_filename{1}(1)~='1'
    mkdir('temp');
    h = waitbar(0/length(old_filename),['Copying dicom file ', num2str(0),'/',num2str(length(old_filename))]);
    for i=1:length(old_filename)
        waitbar(i/length(old_filename),h,['Copying dicom file ', num2str(i),'/',num2str(length(old_filename))])
        old = char(old_filename(i));
        new = num2str(i);
        copyfile(fullfile(pathname,old),fullfile(pathname,['/temp/',new,'.img']));
    end
    cd('temp')
end

%  DICOM info
start_fnamePET = '1.img';  %utile per rinominare i vari DICOM files
start_PET_info = dicominfo(fullfile(pathname,start_fnamePET)); 
PETSliceTot = double(start_PET_info.NumberOfSlices);
if isfield(start_PET_info, 'NumberOfTimeSlices')
    PETTimeTot = double(start_PET_info.NumberOfTimeSlices);
else
    PETTimeTot = 1;
    h = warndlg('!!! ERROR !!! The selected dataset is not 4D');
end
FrameAcquisitionTime = zeros(PETTimeTot,1); 
FrameStartTime = zeros(PETTimeTot,1);
PETData = zeros(start_PET_info.Height,start_PET_info.Width,PETSliceTot,PETTimeTot);

start_code=str2double(start_fnamePET(1:end-4));

% BUILDING 3D+T DATASET
h = waitbar(0/length(old_filename),['Time frame ', num2str(0),'/',num2str(PETTimeTot)]);
for time = 1:PETTimeTot
    waitbar(time/PETTimeTot,h,['Time frame ', num2str(time),'/',num2str(PETTimeTot)])
    start = (PETSliceTot * time)-1 + start_code;
    stop = start - (PETSliceTot -1);
    idx=1;
    for slice = start : -1 : stop
        suffix = num2str(slice);
        fnamePET = [suffix,'.img'];
        PETFrame = double(dicomread([pathname,'/', fnamePET])); %apro file DICOM
        Inf = dicominfo([pathname ,'/', fnamePET]);  %e il corrispondente header
        PETData(:,:,idx,time)=double(PETFrame*Inf.RescaleSlope + Inf.RescaleIntercept);  %N.B. rescaling dei dati!!!
        idx=idx+1;
    end
    
    FrameAcquisitionTime(time) = double(Inf.ActualFrameDuration/1000);
    FrameStartTime(time) = double(Inf.FrameReferenceTime/1000);
end

Time = [FrameStartTime,FrameStartTime+FrameAcquisitionTime]; % start and end time of each time frame
FOV = double(Inf.ReconstructionDiameter);
PixelSpacing = double(Inf.PixelSpacing);
ImageSize = [double(Inf.Rows), double(Inf.Columns)];
SliceThickness = double(Inf.SliceThickness);

% SAVING DATA INFO
PETInfo = struct('slice_num',PETSliceTot, 'frame_num',PETTimeTot, ...
                 'FrameAcquisitionTime',FrameAcquisitionTime, ...
                 'FrameStartTime', FrameStartTime,'time',Time,'FOV',FOV,...
                 'pixel_size',PixelSpacing,'image_size',ImageSize,...
                 'slice_thickness',SliceThickness);
             
save('PETData','PETData')
save('PETInfo','PETInfo')
cd(home)             
=======
function varargout = KMtool(varargin)
% KMTOOL MATLAB code for KMtool.fig
%      KMTOOL, by itself, creates a new KMTOOL or raises the existing
%      singleton*.
%
%      H = KMTOOL returns the handle to a new KMTOOL or the handle to
%      the existing singleton*.
%
%      KMTOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KMTOOL.M with the given input arguments.
%
%      KMTOOL('Property','Value',...) creates a new KMTOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before KMtool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to KMtool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help KMtool

% Last Modified by GUIDE v2.5 08-Mar-2017 10:06:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @KMtool_OpeningFcn, ...
    'gui_OutputFcn',  @KMtool_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);

if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT


% --- Executes just before KMtool is made visible.
function KMtool_OpeningFcn(hObject, eventdata, handles, varargin)

if isfield(handles,'FIT'); handles = rmfield(handles,'FIT'); end
if isfield(handles,'S'); handles = rmfield(handles,'S'); end
axes(handles.TAC_preview) ; cla
axes(handles.FittingPlot) ; cla
axes(handles.ResidualPlot); cla
set(handles.slice_num,'String',' ');
set(handles.pixel_num,'String',' ');
set(handles.sd_mean  ,'String',' ');
set(handles.ki_value ,'String',' ');
set(handles.k1_value ,'String',' ');
set(handles.k2_value ,'String',' ');
set(handles.k3_value ,'String',' ');
set(handles.k4_value ,'String',' ');
set(handles.vb_value ,'String',' ');
set(handles.sa_value ,'String',' ');
set(handles.Ki_tag   ,'String','Ki')
set(handles.FittingErrorBox,'String',' ')
set(handles.text23   ,'String',' ')

if length(varargin) == 1
    
    handles.immagine = varargin{1};
    handles.clims = [min(handles.immagine(:)) max(handles.immagine(:))];
    
    handles.slices = size(handles.immagine,3);
    handles.frames = size(handles.immagine,4);
    
    set(handles.slice_slider,'max',handles.slices);
    set(handles.slice_slider,'SliderStep',[1/(handles.slices-1) 1/(handles.slices-1)]);
    set(handles.time_slider,'max',handles.frames);
    set(handles.time_slider,'SliderStep',[1/(handles.frames-1) 1/(handles.frames-1)]);
    
    handles.idx_slice=round(handles.slices/2);
    handles.idx_time=round(handles.frames/2);
    set(handles.title,'String',...
        ['Time frame: ',num2str(handles.idx_time),'  |  Slice: ',num2str(handles.idx_slice)]);
    set(handles.time_slider,'Value',handles.idx_time);
    set(handles.slice_slider,'Value',handles.idx_slice);
    set(handles.clim_slider,'min',handles.clims(1)+eps);
    set(handles.clim_slider,'max',handles.clims(2));
    %     set(handles.clim_slider,'SliderStep',[0.1 0.1]);
    set(handles.clim_slider,'Value',handles.clims(2));
    plotMainDataset(handles);
    
elseif length(varargin) > 1
    h = warndlg('!!! ERROR !!! Input data format error');
else
    
end

% Choose default command line output for KMtool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



% UIWAIT makes KMtool wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = KMtool_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;


% ------------------------------ SLIDERS ----------------------------------
function slice_slider_Callback(hObject, eventdata, handles)
handles.idx_slice=round(get(hObject,'Value'));
set(handles.title,'String',['Time frame: ',num2str(handles.idx_time),'  |  Slice: ',num2str(handles.idx_slice)]);
plotMainDataset(handles);
guidata(hObject, handles);
function slice_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject,'min',1);
set(hObject,'max',10);
set(hObject,'SliderStep',[1/(10-1) 1/(10-1)]);
set(hObject,'Value',1);
handles.idx_slice=round(get(hObject,'Value'));
guidata(hObject, handles);


function time_slider_Callback(hObject, eventdata, handles)
handles.idx_time=round(get(hObject,'Value'));
set(handles.title,'String',['Time frame: ',num2str(handles.idx_time),'  |  Slice: ',num2str(handles.idx_slice)]);
plotMainDataset(handles);
guidata(hObject, handles);
function time_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject,'min',1);
set(hObject,'max',10);
set(hObject,'SliderStep',[1/(10-1) 1/(10-1)]);
set(hObject,'Value',1);
handles.idx_time=round(get(hObject,'Value'));
guidata(hObject, handles);

function clim_slider_Callback(hObject, eventdata, handles)
handles.clims(2) = get(hObject,'Value');
plotMainDataset(handles);
guidata(hObject, handles);
function clim_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject,'min',1);
set(hObject,'max',10);
set(hObject,'Value',1);
guidata(hObject, handles);


% ------------------------------- LOAD FILE -------------------------------
function filetype_popup_Callback(hObject, eventdata, handles)
type = get(hObject,'Value');
options = get(hObject,'String');
ft = options(type);
handles.filetype = ft{1};
guidata(hObject, handles);
function filetype_popup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'Value',1);
options = get(hObject,'String');
ft = options(1);
handles.filetype = ft{1};
guidata(hObject, handles);

function load_file_Callback(hObject, eventdata, handles)
[filename, pathname] = uigetfile(['*',handles.filetype], ['Load 4D ',handles.filetype,' image']);

if isfield(handles,'FIT'); handles = rmfield(handles,'FIT'); end
if isfield(handles,'S'); handles = rmfield(handles,'S'); end
axes(handles.TAC_preview) ; cla
axes(handles.FittingPlot) ; cla
axes(handles.ResidualPlot); cla
set(handles.slice_num,'String',' ');
set(handles.pixel_num,'String',' ');
set(handles.sd_mean  ,'String',' ');
set(handles.ki_value ,'String',' ');
set(handles.k1_value ,'String',' ');
set(handles.k2_value ,'String',' ');
set(handles.k3_value ,'String',' ');
set(handles.k4_value ,'String',' ');
set(handles.vb_value ,'String',' ');
set(handles.sa_value ,'String',' ');
set(handles.Ki_tag   ,'String','Ki')
set(handles.FittingErrorBox,'String',' ')
set(handles.text23   ,'String',' ')

if strcmp(handles.filetype,'.mat') 
    tmp = load([pathname,filename]);
    whichVariables = fieldnames(tmp); 
    if numel(whichVariables) == 1
        handles.immagine = tmp.(whichVariables{1});
        handles.clims = [min(handles.immagine(:)) max(handles.immagine(:))];
        handles.slices = size(handles.immagine,3);
        handles.frames = size(handles.immagine,4);

        set(handles.slice_slider,'max',handles.slices);
        set(handles.slice_slider,'SliderStep',[1/(handles.slices-1) 1/(handles.slices-1)]);
        set(handles.time_slider,'max',handles.frames);
        set(handles.time_slider,'SliderStep',[1/(handles.frames-1) 1/(handles.frames-1)]);


        handles.idx_slice=round(handles.slices/2);
        handles.idx_time=round(handles.frames/2);
        set(handles.title,'String',...
            ['Time frame: ',num2str(handles.idx_time),'  |  Slice: ',num2str(handles.idx_slice)]);
        set(handles.time_slider,'Value',handles.idx_time);
        set(handles.slice_slider,'Value',handles.idx_slice);

        set(handles.clim_slider,'min',handles.clims(1)+eps);
        set(handles.clim_slider,'max',handles.clims(2));
        set(handles.clim_slider,'Value',handles.clims(2));

        plotMainDataset(handles);

        h = warndlg('Image data loaded!');
    else
        h = warndlg('!!! ERROR !!! Input file format error');
    end
    
elseif strcmp(handles.filetype,'.img') || strcmp(handles.filetype,'.dcm')
    [PETData, PETInfo] = read4Ddicom(pathname);
    handles.immagine = PETData;
    handles.S.time_vec= [PETInfo.time(1,1); PETInfo.time(:,2)]; % mean(Info.time,2);
    handles.S.scant = PETInfo.time;
    handles.clims = [min(handles.immagine(:)) max(handles.immagine(:))];
    handles.slices = size(handles.immagine,3);
    handles.frames = size(handles.immagine,4);

    set(handles.slice_slider,'max',handles.slices);
    set(handles.slice_slider,'SliderStep',[1/(handles.slices-1) 1/(handles.slices-1)]);
    set(handles.time_slider,'max',handles.frames);
    set(handles.time_slider,'SliderStep',[1/(handles.frames-1) 1/(handles.frames-1)]);


    handles.idx_slice=round(handles.slices/2);
    handles.idx_time=round(handles.frames/2);
    set(handles.title,'String',...
        ['Time frame: ',num2str(handles.idx_time),'  |  Slice: ',num2str(handles.idx_slice)]);
    set(handles.time_slider,'Value',handles.idx_time);
    set(handles.slice_slider,'Value',handles.idx_slice);

    set(handles.clim_slider,'min',handles.clims(1)+eps);
    set(handles.clim_slider,'max',handles.clims(2));
    set(handles.clim_slider,'Value',handles.clims(2));

    plotMainDataset(handles);

    h = warndlg('Image data loaded!');
end

guidata(hObject, handles);


% --- Executes on button press in transpose_btn.
function transpose_btn_Callback(hObject, eventdata, handles)
handles.immagine = permute(handles.immagine,[2,1,3,4]);
plotMainDataset(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function title_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','Time frame: --  |  Slice: --');
guidata(hObject, handles);

% ---------------------------- COLORMAP -----------------------------------
function colormap_Callback(hObject, eventdata, handles)
values = cellstr(get(hObject,'String'));
handles.cmap=values{get(hObject,'Value')};
plotMainDataset(handles);
guidata(hObject, handles);

function colormap_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
values = cellstr(get(hObject,'String'));
set(hObject,'Value',7);
handles.cmap=values{get(hObject,'Value')};
guidata(hObject, handles);


% ------------------------- LOAD TIME VECTOR ------------------------------
function LoadTime_btn_Callback(hObject, eventdata, handles)

[filename, pathname] = uigetfile('*.mat', 'Load PETInfo file');

tmp = load([pathname,filename]);
whichVariables = fieldnames(tmp);
if numel(whichVariables) == 1
    Info = tmp.(whichVariables{1});
    handles.S.scant = Info.time;
    handles.S.time_vec= mean(Info.time,2); % [Info.time(1,1); Info.time(:,2)];
    
    h = warndlg('Time vector loaded!') ;
    
    if isfield(handles.S,'TAC')
        axes(handles.TAC_preview)
        cla
        plot(handles.S.time_vec,handles.S.TAC,'*-');
        xlabel('Time(sec)');ylabel('Counts')
        xlim([min(handles.S.time_vec) max(handles.S.time_vec)] )
    end
    
else
    h = warndlg('!!! ERROR !!! Input file format error');
end
guidata(hObject, handles);

% ------------------------- ROI SELECTION ---------------------------------
function roi_btn_Callback(hObject, eventdata, handles)

%plotMainDataset(handles);
handles = externalFigureRoiSelection(hObject,handles);

h = imfreehand; %position = wait(h);
ROI = h.createMask();

k = 0;
for i = 1:size(ROI,1)
    for j = 1:size(ROI,2)
        if ROI(i,j)==0
            continue
        else
            k=k+1;
            TAC(k,:) = squeeze(handles.immagine(i,j,handles.idx_slice,:));
        end
    end
end

set(handles.slice_num,'String',num2str(handles.idx_slice));
set(handles.pixel_num,'String',num2str(size(TAC,1)));
%set(handles.sd_mean,'String',num2str(mean(std(TAC)./mean(TAC),'omitnan')));

handles.S.TAC   = mean(TAC)'; %[0; mean(TAC)'];
handles.S.slice = handles.idx_slice;
handles.S.ROI   = ROI;

axes(handles.TAC_preview)
cla
if ~isfield(handles.S,'time_vec')
    plot(handles.S.TAC,'*-');
    xlabel('Frame number');ylabel('Counts')
    xlim([0 length(handles.S.TAC)])
else
    plot(handles.S.time_vec,handles.S.TAC,'*-');
    xlabel('Time(sec)');ylabel('Counts')
    xlim([min(handles.S.time_vec) max(handles.S.time_vec)] )
end

close(handles.H)
guidata(hObject, handles);

% ---------------------- SEND TACs TO MODELING TOOL -----------------------

function setIF_btn_Callback(hObject, eventdata, handles)

if ~isempty(handles.S)
    if isfield(handles.S,'time_vec')
        handles.FIT.IF     = double(handles.S.TAC);
        handles.FIT.time   = handles.S.time_vec;
        handles.FIT.scant  = handles.S.scant;
        
        if isfield(handles.FIT,'IF_fit')
            handles.FIT = rmfield(handles.FIT,'IF_fit');
            handles.FIT = rmfield(handles.FIT,'IF_params');
            set(handles.IFfit_chk,'Value',0.0);
        end
        set(handles.IFmeas_chk,'Value',1.0);        
        plotFIT(hObject, handles);
    else
        h = warndlg('You are currently using a uniform time vector!') ;
        handles.FIT.IF     = double(handles.S.TAC);
        handles.FIT.scant  = [linspace(0,length(handles.S.TAC)-1,length(handles.S.TAC)); ...
                                linspace(1,length(handles.S.TAC),length(handles.S.TAC))]';
        handles.FIT.time   = mean(handles.FIT.scant,2);
        if isfield(handles.FIT,'IF_fit')
            handles.FIT = rmfield(handles.FIT,'IF_fit');
            handles.FIT = rmfield(handles.FIT,'IF_params');
            set(handles.IFfit_chk,'Value',0.0);
        end
        set(handles.IFmeas_chk,'Value',1.0);        
        plotFIT(hObject, handles);
    end
else
    h = warndlg('You need to draw a ROI first!') ;
end
guidata(hObject, handles);

% --- Executes on button press in setTAC_btn.
function setTAC_btn_Callback(hObject, eventdata, handles)

if ~isempty(handles.S)
    if isfield(handles.S,'time_vec')
        handles.FIT.tissue = double(handles.S.TAC);
        handles.FIT.time   = handles.S.time_vec;
        handles.FIT.scant  = handles.S.scant;
        disp(max(handles.S.scant(:)))
        
        if isfield(handles.FIT,'tissue_fit')
            handles.FIT = rmfield(handles.FIT,'tissue_fit');
            handles.FIT = rmfield(handles.FIT,'tissue_params');
            handles.FIT = rmfield(handles.FIT,'residual');
            set(handles.TACfit_chk,'Value',0.0);
            axes(handles.ResidualPlot)
            cla
        end
        set(handles.TACmeas_chk,'Value',1.0);
        plotFIT(hObject, handles);
    else
        h = warndlg('You are currently using a uniform time vector!') ;
        handles.FIT.tissue = double(handles.S.TAC);
        handles.FIT.scant  = [linspace(0,length(handles.S.TAC)-1,length(handles.S.TAC)); ...
                                linspace(1,length(handles.S.TAC),length(handles.S.TAC))]';
        handles.FIT.time   = mean(handles.FIT.scant,2);
        if isfield(handles.FIT,'tissue_fit')
            handles.FIT = rmfield(handles.FIT,'tissue_fit');
            handles.FIT = rmfield(handles.FIT,'tissue_params');
            handles.FIT = rmfield(handles.FIT,'residual');
            set(handles.TACfit_chk,'Value',0.0);
            axes(handles.ResidualPlot)
            cla
        end
        set(handles.TACmeas_chk,'Value',1.0);
        plotFIT(hObject, handles);
    end
else
    h = warndlg('You need to draw a ROI first!') ;
end
guidata(hObject, handles);




% ------------------------ EXPORT DATA FROM GUI ---------------------------
function saveROI_btn_Callback(hObject, eventdata, handles)
if isfield(handles,'S')
    ROI = handles.S.ROI;
    slice_num = handles.S.slice;
    TAC = handles.S.TAC;
    uisave({'ROI','slice_num','TAC'},'ROI.mat')
end
guidata(hObject, handles);

function saveFIT_btn_Callback(hObject, eventdata, handles)
if isfield(handles,'FIT')
    FIT = handles.FIT;
    uisave({'FIT'},'FIT.mat')
end
guidata(hObject, handles);

% ------------------------------ CHECKBUTTONS -----------------------------
function IFmeas_chk_Callback(hObject, eventdata, handles)
plotFIT(hObject, handles);guidata(hObject, handles);

function IFfit_chk_Callback(hObject, eventdata, handles)
plotFIT(hObject, handles);guidata(hObject, handles);

function TACmeas_chk_Callback(hObject, eventdata, handles)
plotFIT(hObject, handles);guidata(hObject, handles);

function TACfit_chk_Callback(hObject, eventdata, handles)
plotFIT(hObject, handles);guidata(hObject, handles);

% ---------------------------TAC MODEL SELECTION --------------------------
function TACmodels_switch_Callback(hObject, eventdata, handles)

values = cellstr(get(hObject,'String'));
handles.TACmodel=values{get(hObject,'Value')};

if isfield(handles.FIT,'IF_fit')
    IF = handles.FIT.IF_fit;
else
    IF = handles.FIT.IF;
end

switch handles.TACmodel
    
    case {'DCE_MRI 1TC','DCE_MRI 2TCr','DCE_MRI 2TCi'}
        [fit, params, resnorm, residual, output] = fit_dce_models (handles.FIT.tissue,handles.FIT.scant,IF,handles.TACmodel);
        handles.FIT.tissue_fit = fit;
        handles.FIT.tissue_params = params;
        handles.FIT.residual = residual;

        msg1 = struct2cell(output);
        msg2 = fieldnames(output);
        set(handles.FittingErrorBox,'String',msg1(1:6))
        set(handles.text23,'String',msg2(1:6))
        set(handles.TACfit_chk,'Value',1.0);
                    
        plotFIT(hObject, handles);
        plotResidual(hObject,handles)
        
        set(handles.ki_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(7),3)));
        set(handles.k1_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(3),3)));
        set(handles.k2_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(4),3)));
        set(handles.k3_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(5),3)));
        set(handles.k4_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(6),3)));
        set(handles.vb_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(1),3)));
        set(handles.sa_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(2),3)));
        set(handles.Ki_tag,'String','Ki')
    
    case 'Logan'
        
        if isfield(handles.FIT,'residual')
            handles.FIT = rmfield(handles.FIT,'tissue_fit');
            handles.FIT = rmfield(handles.FIT,'tissue_params');
            handles.FIT = rmfield(handles.FIT,'residual');
            axes(handles.ResidualPlot)
            cla
        end
        
        res = logan_plot(IF, handles.FIT.tissue, handles.FIT.scant);
        
        axes(handles.FittingPlot)
        cla
        plot(res.termx, res.termy,'LineStyle','none',...
            'Marker','o','MarkerEdgeColor','black'); hold on;
        plot(res.termx(res.t0:end), res.termy(res.t0:end),'LineStyle','none',...
            'Marker','o','MarkerEdgeColor','black',...
            'MarkerFaceColor','black')
        plot(res.termx, res.kparms(1) + res.kparms(2)*res.termx,'r-')
        xlim([min(res.termx) max(res.termx)])
        ylim([min(res.termy) max(res.termy)])
        xlabel('Int(Cplasma) / Ctissue [min]')
        ylabel('Int(Ctissue) / Ctissue [min]')
        
        handles.FIT.tissue_fit = res.fit;
        handles.FIT.tissue_params = res.kparms;
        
        set(handles.ki_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(2),3)));
        set(handles.k1_value,'String',' ');
        set(handles.k2_value,'String',' ');
        set(handles.k3_value,'String',' ');
        set(handles.k4_value,'String',' ');
        set(handles.vb_value,'String',' ');
        set(handles.sa_value,'String',' ');
        set(handles.Ki_tag,'String','DV')
        
        
    case 'Patlak'
                
        if isfield(handles.FIT,'residual')
            handles.FIT = rmfield(handles.FIT,'tissue_fit');
            handles.FIT = rmfield(handles.FIT,'tissue_params');
            handles.FIT = rmfield(handles.FIT,'residual');
            axes(handles.ResidualPlot)
            cla
        end
        
        res = patlak_plot(IF, handles.FIT.tissue, handles.FIT.scant);
        
        axes(handles.FittingPlot)
        cla
        plot(res.termx, res.termy,'LineStyle','none',...
            'Marker','o','MarkerEdgeColor','black'); hold on;
        plot(res.termx(res.t0:end), res.termy(res.t0:end),'LineStyle','none',...
            'Marker','o','MarkerEdgeColor','black',...
            'MarkerFaceColor','black')
        plot(res.termx, res.kparms(1) + res.kparms(2)*res.termx,'r-')
        xlim([min(res.termx) max(res.termx)])
        ylim([min(res.termy) max(res.termy)])
        xlabel('Int(Cplasma) / Ctissue [min]')
        ylabel('Int(Ctissue) / Ctissue [min]')
        
        handles.FIT.tissue_fit = res.fit;
        handles.FIT.tissue_params = res.kparms;
        
        set(handles.ki_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(2),3)));
        set(handles.k1_value,'String',' ');
        set(handles.k2_value,'String',' ');
        set(handles.k3_value,'String',' ');
        set(handles.k4_value,'String',' ');
        set(handles.vb_value,'String',' ');
        set(handles.sa_value,'String',' ');
        set(handles.Ki_tag,'String','Ki')
        
    case 'Tissue modeling'
        return;
        
    case {'2TCr_analytic','2TCr_analytic_nodelay'}
        [fit, params, resnorm, residual, output] = fit_analytic_models (handles.FIT.tissue,handles.FIT.scant,handles.FIT.IF_params,handles.FIT.IF_fit, handles.TACmodel);
        handles.FIT.tissue_fit = fit;
        handles.FIT.tissue_params = params;
        handles.FIT.residual = residual;

        msg1 = struct2cell(output);
        msg2 = fieldnames(output);
        set(handles.FittingErrorBox,'String',msg1(1:6))
        set(handles.text23,'String',msg2(1:6))
        set(handles.TACfit_chk,'Value',1.0);
                    
        plotFIT(hObject, handles);
        plotResidual(hObject,handles)
        
        set(handles.ki_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(7),3)));
        set(handles.k1_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(3),3)));
        set(handles.k2_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(4),3)));
        set(handles.k3_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(5),3)));
        set(handles.k4_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(6),3)));
        set(handles.vb_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(1),3)));
        set(handles.sa_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(2),3)));
        set(handles.Ki_tag,'String','Ki')
        
    case {'DCE-MRI an_2TCr','DCE-MRI an_2TCi'}
        disp(max(handles.FIT.scant(:)))
        [fit, params, resnorm, residual, output] = fit_dce_analytic_models (handles.FIT.tissue,handles.FIT.scant,handles.FIT.IF_params,handles.FIT.IF_fit, handles.TACmodel);
        handles.FIT.tissue_fit = fit;
        handles.FIT.tissue_params = params;
        handles.FIT.residual = residual;

        msg1 = struct2cell(output);
        msg2 = fieldnames(output);
        set(handles.FittingErrorBox,'String',msg1(1:6))
        set(handles.text23,'String',msg2(1:6))
        set(handles.TACfit_chk,'Value',1.0);
                    
        plotFIT(hObject, handles);
        plotResidual(hObject,handles)
        
        set(handles.ki_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(7),3)));
        set(handles.k1_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(3),3)));
        set(handles.k2_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(4),3)));
        set(handles.k3_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(5),3)));
        set(handles.k4_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(6),3)));
        set(handles.vb_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(1),3)));
        set(handles.sa_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(2),3)));
        set(handles.Ki_tag,'String','Ki')
        
    otherwise
        [fit, params, resnorm, residual, output] = fit_compartmental_models (handles.FIT.tissue,handles.FIT.scant,IF,handles.TACmodel);
        handles.FIT.tissue_fit = fit;
        handles.FIT.tissue_params = params;
        handles.FIT.residual = residual;

        msg1 = struct2cell(output);
        msg2 = fieldnames(output);
        set(handles.FittingErrorBox,'String',msg1(1:6))
        set(handles.text23,'String',msg2(1:6))
        set(handles.TACfit_chk,'Value',1.0);
                    
        plotFIT(hObject, handles);
        plotResidual(hObject,handles)

        set(handles.ki_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(7),3)));
        set(handles.k1_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(3),3)));
        set(handles.k2_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(4),3)));
        set(handles.k3_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(5),3)));
        set(handles.k4_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(6),3)));
        set(handles.vb_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(1),3)));
        set(handles.sa_value,'String',num2str(optimizeParamVisualiztion(handles.FIT.tissue_params(2),3)));
        set(handles.Ki_tag,'String','Ki')
end

guidata(hObject, handles);


function TACmodels_switch_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'Value',1);


% --------------------------- IF MODEL SELECTION --------------------------
function IFmodels_switch_Callback(hObject, eventdata, handles)
handles.IFmodel=get(hObject,'Value');

switch handles.IFmodel
    case 1
        return;
    otherwise
        [Cp, ppCp, ifParams] = fittingFengInputModel2(handles.FIT.IF, handles.FIT.time, handles.IFmodel);
        handles.FIT.IF_fit = Cp;
        handles.FIT.IF_params = ifParams;
end
isTissue = isfield(handles.FIT,'tissue');
isIF = isfield(handles.FIT,'IF');
isTissue_fit = isfield(handles.FIT,'tissue_fit');
isIF_fit = isfield(handles.FIT,'IF_fit');
checks = [isTissue, isIF, isTissue_fit, isIF_fit];
if isTissue;     set(handles.TACmeas_chk,'Value',1);end
if isIF;         set(handles.IFmeas_chk,'Value',1);end
if isTissue_fit; set(handles.TACfit_chk,'Value',1);end
if isIF_fit;     set(handles.IFfit_chk,'Value',1);end
plotFIT(hObject, handles);

guidata(hObject, handles);

function IFmodels_switch_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'Value',1);


% ---------------------------------------------------------------------------- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function k1_value_Callback(hObject, eventdata, handles)
function k1_value_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function vb_value_Callback(hObject, eventdata, handles)
function vb_value_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function k4_value_Callback(hObject, eventdata, handles)
function k4_value_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function k3_value_Callback(hObject, eventdata, handles)
function k3_value_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function k2_value_Callback(hObject, eventdata, handles)
function k2_value_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sa_value_Callback(hObject, eventdata, handles)
function sa_value_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ki_value_Callback(hObject, eventdata, handles)
function ki_value_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function slice_num_Callback(hObject, eventdata, handles)
function slice_num_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pixel_num_Callback(hObject, eventdata, handles)
function pixel_num_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sd_mean_Callback(hObject, eventdata, handles)
function sd_mean_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% -------------------------- UTILITY FUNCTION -----------------------------

function ret = optimizeParamVisualiztion(n,dec)
f = 10.^dec;
ret = round(f*n)/f;

function plotMainDataset(handles)
axes(handles.plot_axes)
imagesc(squeeze(handles.immagine(:,:,handles.idx_slice,handles.idx_time)),handles.clims),colormap(handles.cmap);

function handles = externalFigureRoiSelection(hObject, handles)
handles.H = figure(1);
imagesc(squeeze(handles.immagine(:,:,handles.idx_slice,handles.idx_time)),handles.clims),colormap(handles.cmap),axis equal tight;
guidata(hObject, handles);


function plotFIT (hObject, handles)

isTissue = isfield(handles.FIT,'tissue') && get(handles.TACmeas_chk,'Value');
isIF = isfield(handles.FIT,'IF') && get(handles.IFmeas_chk,'Value');
isTissue_fit = isfield(handles.FIT,'tissue_fit') && get(handles.TACfit_chk,'Value');
isIF_fit = isfield(handles.FIT,'IF_fit') && get(handles.IFfit_chk,'Value');


if isTissue;     set(handles.TACmeas_chk,'Value',1);end
if isIF;         set(handles.IFmeas_chk,'Value',1);end
if isTissue_fit; set(handles.TACfit_chk,'Value',1);end
if isIF_fit;     set(handles.IFfit_chk,'Value',1);end

if isTissue; TAC = handles.FIT.tissue;
else TAC = zeros(length(handles.FIT.time))*NaN;end;

if isTissue_fit; TAC_fit = handles.FIT.tissue_fit;
else TAC_fit = zeros(length(handles.FIT.time))*NaN;end;

if isIF; IF = handles.FIT.IF;
else IF = zeros(length(handles.FIT.time))*NaN;end;

if isIF_fit; IF_fit = handles.FIT.IF_fit;
else IF_fit = zeros(length(handles.FIT.time))*NaN;end;

axes(handles.FittingPlot)
cla
plot(handles.FIT.time,TAC,'*-c'); hold on
plot(handles.FIT.time,IF,'.-m');
plot(handles.FIT.time,TAC_fit,'-b');
plot(handles.FIT.time,IF_fit,'-r'); hold off
xlim([min(handles.FIT.time) max(handles.FIT.time)] )
xlabel('Time(sec)');ylabel('Counts')

guidata(hObject, handles);

function plotResidual(hObject,handles)
axes(handles.ResidualPlot)
cla
plot(handles.FIT.time,zeros(length(handles.FIT.time)),'-r'); hold on
plot(handles.FIT.time,handles.FIT.residual,'*-k');hold off
xlim([min(handles.FIT.time) max(handles.FIT.time)] )
xlabel('Time(sec)');ylabel('Residual')

guidata(hObject, handles);


function [PETData, PETInfo] = read4Ddicom(pathname)

home = pwd;
cd(pathname);
% COPY AND RENAME FILES
directory = dir(pathname); % list all file in directory
old_filename = {directory.name}; % Get the filenames
old_filename = old_filename(3:end);

if old_filename{1}(1)~='1'
    mkdir('temp');
    h = waitbar(0/length(old_filename),['Copying dicom file ', num2str(0),'/',num2str(length(old_filename))]);
    for i=1:length(old_filename)
        waitbar(i/length(old_filename),h,['Copying dicom file ', num2str(i),'/',num2str(length(old_filename))])
        old = char(old_filename(i));
        new = num2str(i);
        copyfile(fullfile(pathname,old),fullfile(pathname,['/temp/',new,'.img']));
    end
    cd('temp')
end

%  DICOM info
start_fnamePET = '1.img';  %utile per rinominare i vari DICOM files
start_PET_info = dicominfo(fullfile(pathname,start_fnamePET)); 
PETSliceTot = double(start_PET_info.NumberOfSlices);
if isfield(start_PET_info, 'NumberOfTimeSlices')
    PETTimeTot = double(start_PET_info.NumberOfTimeSlices);
else
    PETTimeTot = 1;
    h = warndlg('!!! ERROR !!! The selected dataset is not 4D');
end
FrameAcquisitionTime = zeros(PETTimeTot,1); 
FrameStartTime = zeros(PETTimeTot,1);
PETData = zeros(start_PET_info.Height,start_PET_info.Width,PETSliceTot,PETTimeTot);

start_code=str2double(start_fnamePET(1:end-4));

% BUILDING 3D+T DATASET
h = waitbar(0/length(old_filename),['Time frame ', num2str(0),'/',num2str(PETTimeTot)]);
for time = 1:PETTimeTot
    waitbar(time/PETTimeTot,h,['Time frame ', num2str(time),'/',num2str(PETTimeTot)])
    start = (PETSliceTot * time)-1 + start_code;
    stop = start - (PETSliceTot -1);
    idx=1;
    for slice = start : -1 : stop
        suffix = num2str(slice);
        fnamePET = [suffix,'.img'];
        PETFrame = double(dicomread([pathname,'/', fnamePET])); %apro file DICOM
        Inf = dicominfo([pathname ,'/', fnamePET]);  %e il corrispondente header
        PETData(:,:,idx,time)=double(PETFrame*Inf.RescaleSlope + Inf.RescaleIntercept);  %N.B. rescaling dei dati!!!
        idx=idx+1;
    end
    
    FrameAcquisitionTime(time) = double(Inf.ActualFrameDuration/1000);
    FrameStartTime(time) = double(Inf.FrameReferenceTime/1000);
end

Time = [FrameStartTime,FrameStartTime+FrameAcquisitionTime]; % start and end time of each time frame
FOV = double(Inf.ReconstructionDiameter);
PixelSpacing = double(Inf.PixelSpacing);
ImageSize = [double(Inf.Rows), double(Inf.Columns)];
SliceThickness = double(Inf.SliceThickness);

% SAVING DATA INFO
PETInfo = struct('slice_num',PETSliceTot, 'frame_num',PETTimeTot, ...
                 'FrameAcquisitionTime',FrameAcquisitionTime, ...
                 'FrameStartTime', FrameStartTime,'time',Time,'FOV',FOV,...
                 'pixel_size',PixelSpacing,'image_size',ImageSize,...
                 'slice_thickness',SliceThickness);
             
save('PETData','PETData')
save('PETInfo','PETInfo')
cd(home)             
>>>>>>> 5fa9534adcf332ff010dadcf5e9217186ab53ed8
