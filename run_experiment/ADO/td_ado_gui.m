function varargout = td_ado_gui(varargin)
    %TD_ADO_GUI MATLAB code file for td_ado_gui.fig
    %      TD_ADO_GUI, by itself, creates a new TD_ADO_GUI or raises the existing
    %      singleton*.
    %
    %      H = TD_ADO_GUI returns the handle to a new TD_ADO_GUI or the handle to
    %      the existing singleton*.
    %
    %      TD_ADO_GUI('Property','Value',...) creates a new TD_ADO_GUI using the
    %      given property value pairs. Unrecognized properties are passed via
    %      varargin to td_ado_gui_OpeningFcn.  This calling syntax produces a
    %      warning when there is an existing singleton*.
    %
    %      TD_ADO_GUI('CALLBACK') and TD_ADO_GUI('CALLBACK',hObject,...) call the
    %      local function named CALLBACK in TD_ADO_GUI.M with the given input
    %      arguments.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % THIS SCRIPT IS SUPPOSE TO BE INITIALIZED FROM ../td_demo.m DO NOT RUN
    % FORM HERE!

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @td_ado_gui_OpeningFcn, ...
                       'gui_OutputFcn',  @td_ado_gui_OutputFcn, ...
                       'gui_LayoutFcn',  [], ...
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
end



% --- Executes just before td_ado_gui is made visible.
function td_ado_gui_OpeningFcn(hObject, eventdata, handles, varargin)
    global rep id staging save_dir
    
    % import libraries
    addpath(genpath('../../GPStuff-4.7'));
    addpath(genpath('../ADO'));
    rep=1; % repeat this condition after finishing will be set to -1 automatically after one repetition
    machine='A'; % computer identifier
   
    staging=true; % represents whether this is the first condition to be run
    save_dir='data';

    % Give participant an id that hasn't been used
    id=1; 
    exists = true;
    while exists
        name= strcat(save_dir,"/",num2str(id),machine);
        if exist(strcat(name,'_1'),'file')==7 || exist(strcat(name,'_2'),'file')==7
            id=id+1;
        else
            exists=false;
        end
    end
    id = [num2str(id),machine,'_1'];
    initialize_experiment(handles);
    if ~staging
        setup_experiment(handles);
        set(handles.pushbutton3, 'visible','on');
    end

    % Choose default command line output for td_ado_gui
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);
end

function initialize_experiment(handles)
    global yal xal experiments trial
    experiments = 25;
    yal= [];
    xal = [];
    trial = 0;
    debug =false;
    if ~debug
        set(handles.text4,'visible','off');
    end
    %rng(80);
end

function varargout = td_ado_gui_OutputFcn(hObject, eventdata, handles)
    varargout{1} = handles.output;
end

% --- Run ADO
function next_design = run_experiment(y,handles)
    global trial id save_dir
    if trial == 0
        vVar = td(trial,1,id,save_dir);
    else
        vVar = td(trial,y(end),id,save_dir);
    end
    next_design = vVar.design;
end

% --- Disable some buttons while handeling data
function pushbutton2_Callback(hObject, eventdata, handles)
    global trial yal b2
    set(handles.pushbutton2,'Enable','off');
    set(handles.pushbutton3,'Enable','off');
    if trial~=0
        yal = [yal; b2];
    else
        set(handles.pushbutton3, 'visible','on');
    end
    choose_experiment(handles);
end

% --- Decide what to do next
function choose_experiment(handles)
    global trial yal xal experiments rep id staging;
    if trial<experiments
        setup_experiment(handles);
    else
        if rep==1
            rep=-1;
            id = [id(1:end-2),'_2'];
            initialize_experiment(handles);
            setup_experiment(handles);
        else
            if staging
                close(td_ado_gui);
                addpath('..');
                run('td_gui');
            else
                set(handles.pushbutton2,'visible','off');
                set(handles.pushbutton3,'visible','off');
                set(handles.text6,'visible','on');
            end
        end
    end
end

% --- process user choice
function pushbutton3_Callback(hObject, eventdata, handles)
    global b3 yal

    set(handles.pushbutton2,'Enable','off');
    set(handles.pushbutton3,'Enable','off');
    yal = [yal; b3];
    choose_experiment(handles);
end

% --- Show Instructions
function pushbutton4_Callback(hObject, eventdata, handles)
    f = warndlg(['In this study you will be asked about your preferences regarding different amounts of money'...
    ' to be delivered at different times in the future.  You will have several seconds to decide between the '...
    'two options provided on the screen, so there is no rush to respond. If you prefer the option on the left.'...
    'you should press the left button to indicate your choice. If you prefer the option on the right, please '...
    'press the right button. After you have responded. There will then be a brief pause before the next pair of '...
    'alternatives is presented. You may be asked to choose between the same set of options more than once. '...
    'There are no wrong answers in this experiment. Please indicate your preference as honestly as you can, '...
    'and make your best guess when you are not sure.'],'Instructions');
end

% --- Update experiment prompt
function setup_experiment(handles)
    global trial yal xal b2 b3

    h = waitbar(0,'Loading, please wait (this may take up to 20 seconds)');
    next_design = run_experiment(yal,handles);
    close(h)
    xal = [xal; next_design];

    trial=trial+1;
    time_string = parse_time(next_design(4));
    money_string = sprintf('$%0.0f',next_design(1));
    rnum = randi(2);
    if rnum==1
        set(handles.pushbutton2,'string',[money_string, ' Now']);
        set(handles.pushbutton3,'string',['$800 in ', time_string]);
        b2=0;
        b3=1;
    else
        set(handles.pushbutton3,'string',[money_string, ' Now']);
        set(handles.pushbutton2,'string',['$800 in ', time_string]);
        b2=1;
        b3=0;
    end
    set(handles.text4,'string',['Trial: ', num2str(trial)]);
    set(handles.pushbutton2,'Enable','on');
    set(handles.pushbutton3,'Enable','on');
end

% --- Convert time from number to text
function time = parse_time(logtime)
    days = logtime*7;
    if days>364
       years = round(days/365);
       if years~=1
           time = [num2str(years),' years'];
       else
           time = [num2str(years),' year'];
       end

       if abs(years-days/365)>0.1
           time = [time, ' and a half'];
       end
    elseif days<21,
       if days==1
           time = [num2str(round(days)), ' day'];
       else
            time = [num2str(round(days)), ' days'];
       end
    elseif days<59
        weeks = round(days/7);
        time = [num2str(weeks),' weeks'];
        if abs(weeks-days/7)>0.1
           time = [time, ' and a half'];
        end
    else
        months = round(days/30);
        time = [num2str(months),' months'];
        if abs(months-days/30)>0.1
            time = [time, ' and a half'];
        end
    end
end
