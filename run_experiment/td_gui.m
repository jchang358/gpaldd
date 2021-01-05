function varargout = td_gui(varargin)
    % TD_GUI MATLAB code for td_gui.fig
    %      TD_GUI, by itself, creates a new TD_GUI or raises the existing
    %      singleton*.
    %
    %      H = TD_GUI returns the handle to a new TD_GUI or the handle to
    %      the existing singleton*.
    %
    %      TD_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in TD_GUI.M with the given input arguments.
    %
    %      TD_GUI('Property','Value',...) creates a new TD_GUI or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before td_gui_OpeningFcn gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to td_gui_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @td_gui_OpeningFcn, ...
                       'gui_OutputFcn',  @td_gui_OutputFcn, ...
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
end

% --- Executes just before td_gui is made visible.
function td_gui_OpeningFcn(hObject, eventdata, handles, varargin)
    global rep machine staging
    addpath(genpath('../GPStuff-4.7'));
    rep=1;
    machine='D';
    staging=false;
    initialize_experiment(handles);
    setup_experiment(handles);
    set(handles.pushbutton3, 'visible','on');
    % Choose default command line output for td_gui
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);
end

% --- Initialize GPAL
function initialize_experiment(handles)
    global gp opt va_t m experiments trial invTran
    global xseed yseed xt xopt xal yal xt1 xt2 debug Tran
    %initialize all hyperparameters
    Tran = @(x) log10(x);
    debug = false;
    invTran = @(x) 10.^x;
    trial = 0;
    experiments = 40;
    
    %grid bounds
    lb = [0.143 10];
    ub = [520 790];
    
    %margin init
    va_t = 0.5;
    m = -1;
    nseed = 15;
    
    %rng(80);
    
    if ~debug
        set(handles.text4,'visible','off');
    end
    ASS_seed =linspace(Tran(lb(2)),Tran(ub(2)),nseed);
    TLL_seed = linspace(Tran(lb(1)),Tran(ub(1)),nseed);
    TLL_space = [  0.143,  0.286,   0.429,  0.571,  0.714,   0.857,    1,...
                   1.143,  1.286,   1.428,  1.571,  1.714,   1.857,    2,...
                   2.286,  2.571,   3,      3.5,    4,       4.5,      5,...
                   6,      7,       8.57,   10.714, 12.86,   15,       17.14,...
                   21.43,  25.714,  30,     34.286, 38.57,   42.86,    47.14,...
                   52.14,  78.214,  104.286,130.36, 156.43,  182.5,    208.57,...
                   260.714,312.86,  365,    417.14, 469.286, 521.43];

    ASS_space =[lb(2):10:ub(2)]; 

    ASS_fullgrid =linspace(lb(2),ub(2),500);
    TLL_fullgrid = linspace(lb(1),ub(1),500);
    [xt1, xt2] = meshgrid((Tran(TLL_fullgrid))', (Tran(ASS_fullgrid))');

    [xopt1, xopt2] = meshgrid((Tran(TLL_space))', (Tran(ASS_space))');

    [x1seed, x2seed] = meshgrid(TLL_seed', ASS_seed');
    xseed = [x1seed(:) x2seed(:)];
    yseed = 2*(xseed(:,1)+xseed(:,2)<4)-1;

    bt = va_t*sqrt(m^2+1)+2-m*2;
    bu = 2-m*2-va_t*sqrt(m^2+1);
    margin = (0<m*xseed(:,1)+bt-xseed(:,2)) & (0>m*xseed(:,1)+bu-xseed(:,2));
    
    xseed(margin,:)=[];
    yseed(margin)=[];
    
    % Test data
    xt=[xt1(:) xt2(:)];
    xopt = [xopt1(:) xopt2(:)];
    % Create likelihood function
    lik = lik_probit();

    % Create covariance functions
    gpcf = gpcf_sexp('lengthScale', [0.9 0.9], 'magnSigma2', 10);

    %initialize GP
    lik = lik_probit();
    pl = prior_t();
    pm = prior_loggaussian('mu',0,'s2',10);
    gpcf = gpcf_sexp('lengthScale', [0.9 0.9], 'magnSigma2', 10);
    gpcf = gpcf_sexp(gpcf, 'lengthScale_prior', pl,'magnSigma2_prior', pm); %

    % Create the GP structure (type is by default FULL)
    gp = gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-9);
    % Set the approximate inference method 
    gp = gp_set(gp, 'latent_method', 'EP');
    gp.code=[1];
    gp.monotonicity=0;
    opt=optimset('TolFun',1e-3,'TolX',1e-3);
    xal = [];
    yal = [];
end

% --- Run GPAL
function next_design = run_experiment(x,y,handles)
    global gp opt xopt trial xt xt1 xt2 debug invTran xal yal rep Tran staging experiments
    Cent = sqrt(pi*log(2)/2);
    Ent = @(p) -p.*log(p)-(1-p).*log(1-p);

    gp=gp_optim(gp,x,y,'opt',opt);

    [Eft_la, Varft_la, lpyt_la, Eyt_la, Varyt_la] = ...
            gp_pred(gp, x, y, xopt, 'yt', ones(size(xopt,1),1) );

    efp =(Cent*exp(-Eyt_la./(2*(Varyt_la+Cent^2))))./sqrt(Varyt_la+Cent^2);
    Entropy = Ent(normcdf(Eyt_la./sqrt(Varyt_la+1)))-efp;
    % Start monotonicity detection if we are close to the last trials
    if trial>experiments-5 && any(gp.code~=0),
        gpam=gp;gpam.xv=x(2:2:end,:);
        gpam.nvd =1;
        if ~isfield(gpam, 'lik_mono') || ~ismember(gpam.lik_mono.type, {'Probit', 'Logit'}) 
            gpam.lik_mono=lik_probit();
        end
        gpam.derivobs=1;
        gpam.lik_mono.nu=1e-6;
        gpam=gp_set(gpam,'latent_method','EP');

        nblocks = 10;
        [tmp,itst]=cvit(size(xopt,1),nblocks);
        gpam.nvd = [-1,-2];
        Ef=zeros(size(xopt,1),length(gpam.nvd));
        error = false;
        try
            for i=1:nblocks
                % Predict in blocks to save memory
                Ef(itst{i},:)=gpep_predgrad(gpam,x,y,xopt(itst{i},:));
            end
        catch ME
            error = true;
            gp.code = 0;
        end
        if ~error && any(any(Ef>0))
            gp.monotonicity = 1;
            gp.code = [gp.code 2];
            Entropy = Entropy .* -1.* max(Ef,[],2);
        else
            gp.code = [gp.code 0];
        end
    end

    [~,next_idx] = max(Entropy);%varyt_la

    next_design = xopt(next_idx,:);
end

% --- Convert time from number to text
function time = parse_time(logtime)
    global invTran
    days = invTran(logtime)*7;
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

function varargout = td_gui_OutputFcn(hObject, eventdata, handles) 
    varargout{1} = handles.output;
end

% --- Log user choice
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

% --- GPAL post processing
function setup_experiment(handles)
    global trial invTran yal xal va_t m xseed yseed b2 b3 gp
    x = [xseed; xal];
    y = [yseed; yal];
    gp.anchors = 1:length(yseed);

    h = waitbar(0,'Loading, please wait (this may take up to 20 seconds)');
    next_design = run_experiment(x,y,handles);
    close(h)
    xal = [xal; next_design];
    bt = va_t*sqrt(m^2+1)+2-m*2;
    bu = 2-m*2-va_t*sqrt(m^2+1);
    new_dist_top = abs(m*next_design(1)-next_design(2)+bt)/sqrt(m^2+1);
    new_dist_bot = abs(m*next_design(1)-next_design(2)+bu)/sqrt(m^2+1);

    change = new_dist_top>va_t || new_dist_bot>va_t||trial==1;

    if change
        coeff = polyfit(xal(:,1), xal(:,2), 1);
        if coeff(1)<0 && trial>10
            m=coeff(1);
        end
        bt = va_t*sqrt(m^2+1)+next_design(2)-m*next_design(1);
        bu = next_design(2)-m*next_design(1)-va_t*sqrt(m^2+1);
        margin = (0<m*xseed(:,1)+bt-xseed(:,2)) & (0>m*xseed(:,1)+bu-xseed(:,2));
        xseed(margin,:)=[];
        yseed(margin)=[];
    end
    
    trial=trial+1;
    time_string = parse_time(next_design(1));
    money_string = sprintf('$%0.0f',invTran(next_design(2)));
    rnum = randi(2);
    if rnum==1
        set(handles.pushbutton2,'string',[money_string, ' Now']);
        set(handles.pushbutton3,'string',['$800 in ', time_string]);
        b2=-1;
        b3=1;
    else
        set(handles.pushbutton3,'string',[money_string, ' Now']);
        set(handles.pushbutton2,'string',['$800 in ', time_string]);
        b2=1;
        b3=-1;
    end
    set(handles.text4,'string',['Trial: ', num2str(trial)]);
    set(handles.pushbutton2,'Enable','on');
    set(handles.pushbutton3,'Enable','on');
end    

% --- Choose action after choice
function choose_experiment(handles)
    global trial yal xal experiments rep gp xopt xseed yseed margin_bot margin_top staging;
    if trial<experiments
        setup_experiment(handles);
    else
        if rep==1
            save_results(xal,yal);
            rep=-1;
            clear gp xopt xal yal xseed yseed margin_bot margin_top;
            initialize_experiment(handles);
            setup_experiment(handles);
        else
            save_results(xal,yal);
            if staging
                addpath(genpath('ADO'));
                close(td_gui);
                run('td_ado_gui');

            else
                set(handles.pushbutton2,'visible','off');
                set(handles.pushbutton3,'visible','off');
                set(handles.text6,'visible','on');
            end
        end
    end
end    



% --- Disable some buttons while handeling data
function pushbutton3_Callback(hObject, eventdata, handles)
    global yal b3
    set(handles.pushbutton2,'Enable','off');
    set(handles.pushbutton3,'Enable','off');
    yal = [yal; b3];
    choose_experiment(handles);
end

% --- Save results into data folder
function save_results(xal,yal)
    global gp xseed yseed rep machine
    exists = true;
    id = 1;
    save_dir='data/';
    suffix ='';
    if rep==1
        suffix=[machine,'_1'];
    elseif rep==-1
        suffix = [machine,'_2'];
    end
    while exists
        name= [save_dir,num2str(id),suffix,'.mat'];
        if exist(name,'file')==2
            id=id+1;
        else
            save(name,'xal','yal','gp','xseed','yseed');
            exists=false;
        end
    end
end


% --- Show instructions
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
