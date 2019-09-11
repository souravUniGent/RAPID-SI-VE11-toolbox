function[]= preprocessing_new()

% This function is used for pre-processing the data: Phase correction, Apodisation

% Author: Sourav Bhaduri
% University of Ghent
% email: Sourav.Bhaduri@UGent.be
% Feb 2018; Last revision: Feb 2018


global data_value;
global preprocess_handles;
global details_para;
global apod_para;
global pha_para;
global ind_allvox;
global pha_cor_data;
global default_parameters;
global zerofill_value;
global info_struct;
global cur_vox;
global preprocess_display;
global FID_FD_copy;
global met_tick;
global mainGUI_handles;

FID_FD_copy = data_value.FID_FD;
default_parameters.pha_0=0;
default_parameters.pha_1=0;
default_parameters.apod_para_lor=0;
default_parameters.apod_para_gaus=0;
default_parameters.apod_para_sigm=0;
default_parameters.apod_para_type=1;
ind_allvox=1;
met_tick=0;
info_struct = [];
details_para.selected_voxels=1:size(data_value.FID_FD,2);
details_para.num_vox=size(data_value.FID_FD,2);
details_para.N=size(data_value.FID_FD,1);
cur_vox = details_para.selected_voxels(1,1);
details_para.disp_ppm_seg=[-1+2.01 4+2.01];
details_para.disp_freq_seg(1) = ((details_para.disp_ppm_seg(1) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.disp_seg(1) = ceil((details_para.disp_freq_seg(1)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.disp_freq_seg(2) = ((details_para.disp_ppm_seg(2) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.disp_seg(2) = ceil((details_para.disp_freq_seg(2)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
pha_para.temp_pha0=default_parameters.pha_0*ones(details_para.num_vox,1);
pha_para.temp_pha1=default_parameters.pha_1*ones(details_para.num_vox,1);
pha_para.pha_0=default_parameters.pha_0*ones(details_para.num_vox,1);
pha_para.pha_1=default_parameters.pha_1*ones(details_para.num_vox,1);
apod_para.type = default_parameters.apod_para_type*ones(details_para.num_vox,1);
apod_para.lor = default_parameters.apod_para_lor*ones(details_para.num_vox,1);
apod_para.gaus = default_parameters.apod_para_gaus*ones(details_para.num_vox,1);
apod_para.sigm = default_parameters.apod_para_sigm*ones(details_para.num_vox,1);
apod_para.temp_lor = default_parameters.apod_para_lor*ones(details_para.num_vox,1);
apod_para.temp_gaus = default_parameters.apod_para_gaus*ones(details_para.num_vox,1);
apod_para.temp_sigm = default_parameters.apod_para_sigm*ones(details_para.num_vox,1);
pha_cor_data =  data_value.FID_TD;
if(isempty(apod_para)) % if pre-processing is not done previously load default pre-processing parameters
    
    apod_para.type = default_parameters.apod_para_type*ones(details_para.num_vox,1);
    apod_para.lor = default_parameters.apod_para_lor*ones(details_para.num_vox,1);
    apod_para.gaus = default_parameters.apod_para_gaus*ones(details_para.num_vox,1);
    apod_para.sigm = default_parameters.apod_para_sigm*ones(details_para.num_vox,1);
    apod_para.function_weight = ones(size(data_value.FID_FD));
    pha_para.pha_0 = default_parameters.pha_0*ones(details_para.num_vox,1);
    pha_para.pha_1 = default_parameters.pha_1*ones(details_para.num_vox,1);
    pha_para.temp_pha0=default_parameters.pha_0*ones(details_para.num_vox,1);
    pha_para.temp_pha1=default_parameters.pha_1*ones(details_para.num_vox,1);
    pha_cor_data =  data_value.FID_TD;
else % if pre-processing is done use the current pre-processing parameters 
    % and perform pre-processing using those values 
    for i = 1:length(details_para.selected_voxels)
        vox = details_para.selected_voxels(1,i);
%         pha_cor_data(:,vox) = phase_correction(FID_FD_copy(:,vox),pha_para.pha_0(vox),pha_para.pha_1(vox),details_para.fres); %phase correction is done on frequency domain non-fftshifted data
%         [temp2, function_weight] = apodize(pha_cor_data(:,vox),apod_para.type(vox),[details_para.Fs,apod_para.lor(vox),apod_para.gaus(vox),apod_para.sigm(vox)]);
%         sht_sig = fftshift(fft(temp2),1);
%         apod_para.function_weight(:, vox) = function_weight;
%         data_value.FID_TD(:,vox) = temp2;
%         data_value.FID_FD(:,vox) = sht_sig;
        apod_para.type=1*ones(details_para.num_vox,1);
    end
end
x_max = 850;
y_max = 600;
selected_vox = details_para.selected_voxels;
vox_no = cell(1,length(selected_vox));
for i = 1:length(selected_vox)
vox_no{i} =  num2str(selected_vox(1,i));
end

% Cerate GUI
preprocess_handles.fig = figure('Visible','off','Position',[100,100,x_max,y_max],'Color',[0.68,0.92,1],'Name','RAPID 1.0','NumberTitle', 'off');
set(preprocess_handles.fig,'MenuBar','none'); 
% preprocess_handles.save_figure = uimenu(preprocess_handles.fig,'Label','Save Figure','Callback',@save_preprocess_fig_callback);

preprocess_handles.h1 = uipanel(preprocess_handles.fig,'Title','Voxel No','Units','normalized','Position',[15/x_max 190/y_max 70/x_max 270/y_max],...
    'BackgroundColor','white','FontSize',8, 'FontWeight','Bold');
% preprocess_handles.voxel_list = uicontrol(preprocess_handles.h1,'Style', 'listbox','Units','normalized','Position', [.05 .027 .9 .97],'String', vox_no,...
%     'Max', details_para.num_vox, 'Min', 1, 'Callback', @voxel_list_callback, 'Interruptible','off','BusyAction','cancel'); 
preprocess_handles.voxel_list = uicontrol(preprocess_handles.h1,'Style', 'listbox','Units','normalized','Position', [.05 .027 .9 .97],'String', vox_no,...
    'value',1, 'Callback', @voxel_list_callback, 'Interruptible','off','BusyAction','cancel'); 

preprocess_handles.axis = axes('Units','Pixels','Units','normalized','Position',[120/x_max,200/y_max,420/x_max,250/y_max]);
% preprocess_handles.axis_text = uicontrol(preprocess_handles.fig,'Style','text','Units','normalized','Position',[220/x_max 460/y_max 130/x_max 15/y_max],...
%     'String','Signal Spectrum','BackgroundColor',[0.68,0.92,1],'FontSize',9,'FontWeight','Bold');

preprocess_handles.TF_FD_display = uibuttongroup(preprocess_handles.fig,'Units','normalized','Position',[0.1 0.92 .4 0.06],...
    'FontSize',9,'BackgroundColor',[0.68,0.92,1],'BorderType','none','SelectionChangeFcn',@TD_FD_display_callback);
preprocess_handles.FD_disp_real = uicontrol(preprocess_handles.TF_FD_display,'Style','radiobutton','Units','normalized','Position',[.01 .1 .3 .8],'String',...
    'FD signal real','FontWeight','Bold','BackgroundColor',[0.68,0.92,1]);
preprocess_handles.FD_disp_abs = uicontrol(preprocess_handles.TF_FD_display,'Style','radiobutton','Units','normalized','Position',[.35 .1 .3 .8],'String',...
    'FD signal abs','FontWeight','Bold','BackgroundColor',[0.68,0.92,1]);
preprocess_handles.TD_disp = uicontrol(preprocess_handles.TF_FD_display,'Style','radiobutton','Units','normalized','String','TD signal','Position',...
    [.71 .1 .3 .8],'FontWeight','Bold','BackgroundColor',[0.68,0.92,1]);


preprocess_handles.display_panel = uipanel(preprocess_handles.fig,'Title','Display Properties','Position',[.1 .8 .27 .1],'BackgroundColor','white','FontSize',9,'FontWeight','Bold');
preprocess_handles.seg_min_text = uicontrol(preprocess_handles.display_panel,'Style','text','Units','normalized','Position',[0.4,0.7,0.2,0.1*2.5],'String',...
    'Min','FontSize',8,'FontWeight','Bold','BackgroundColor','white');
preprocess_handles.seg_max_text = uicontrol(preprocess_handles.display_panel,'Style','text','Units','normalized','Position',[0.7,0.7,0.2,0.1*2.5],'String',...
    'Max','FontSize',8,'FontWeight','Bold','BackgroundColor','white');
preprocess_handles.display_seg_mim = uicontrol(preprocess_handles.display_panel,'Style','edit','Units','normalized','Position',[0.45,0.2,0.15,0.15*2],'String','0',...
    'Callback',{@display_seg_mim_Callback},'Interruptible','off','BusyAction','cancel'); 
preprocess_handles.display_seg_max = uicontrol(preprocess_handles.display_panel,'Style','edit','Units','normalized','Position',[0.75,0.2,0.15,0.15*2],'String','0',...
    'Callback',{@display_seg_max_Callback},'Interruptible','off','BusyAction','cancel'); 
% preprocess_handles.zerofill_text = uicontrol(preprocess_handles.fig,'Style','text','Units','normalized','Position',[557/x_max 460/y_max 70/x_max 15/y_max],...
%     'String','Zerofill','BackgroundColor',[0.68,0.92,1],'FontSize',9,'FontWeight','Bold');
% preprocess_handles.zerofill_options = uicontrol(preprocess_handles.fig,'Style','popup','String','No Fill|1 times|2 times|3 times','Units','normalized',...
%     'Position',[555/x_max 405/y_max 80/x_max 15/y_max],'Value',zerofill_value,'Interruptible','off','BusyAction','cancel','Callback', @zerofill_options_callback); 

preprocess_handles.file_panel = uipanel(preprocess_handles.fig,'Title','','Position',[.5 .9 .43 .1],'BackgroundColor','white','FontSize',5,'FontWeight','Bold');
preprocess_handles.main_text = uicontrol(preprocess_handles.file_panel,'Style','text','Units','normalized','Position',[0.005,0.65,0.05,0.1*2.5],'String',...
    'M','FontSize',6,'FontWeight','Bold','BackgroundColor','white');
preprocess_handles.ref_text = uicontrol(preprocess_handles.file_panel,'Style','text','Units','normalized','Position',[0.005,0.1,0.05,0.1*2.5],'String',...
    'R','FontSize',6,'FontWeight','Bold','BackgroundColor','white');
preprocess_handles.main_val = uicontrol(preprocess_handles.file_panel,'Style','edit','Units','normalized','Position',[0.05,0.65,0.95,0.15*2],'String',data_value.file_name{1},...
    'Callback',{@file_main_Callback},'Interruptible','off','BusyAction','cancel'); 
preprocess_handles.ref_val = uicontrol(preprocess_handles.file_panel,'Style','edit','Units','normalized','Position',[0.05,0.1,0.95,0.15*2],'String',data_value.file_name{2},...
    'Callback',{@file_ref_Callback},'Interruptible','off','BusyAction','cancel'); 
% apodization panel
preprocess_handles.h2 = uipanel(preprocess_handles.fig,'Title','Apodization','Units','normalized','Position',[590/x_max 40/y_max 150/x_max 350/y_max],...
    'BackgroundColor','white','FontSize',9, 'FontWeight','Bold');
preprocess_handles.Apodize = uicontrol(preprocess_handles.h2,'Style','text','Units','normalized','String','Apodization Type','Units','normalized',...
    'Position',[0.1 0.65 0.8 0.3],'BackgroundColor','white','FontSize',8,'FontWeight','Bold');
preprocess_handles.Apodize_type = uicontrol(preprocess_handles.h2,'Style','popup','String','Exponential|Gauss|Gauss-Exponential|sigmoid','Units','normalized','Position',...
    [0.1 0.59 0.8 0.3],'Value',apod_para.type(cur_vox), 'Callback',{@Apodize_type_Callback});

preprocess_handles.Apodize_Lor_text = uicontrol(preprocess_handles.h2, 'Style','text','Units','normalized','Position',[0.1 0.7 0.8 0.07],'String',...
    'Exponential','FontSize',8,'FontWeight','Bold','BackgroundColor','white');
preprocess_handles.Apodize_Lor_edit = uicontrol(preprocess_handles.h2,'Style','edit','Units','normalized','Position',[0.1 0.65 0.8 0.07],'String',...
    num2str(apod_para.lor(cur_vox)),'Callback',{@lor_Callback}); 
preprocess_handles.Apodize_Lor_slide = uicontrol(preprocess_handles.h2,'Style', 'slider','Min',0,'Max',100,'Value',apod_para.lor(cur_vox),'Units','normalized',...
    'Position', [0.1 0.58 0.8 0.07],'SliderStep',[0.001,0.01],'Callback', {@slider_lor_Callback}); 
%     
preprocess_handles.Apodize_Gaus_text = uicontrol(preprocess_handles.h2,'Style','text','Units','normalized','Position',[0.1 0.45 0.8 0.07],'String',....
    'Gauss','FontSize',8,'FontWeight','Bold','BackgroundColor','white');
preprocess_handles.Apodize_Gaus_edit = uicontrol(preprocess_handles.h2,'Style','edit','Units','normalized','Position',[0.1 0.4 0.8 0.07],...
    'String',num2str(apod_para.gaus(cur_vox)),'Callback',{@gaus_Callback});  
preprocess_handles.Apodize_Gaus_slide  = uicontrol(preprocess_handles.h2,'Style', 'slider','Min',0,'Max',1000,'Value',apod_para.gaus(cur_vox),'Units','normalized',...
    'Position', [0.1 0.33 0.8 0.07],'SliderStep',[0.0001,0.001],'Callback', {@slider_gaus_Callback});  

preprocess_handles.Apodize_Sig_text = uicontrol(preprocess_handles.h2,'Style','text','Units','normalized','Position',[0.1 0.2 0.8 0.07],'String',...
    'Sigmoid','FontSize',8,'FontWeight','Bold','BackgroundColor','white');
preprocess_handles.Apodize_Sigm_edit = uicontrol(preprocess_handles.h2,'Style','edit','Units','normalized','Position',[0.1 0.15 0.8 0.07],...
    'String',num2str(apod_para.sigm(cur_vox)),'Callback',{@sigm_Callback});  
preprocess_handles.Apodize_Sigm_slide = uicontrol(preprocess_handles.h2,'Style', 'slider','Min',0,'Max',500,'Value',apod_para.sigm(cur_vox),'Units','normalized',...
    'Position', [0.1 0.08 0.8 0.07],'SliderStep',[0.001,0.01],'Callback', {@slider_sigm_Callback});


% phase correction panel
preprocess_handles.h3 = uipanel(preprocess_handles.fig,'Title','Phase correction','Units','normalized','Position',[20/x_max 30/y_max 520/x_max 150/y_max],...
    'BackgroundColor','white','FontSize',9, 'FontWeight','Bold');
preprocess_handles.pha_corr0_text = uicontrol(preprocess_handles.h3,'Style','text','Units','normalized','Position',[0.1 0.5 0.2 .2],'String',...
    'Zero-order phase correction ','FontSize',8,'FontWeight','Bold','BackgroundColor','white');
preprocess_handles.pha_corr1_text = uicontrol(preprocess_handles.h3,'Style','text','Units','normalized','Position',[0.1 0.2 0.2 0.2],'String',...
    'First-order phase correction','FontSize',8,'FontWeight','Bold','BackgroundColor','white');
preprocess_handles.pha_corr_value_text = uicontrol(preprocess_handles.h3,'Style','text','Units','normalized','Position',[0.31 0.61 0.2 .2],'String',...
    'Value','FontSize',8,'FontWeight','Bold','BackgroundColor','white');
preprocess_handles.pha_corr_slider_text = uicontrol(preprocess_handles.h3,'Style','text','Units','normalized','Position',[0.66 0.61 0.2 0.2],'String',...
    'Slider','FontSize',8,'FontWeight','Bold','BackgroundColor','white');
preprocess_handles.pha_corr0_value = uicontrol(preprocess_handles.h3,'Style','edit','Units','normalized','Position',[0.34 0.5 .15 .2],...
    'String',num2str(pha_para.pha_0(cur_vox)),'Callback',{@pha0_val_Callback});
preprocess_handles.pha_corr1_value = uicontrol(preprocess_handles.h3,'Style','edit','Units','normalized','Position',[0.34 0.2 .15 .2],...
    'String',num2str(pha_para.pha_1(cur_vox)),'Callback',{@pha1_val_Callback});
preprocess_handles.pha_corr0_slide = uicontrol(preprocess_handles.h3,'Style','slider','Min',-pi,'Max',pi,'Value',pha_para.pha_0(cur_vox),'Units','normalized','Position',...
    [0.55 0.5 .42 .2],'SliderStep',[0.001,0.01],'Callback', {@pha0_slider_Callback});  
preprocess_handles.pha_corr1_slide = uicontrol(preprocess_handles.h3,'Style','slider','Min',-10,'Max',10,'Value',pha_para.pha_1(cur_vox),'Units','normalized','Position',...
    [0.55 0.2 .42 .2],'SliderStep',[0.001,0.01],'Callback', {@pha1_slider_Callback});

preprocess_handles.close = uicontrol(preprocess_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[680/x_max 10/y_max 100/x_max 25/y_max],...
'String', 'OK/Close','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@close_Callback);
% preprocess_handles.met_map = uicontrol(preprocess_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[560/x_max 10/y_max 100/x_max 25/y_max],...
% 'String', 'metabolite map','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@metabolite_map_Callback);

preprocess_handles.buttongroup = uibuttongroup(preprocess_handles.fig,'Units','normalized','Position',[655/x_max 410/y_max 130/x_max 70/y_max],...
    'BackgroundColor','white','FontSize',9, 'FontWeight','Bold');%,'SelectionChangeFcn',@ind_allvox_value);

preprocess_handles.all_vox = uicontrol('parent',preprocess_handles.buttongroup,'Style','radiobutton','String','All Voxel',...
    'Units','normalized','Position',[0.1 0.2 0.7 0.3],'BackgroundColor','white');
preprocess_handles.ind_vox = uicontrol('parent',preprocess_handles.buttongroup,'Style','radiobutton','String','Individual Voxel',...
    'Units','normalized','Position',[0.1 0.6 0.8 0.3],'BackgroundColor','white');
preprocess_handles.water_suppress = uicontrol(preprocess_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[240/x_max 150/y_max 100/x_max 25/y_max],...
'String', 'Water suppress','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@water_suppress_Callback);

% preprocess_handles.wavelet = uicontrol(preprocess_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[240/x_max 140/y_max 100/x_max 25/y_max],...
% 'String', 'Denoise','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@denoise_Callback);
set(preprocess_handles.buttongroup,'SelectionChangeFcn',@ind_allvox_value);



if ind_allvox == 1
    set(preprocess_handles.buttongroup,'SelectedObject',preprocess_handles.ind_vox);
else
    set(preprocess_handles.buttongroup,'SelectedObject',preprocess_handles.all_vox);
end

if (apod_para.type(selected_vox(1,1)) == 1)
    set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'on');
    set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'on');
    set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'off');
    set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'off');
    set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'off');
    set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'off');
elseif (apod_para.type(selected_vox(1,1)) == 2)
    set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'off');
    set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'off');
    set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'off');
    set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'off');
    set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'on');
    set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'on');
elseif (apod_para.type(selected_vox(1,1)) == 3)
    set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'on');
    set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'off');
    set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'on');
    set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'off');
    set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'on');
    set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'on');
else
    set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'off');
    set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'on');
    set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'off');
    set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'on');
    set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'off');
    set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'off');
end

    

set(preprocess_handles.fig,'Visible','on');
set(preprocess_handles.display_seg_mim, 'Enable', 'on','string',num2str(details_para.disp_ppm_seg(1)))
set(preprocess_handles.display_seg_max, 'Enable', 'on','string',num2str(details_para.disp_ppm_seg(2)))
set(preprocess_handles.main_val, 'Enable', 'off')
set(preprocess_handles.ref_val, 'Enable', 'off')

 
preprocess_display = 1;

% Display the selected voxel 
display_data()
waitfor(preprocess_handles.fig);

function voxel_list_callback(hObject, eventdata)
% Callback for 'voxel no.' get the selected voxel and update the plot according to the voxel selected
global details_para;
global preprocess_handles;
global apod_para;
global pha_para;
global cur_vox;

% get the selected voxel
preprocess_sel_vox = get(hObject,'Value');

if(isempty(preprocess_sel_vox))
    preprocess_sel_vox = 1;
    set(preprocess_handles.voxel_list,'value',1);
end

selected_vox = details_para.selected_voxels;
cur_vox = selected_vox(1,preprocess_sel_vox);
display_data()


% Update the GUI based on the selected voxel
temp = apod_para.type(cur_vox);
set(preprocess_handles.Apodize_type,'value',temp);
switch(temp)
  case 1
        set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'on','string',num2str(apod_para.lor(cur_vox)));
        set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'on','value',apod_para.lor(cur_vox));
        set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'off');
        set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'off');
        
  case 2
        set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'on','string',num2str(apod_para.gaus(cur_vox)));
        set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'off');
        set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'on','value',apod_para.gaus(cur_vox));
        set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'off');
  case 3
        set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'on','string',num2str(apod_para.lor(cur_vox)));
        set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'on','string',num2str(apod_para.gaus(cur_vox)));
        set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'on','value',apod_para.lor(cur_vox));
        set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'on','value',apod_para.gaus(cur_vox));
        set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'off');
  case 4
        set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'on','string',num2str(apod_para.sigm(cur_vox)));
        set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'off');
        set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'off');
        set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'on','value',apod_para.sigm(cur_vox));
end

set(preprocess_handles.pha_corr0_slide,'value',pha_para.pha_0(cur_vox));
set(preprocess_handles.pha_corr1_slide,'value',pha_para.pha_1(cur_vox));
set(preprocess_handles.pha_corr0_value,'string',num2str(pha_para.pha_0(cur_vox)));
set(preprocess_handles.pha_corr1_value,'string',num2str(pha_para.pha_1(cur_vox)));
% set(preprocess_handles.display_seg_mim, 'Enable', 'on','string',num2str(details_para.disp_ppm_seg(1)))
% set(preprocess_handles.display_seg_max, 'Enable', 'on','string',num2str(details_para.disp_ppm_seg(2)))


function zerofill_options_callback(hObject, eventdata)
%Callback for 'Zero-fill' drop down box
global zerofill_value;
global details_para;

% update the selected zero-fill value, But dont zero fill the signal.
% Zero-filling is done once the pre-processing GUI is closed
zerofill_value = get(hObject,'Value');
if(zerofill_value ~= 1)
    [a,b] = size(data_value.FID_TD);
    data_value.FID_TD = [data_value.FID_TD;zeros((zerofill_value-1)*a,b)]; %add zeros in the end
    data_value.FID_FD = fftshift(fft(data_value.FID_TD),1);
    
    %update parameters after zerofilling
    N = (zerofill_value)*(a); 
    n = 0:N-1;
    fres = (-details_para.Fs/2 + (details_para.Fs/N)*n);
    ref = details_para.ref;
    ppm = (fres)*(1E6/details_para.Tf);
    t = n/details_para.Fs;
    details_para.N = N;
    details_para.ppm = ppm;
    details_para.fres = fres;
    details_para.t = t;
    details_para.ppm_referenced = details_para.ref + ppm;
end
display_data();



function Apodize_type_Callback(hObject, eventdata)
%Callback for 'Apodisation type' 
global preprocess_handles;
global apod_para;
global ind_allvox;
global details_para;
global cur_vox;
% global data_value;

% get the selected apodization type
temp = get(hObject,'Value');
if(ind_allvox == 1)
    apod_para.type(cur_vox) = temp;
else
    apod_para.type(details_para.selected_voxels(1,:)) = temp;
end

% enable/disable based on the selected type
switch(temp)
  case 1 %Lorentz
        set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'on');
        set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'on');
        set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'off');
        set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'off');
  case 2 %Gaussian
        set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'on');
        set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'off');
        set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'on');
        set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'off');
  case 3 %Lor-Gauss
        set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'on');
        set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'on');
        set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'on');
        set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'on');
        set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'off');
  case 4 %sigmoid
        set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'on');
        set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'off');
        set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'off');
        set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'on');
end
do_and_display_apod()% apodize and display the pre-processed signal
  
function lor_Callback(hObject, eventdata)
%Callback for 'Lorentz' edit box  
global preprocess_handles;
global apod_para;
global ind_allvox;
global details_para;
global cur_vox;

% read and check the Lorentz apod para from the edit box 
temp = str2double(get(hObject,'string'));
if((temp < 0) || (temp > 100) || (isnan(temp)))
    errordlg('Invalid value');
    set(hObject,'string',num2str(apod_para.lor(cur_vox)))
    return;
end
if(ind_allvox == 1)
    apod_para.lor(cur_vox) = temp;
else
    apod_para.lor(details_para.selected_voxels(1,:)) = temp;
end
set(preprocess_handles.Apodize_Lor_slide,'value',apod_para.lor(cur_vox));

do_and_display_apod(); % apodize and display the pre-processed signal

function slider_lor_Callback(hObject, eventdata)
%Callback for 'Lorentz' slider  
global preprocess_handles;
global apod_para;
global ind_allvox;
global details_para;
global cur_vox;

% read and ckeck the Lorentz apod para from the Lorentz slider 
temp = get(hObject,'Value');
if(ind_allvox == 1)
    apod_para.temp_lor(cur_vox)=temp - apod_para.lor(cur_vox);
    apod_para.lor(cur_vox) = temp;
else
    apod_para.temp_lor(details_para.selected_voxels(1,:))=temp- apod_para.lor(details_para.selected_voxels(1,:));
    apod_para.lor(details_para.selected_voxels(1,:)) = temp;  
end
set(preprocess_handles.Apodize_Lor_edit,'string',num2str(apod_para.lor(cur_vox)));
do_and_display_apod(); % apodize and display the pre-processed signal

function gaus_Callback(hObject, eventdata)
%Callback for 'Gaussian' edit box 
global preprocess_handles;
global apod_para;
global cur_vox;
global ind_allvox;
global details_para;

% read and check the Gaussian apod para from the edit box 
temp = str2double(get(hObject,'string'));
if((temp < 0) || (temp > 1000) || (isnan(temp)))
    errordlg('Invalid value');
    set(hObject,'string',num2str(apod_para.gaus(cur_vox)))
    return;
end
if(ind_allvox == 1)
    apod_para.gaus(cur_vox) = temp;
else
    apod_para.gaus(details_para.selected_voxels(1,:)) = temp;
end
set(preprocess_handles.Apodize_Gaus_slide,'value',apod_para.gaus(cur_vox));
do_and_display_apod();% apodize and display the pre-processed signal


function slider_gaus_Callback(hObject, eventdata)
%Callback for 'Gaussian' slider 
global preprocess_handles;
global apod_para;
global cur_vox;
global ind_allvox;
global details_para;

% read and check the Gaussian apod para from the slider 
temp = get(hObject,'Value');
if(ind_allvox == 1)
    apod_para.temp_gaus(cur_vox)=temp - apod_para.gaus(cur_vox);
    apod_para.gaus(cur_vox) = temp;
else
    apod_para.temp_gaus(details_para.selected_voxels(1,:))=temp- apod_para.gaus(details_para.selected_voxels(1,:));
    apod_para.gaus(details_para.selected_voxels(1,:)) = temp;  
end

set(preprocess_handles.Apodize_Gaus_edit,'string',num2str(apod_para.gaus(cur_vox)));
do_and_display_apod();% apodize and display the pre-processed signal

% 
function sigm_Callback(hObject, eventdata)
%Callback for 'Sigmoid' edit box 
global preprocess_handles;
global apod_para;
global cur_vox;
global ind_allvox;
global details_para;

% read and check the sigmoid apod para from the edit box 
temp = str2double(get(hObject,'string'));
if((temp < 0)|| (temp > 500) || (isnan(temp)))
    errordlg('Invalid value');
    set(hObject,'string',num2str(apod_para.sigm(cur_vox)))
    return;
end
if(ind_allvox == 1)
    apod_para.sigm(cur_vox) = temp;
else
    apod_para.sigm(details_para.selected_voxels(1,:)) = temp;
end
set(preprocess_handles.Apodize_Sigm_slide,'value',apod_para.sigm(cur_vox));
do_and_display_apod();% apodize and display the pre-processed signal


function slider_sigm_Callback(hObject, eventdata)
%Callback for 'Sigmoid' slider 
global preprocess_handles;
global apod_para;
global cur_vox;
global ind_allvox;
global details_para;

% read and check the sigmoid apod para from the slider
temp = get(hObject,'Value');
if(ind_allvox == 1)
    apod_para.temp_sigm(cur_vox)=temp - apod_para.sigm(cur_vox);
    apod_para.sigm(cur_vox) = temp;
else
    apod_para.temp_sigm(details_para.selected_voxels(1,:))=temp- apod_para.sigm(details_para.selected_voxels(1,:));
    apod_para.sigm(details_para.selected_voxels(1,:)) = temp;  
end
set(preprocess_handles.Apodize_Sigm_edit,'string',num2str(apod_para.sigm(cur_vox)));
do_and_display_apod();% apodize and display the pre-processed signal

% 
function pha0_val_Callback(hObject, eventdata)
%Callback for 'zero-order' edit box 
global pha_para;
global cur_vox;
global ind_allvox;
global preprocess_handles;
global details_para;

% read and check the zero-order phase value from the edit box 
phi_0 = str2double(get(hObject,'string'));
if((phi_0 < -pi)|| (phi_0 > pi) || (isnan(phi_0)))
    errordlg('Invalid value');
    set(hObject,'string',num2str(pha_para.pha_0(cur_vox)))
    return;
end

if(ind_allvox == 1)
    pha_para.pha_0(cur_vox) = phi_0;
else
    pha_para.pha_0(details_para.selected_voxels(1,:)) = phi_0;    
end
set(preprocess_handles.pha_corr0_slide,'value',pha_para.pha_0(cur_vox));
do_and_display_phase() % Phase-correct and display the pre-processed signal



% 
function pha1_val_Callback(hObject, eventdata)
%Callback for 'first-order' edit box 
global pha_para;
global cur_vox;
global ind_allvox;
global preprocess_handles;
global details_para;

% read and check the first-order phase value from the edit box 
phi_1 = str2double(get(hObject,'string'));
if((phi_1 < -10)|| (phi_1 > 10) || (isnan(phi_1)))
    errordlg('Invalid value');
    set(hObject,'string',num2str(pha_para.pha_1(cur_vox)))
    return;
end
if(ind_allvox == 1)
    pha_para.pha_1(cur_vox) = phi_1;
else
    pha_para.pha_1(details_para.selected_voxels(1,:)) = phi_1;    
end
set(preprocess_handles.pha_corr1_slide,'value',pha_para.pha_1(cur_vox));
do_and_display_phase() % Phase-correct and display the pre-processed signal

function file_main_Callback(hObject, eventdata)
%Callback for 'first-order' edit box 
global pha_para;
global cur_vox;
global ind_allvox;
global preprocess_handles;
global details_para;
 

function file_ref_Callback(hObject, eventdata)
%Callback for 'first-order' edit box 
global pha_para;
global cur_vox;
global ind_allvox;
global preprocess_handles;
global details_para;
% 

function pha0_slider_Callback(hObject, eventdata)
%Callback for 'zero-order' slider 
global preprocess_handles;
global pha_para;
global cur_vox;
global ind_allvox;
global details_para;

% read and check the zero-order phase value from the slider 
phi_0 = get(hObject,'Value');
if(ind_allvox == 1)
    pha_para.temp_pha0(cur_vox)=phi_0- pha_para.pha_0(cur_vox);
    pha_para.pha_0(cur_vox) = phi_0;
else
    pha_para.temp_pha0(details_para.selected_voxels(1,:))=phi_0- pha_para.pha_0(details_para.selected_voxels(1,:));
    pha_para.pha_0(details_para.selected_voxels(1,:)) = phi_0;  
end
set(preprocess_handles.pha_corr0_value,'string',num2str(pha_para.pha_0(cur_vox)));
do_and_display_phase();% Phase-correct and display the pre-processed signal



function pha1_slider_Callback(hObject, eventdata)
%Callback for 'first-order' slider 
global pha_para;
global cur_vox;
global ind_allvox;
global preprocess_handles;
global details_para;

% read and check the first-order phase value from the slider 
phi_1 = get(hObject,'Value');
if(ind_allvox == 1)
    pha_para.temp_pha1(cur_vox)=phi_1- pha_para.pha_1(cur_vox);
    pha_para.pha_1(cur_vox) = phi_1;
else
    pha_para.temp_pha1(details_para.selected_voxels(1,:))=phi_1- pha_para.pha_1(details_para.selected_voxels(1,:));
    pha_para.pha_1(details_para.selected_voxels(1,:)) = phi_1; 
end
set(preprocess_handles.pha_corr1_value,'string',num2str(pha_para.pha_1(cur_vox)));
do_and_display_phase()% Phase-correct and display the pre-processed signal




function ind_allvox_value(source,eventdata)
% call back for individual voxel and all voxel button group
global ind_allvox;

% select whether to apply pre-processing on all the voxels or in the selected voxel 
temp = get(eventdata.NewValue,'string');
if(strcmp(temp,'Individual Voxel'))
    ind_allvox = 1;
else
    ind_allvox = 2;
end

function do_and_display_phase()
% This function performs phase correction and displays the signal
global ind_allvox;
global cur_vox;
global pha_para;
global pha_cor_data;
global data_value;
global apod_para;
global details_para;
global preprocess_handles;
global FID_FD_copy;

if(ind_allvox == 1) % perform phase correction on the selected vovel
    %     do phase correction
    [td_sig] = phase_correction(data_value.FID_FD(:,cur_vox),pha_para.temp_pha0(cur_vox),pha_para.temp_pha1(cur_vox),details_para.fres);
    pha_cor_data(:,cur_vox) = td_sig;
%     pha_para.pha_0(cur_vox)=order;
%     do apodization
    [temp2, function_weight] = apodize(td_sig,apod_para.type(cur_vox),[details_para.Fs,apod_para.lor(cur_vox),apod_para.gaus(cur_vox),apod_para.sigm(cur_vox)]);
    sht_sig = fftshift(fft(temp2));
    data_value.FID_TD(:,cur_vox) =(temp2);
    data_value.FID_FD(:,cur_vox) = (sht_sig);
    apod_para.function_weight(:, cur_vox) = function_weight;
else % perform phase correction on the all vovel
    %     do phase correction
%     td_sig = phase_correction(data_value.FID_FD_org,pha_para.pha_0(preprocess_sel_vox),pha_para.pha_1(preprocess_sel_vox));
    for i=1:length(details_para.selected_voxels(1,:))
        td_sig = phase_correction(data_value.FID_FD(:,i),pha_para.temp_pha0(cur_vox),pha_para.temp_pha1(cur_vox),details_para.fres);
        pha_cor_data(:,i) = td_sig;
    %     do apodization
%     temp2=td_sig;
        [temp2, function_weight]= apodize(td_sig,apod_para.type(cur_vox),[details_para.Fs,apod_para.lor(cur_vox),apod_para.gaus(cur_vox),apod_para.sigm(cur_vox)]);
        sht_sig = fftshift(fft(temp2));
%     apod_para.function_weight(:, details_para.selected_voxels(1,:)) = function_weight;
        data_value.FID_TD(:,i) = temp2;
        data_value.FID_FD(:,i) = sht_sig;  
    end
end
display_data()
set(preprocess_handles.pha_corr0_value,'string',num2str(pha_para.pha_0(cur_vox)));
set(preprocess_handles.Apodize_Sigm_edit,'string',num2str(apod_para.sigm(cur_vox)));

function do_and_display_apod()
% This function performs apodization and displays the signal
global details_para;
global preprocess_handles;
global data_value;
global cur_vox;
global apod_para;
global pha_cor_data;
global ind_allvox;


if(ind_allvox == 1) % perform apodization on the selected vovel
    temp1 = pha_cor_data(:,cur_vox);
        %     do apodization
    [temp2, function_weight] = apodize(temp1,apod_para.type(cur_vox),[details_para.Fs,apod_para.lor(cur_vox),apod_para.gaus(cur_vox),apod_para.sigm(cur_vox)]);
    sht_sig = fftshift(fft(temp2),1);
    data_value.FID_TD(:,cur_vox) = temp2;
    data_value.FID_FD(:,cur_vox) = sht_sig;
    apod_para.function_weight(:, cur_vox) = function_weight;
else %perform apodization on the all vovel
    temp1 = pha_cor_data(:,details_para.selected_voxels(1,:));
        %     do apodization
    [temp2, function_weight] = apodize(temp1,apod_para.type(cur_vox),[details_para.Fs,apod_para.lor(cur_vox),apod_para.gaus(cur_vox),apod_para.sigm(cur_vox)]);
    sht_sig = fftshift(fft(temp2),1);
    apod_para.function_weight(:, details_para.selected_voxels(1,:)) = function_weight;
    data_value.FID_TD(:,details_para.selected_voxels(1,:)) = temp2;
    data_value.FID_FD(:,details_para.selected_voxels(1,:)) = sht_sig;        
end

display_data()
set(preprocess_handles.Apodize_Sigm_edit,'string',num2str(apod_para.sigm(cur_vox)));


function TD_FD_display_callback(hObject, eventdata)
global preprocess_display;

temp = get(get(hObject,'SelectedObject'),'String');

if strcmp(temp,'FD signal real')
    preprocess_display = 1;
elseif strcmp(temp,'FD signal abs')
    preprocess_display = 2;
else
    preprocess_display = 0;
end
display_data()

function display_data()
global preprocess_display;
global preprocess_handles;
global cur_vox;
global details_para;
global data_value;
global apod_para;

seg = details_para.disp_seg;
ppm = details_para.ppm_referenced;

if (preprocess_display == 1)
    cla(preprocess_handles.axis);
    plot_var = real(data_value.FID_FD(seg(1):seg(2),cur_vox));
    plot(preprocess_handles.axis,ppm(seg(1):seg(2))',plot_var);
    max_p = max(plot_var);
    min_p = min(plot_var);
    axis(preprocess_handles.axis,[ppm(seg(1)), ppm(seg(2)),min_p - 0.1*max_p, max_p + 0.1*max_p])
    set(preprocess_handles.axis,'XDir','reverse');
elseif (preprocess_display == 2)
    cla(preprocess_handles.axis);
    plot_var = abs(data_value.FID_FD(seg(1):seg(2),cur_vox));
    plot(preprocess_handles.axis,ppm(seg(1):seg(2))',plot_var);
    max_p = max(plot_var);
    min_p = min(plot_var);
    axis(preprocess_handles.axis,[ppm(seg(1)), ppm(seg(2)),min_p - 0.1*max_p, max_p + 0.1*max_p])
    set(preprocess_handles.axis,'XDir','reverse');
else
    cla(preprocess_handles.axis);
    plot_var = real(data_value.FID_TD(:,cur_vox));
%     plot_var= ifftshift(real(data_value.FID_TD(:,cur_vox))+flipud(real(data_value.FID_TD(:,cur_vox))));
%     filter = apod_para.function_weight(:,cur_vox);
    max_p = max(abs(plot_var));
    min_p = max(abs(plot_var));
%     if (max(filter)> max_p)
%         max_p = max(filter);
%         min_p = -max_p;
%     else
%         max_p = max(plot_var);
%         min_p = -max_p;
%     end    
    plot(preprocess_handles.axis,details_para.t',plot_var);
%     hold on;
%     plot(preprocess_handles.axis,details_para.t',filter,'r');
%     axis(preprocess_handles.axis,[details_para.t(1), details_para.t(details_para.N),min_p , max_p])
%     hold off
end


function save_preprocess_fig_callback(hObject, eventdata)
%Callback for 'save figure'
global preprocess_handles;
global log_file;
ini_path = log_file.saving_path;
[temp,s_path] = uiputfile(ini_path);
if(temp == 0)
    return;
else
    dot_pos = strfind(temp, '.');
    s_fln = temp(1:dot_pos-1);
    file_name = strcat(s_path,s_fln);
    saveas(preprocess_handles.fig,file_name,'png');
    log_file.saving_path = s_path;
end



function close_Callback(hObject, eventdata)
%Callback for 'close' button
global preprocess_handles;
global met_tick;
global mainGUI_handles;
if (met_tick==1 && ~isempty(preprocess_handles.met_fig))
    close(preprocess_handles.fig)
    close(preprocess_handles.met_fig);
else
     close(preprocess_handles.fig);
end
%-----------------------------------------------------------------------------------------------------------------------------------
function water_suppress_Callback(hObject, eventdata)

global data_value;
global details_para;
global preprocess_handles;
global preprocess_sel_vox;
global FID_FD_copy;
global pha_cor_data;

h = waitbar(0,'Please wait...');
for i = 1:details_para.num_vox
    [ fit_sig,amp,Ind_sig,freq] = HSVDPK_new(data_value.FID_TD(:,i),details_para.fres,details_para.Fs,10);
    ind=find((freq>=-80 & freq<=80));
    half_sig2=sum(Ind_sig(:,ind),2);
    data_value.FID_TD(:,i)=data_value.FID_TD(:,i)- half_sig2;
    data_value.FID_FD(:,i)=fftshift(fft(data_value.FID_TD(:,i)));
    waitbar(i/details_para.num_vox);
end
close(h);
% FID_FD_copy=data_value.FID_FD;
pha_cor_data(:,details_para.selected_voxels(1,:))=data_value.FID_TD;
% seg = details_para.disp_seg;
% ppm = details_para.ppm_referenced;
% cla(preprocess_handles.axis);
% plot_var = real(data_value.FID_FD(seg(1):seg(2)));
% plot(preprocess_handles.axis,ppm(seg(1):seg(2))',plot_var);
% max_p = max(plot_var);
% mim_p = min(plot_var);
% axis(preprocess_handles.axis,[ppm(seg(1)), ppm(seg(2)),mim_p - 0.1*max_p, max_p + 0.1*max_p])
% set(preprocess_handles.axis,'XDir','reverse')
display_data()
%------------------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------------------
function metabolite_map_Callback(hObject, eventdata)

global data_value;
global details_para;
global preprocess_handles;
global preprocess_sel_vox;
global FID_FD_copy;
global met_tick;

for i = 1:details_para.num_vox
    [ fit_sig,amp,Ind_sig,freq] = HSVDPK_new(data_value.FID_TD(:,i),details_para.fres,details_para.Fs,10);
    ind=find((freq>=-80 & freq<=80));
    half_sig2=sum(Ind_sig(:,ind),2);
    data_value.FID_TD(:,i)=data_value.FID_TD(:,i)- half_sig2;
    data_value.FID_FD(:,i)=fftshift(fft(data_value.FID_TD(:,i)));
end
PE=details_para.PE;
RE=details_para.RE;
if round(details_para.Fs)==2000
    seg_NAA=160:180;
    seg_Lac=130:159;%210:220;
    seg_Cho=209:218;%205:210;%200:207;%205:210;%205:210;%205:210;%200:205;198:210;
    seg_Cr=196:205;%200:205;%200:205;198:210;
elseif round(details_para.Fs)==1600
    seg_NAA=135:165;
    seg_Lac=110:130;
    seg_Cho=185:200;
elseif round(details_para.Fs)==1500
    seg_NAA=60:95;
    seg_Lac=30:50;
    seg_Cho=145:165;
    seg_Cr=145:165;
elseif round(details_para.Fs)==2500
    seg_NAA=186:200;
    seg_Lac=165:185;
    seg_Cho=210:225;
end
for i=1:(PE/1)*RE
    temp1=abs(data_value.FID_FD(:,i));
    NAA_CSI(i)=max(temp1(seg_NAA))/1;
    Lac_CSI(i)=max(temp1(seg_Lac))/1;
    Cho_CSI(i)=max(temp1(seg_Cho))/1;
    Cr_CSI(i)=max(temp1(seg_Cr));

end
l=1;
for i=1:PE/1
    for j=1:RE
        NAA_map_csi(i,j)=NAA_CSI(l);
        Lac_map_csi(i,j)=Lac_CSI(l);
        Cho_map_csi(i,j)=Cho_CSI(l);
        Cr_map_csi(i,j)=Cr_CSI(l);
        l=l+1;
    end
end
met_tick=1;
preprocess_handles.met_fig=figure();
set(preprocess_handles.met_fig,'Name','RAPID 1.0','NumberTitle', 'off');
subplot(2,2,1)
imagesc(NAA_map_csi);%./Cr_map_csi);
title('NAA map');
colorbar
subplot(2,2,2)
imagesc(Lac_map_csi);%./Cr_map_csi);
title('Lac map');
colorbar
subplot(2,2,3)
imagesc(Cho_map_csi);%./Cr_map_csi);
title('Cho map');
colorbar
subplot(2,2,4)
imagesc(Cr_map_csi);%./Cr_map_csi);
title('Cr map');
colorbar

%------------------------------------------------------------------------------------------------------------------------------
function display_seg_mim_Callback(hObject, eventdata)

global details_para;

temp = get(hObject,'string');
temp1 = str2double(temp);
if(temp1<details_para.ppm_referenced(1) || temp1>details_para.ppm_referenced(details_para.N) || (temp1>details_para.disp_ppm_seg(2))||(isnan(temp1)))
    errordlg('Invalid value');
    set(hObject,'string',num2str(details_para.disp_ppm_seg(1)))
else
    details_para.disp_ppm_seg(1) = temp1; %in ppm
    details_para.disp_freq_seg(1) = ((details_para.disp_ppm_seg(1) - details_para.ref)*details_para.Tf)/1E6; %ppm to frequency
    details_para.disp_seg(1) = floor((details_para.disp_freq_seg(1)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
    
    display_data();
end
%------------------------------------------------------------------------------------------------------------------------------
function display_seg_max_Callback(hObject, eventdata)
global details_para;

temp = get(hObject,'string');
temp1 = str2double(temp); 
if(temp1<details_para.ppm_referenced(1) || temp1>details_para.ppm_referenced(details_para.N) || (temp1<details_para.disp_ppm_seg(1))||(isnan(temp1)))
    errordlg('Invalid value');
    set(hObject,'string',num2str(details_para.disp_ppm_seg(2)))
else
    details_para.disp_ppm_seg(2) = temp1;%ppm
    details_para.disp_freq_seg(2) = ((details_para.disp_ppm_seg(2) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
    details_para.disp_seg(2) = ceil((details_para.disp_freq_seg(2)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
    
    display_data();
end
%------------------------------------------------------------------------------------------------------------------------------