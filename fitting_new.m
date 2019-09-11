function[]= fitting_new()

% This function is used for fitting the data: HLSVD based 


% Author: Sourav Bhaduri
% University of Ghent
% email: Sourav.Bhaduri@UGent.be
% Feb 2018; Last revision: Feb 2018

global data_value;
global fitting_handles;
global details_para;
global apod_para;
global fit_para;
global ind_allvox;
global pha_cor_data;
global default_parameters;
global zerofill_value;
global info_struct;
global cur_vox;
global fitting_display;
global FID_FD_copy;
global met_tick;

FID_FD_copy = data_value.FID_FD;
default_parameters.pha_0=1.9;
default_parameters.pha_1=2.1;
default_parameters.pha_0_NAA=1.9;
default_parameters.pha_1_NAA=2.3;
default_parameters.pha_0_Cr=2.9;
default_parameters.pha_1_Cr=3.25;
default_parameters.pha_0_Cho=3.3;
default_parameters.pha_1_Cho=3.45;
default_parameters.pha_0_Lac=1.2;
default_parameters.pha_1_Lac=1.6;
default_parameters.pha_0_mi=3.5;
default_parameters.pha_1_mi=3.7;
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
details_para.disp_fit(1)=default_parameters.pha_0;
details_para.disp_fit(2)=default_parameters.pha_1;
details_para.fit_flag=0;
details_para.fit_method=0;
if(isempty(apod_para)) % if pre-processing is not done previously load default pre-processing parameters
    
    apod_para.type = default_parameters.apod_para_type*ones(details_para.num_vox,1);
    apod_para.lor = default_parameters.apod_para_lor*ones(details_para.num_vox,1);
    apod_para.gaus = default_parameters.apod_para_gaus*ones(details_para.num_vox,1);
    apod_para.sigm = default_parameters.apod_para_sigm*ones(details_para.num_vox,1);
    apod_para.function_weight = ones(size(data_value.FID_FD));
    fit_para.pha_0 = default_parameters.pha_0*ones(details_para.num_vox,1);
    fit_para.pha_1 = default_parameters.pha_1*ones(details_para.num_vox,1);
    pha_cor_data =  data_value.FID_TD;
    fit_para.order=15*ones(details_para.num_vox,1);
else % if pre-processing is done use the current pre-processing parameters 
    % and perform pre-processing using those values 
    for i = 1:length(details_para.selected_voxels)
        vox = details_para.selected_voxels(1,i);
        apod_para.type=1*ones(details_para.num_vox,1);
        fit_para.pha_0 = default_parameters.pha_0*ones(details_para.num_vox,1);
        fit_para.pha_1 = default_parameters.pha_1*ones(details_para.num_vox,1);
        fit_para.pha_0_NAA = default_parameters.pha_0_NAA*ones(details_para.num_vox,1);
        fit_para.pha_1_NAA = default_parameters.pha_1_NAA*ones(details_para.num_vox,1);
        fit_para.pha_0_Cr = default_parameters.pha_0_Cr*ones(details_para.num_vox,1);
        fit_para.pha_1_Cr = default_parameters.pha_1_Cr*ones(details_para.num_vox,1);
        fit_para.pha_0_Cho = default_parameters.pha_0_Cho*ones(details_para.num_vox,1);
        fit_para.pha_1_Cho = default_parameters.pha_1_Cho*ones(details_para.num_vox,1);
        fit_para.pha_0_Lac = default_parameters.pha_0_Lac*ones(details_para.num_vox,1);
        fit_para.pha_1_Lac = default_parameters.pha_1_Lac*ones(details_para.num_vox,1);
        fit_para.pha_0_mi = default_parameters.pha_0_mi*ones(details_para.num_vox,1);
        fit_para.pha_1_mi = default_parameters.pha_1_mi*ones(details_para.num_vox,1);
        fit_para.order=15*ones(details_para.num_vox,1);
    end
end
details_para.fit_freq_NAA(1) = ((fit_para.pha_0_NAA(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_NAA(2) = ((fit_para.pha_1_NAA(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_Cr(1) = ((fit_para.pha_0_Cr(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_Cr(2) = ((fit_para.pha_1_Cr(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_Cho(1) = ((fit_para.pha_0_Cho(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_Cho(2) = ((fit_para.pha_1_Cho(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_Lac(1) = ((fit_para.pha_0_Lac(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_Lac(2) = ((fit_para.pha_1_Lac(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_mi(1) = ((fit_para.pha_0_mi(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_mi(2) = ((fit_para.pha_1_mi(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_disp_seg_NAA(1) = floor((details_para.fit_freq_NAA(1)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_NAA(2) = ceil((details_para.fit_freq_NAA(2)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_Cr(1) = floor((details_para.fit_freq_Cr(1)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_Cr(2) = ceil((details_para.fit_freq_Cr(2)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_Cho(1) = floor((details_para.fit_freq_Cho(1)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_Cho(2) = ceil((details_para.fit_freq_Cho(2)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_Lac(1) = floor((details_para.fit_freq_Lac(1)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_Lac(2) = ceil((details_para.fit_freq_Lac(2)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_mi(1) = floor((details_para.fit_freq_mi(1)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_mi(2) = ceil((details_para.fit_freq_mi(2)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
x_max = 850;
y_max = 600;
selected_vox = details_para.selected_voxels;
vox_no = cell(1,length(selected_vox));
for i = 1:length(selected_vox)
vox_no{i} =  num2str(selected_vox(1,i));
end


% Cerate GUI
fitting_handles.fig = figure('Visible','off','Position',[100,100,x_max,y_max],'Color',[0.68,0.92,1],'Name','RAPID 1.0','NumberTitle', 'off');
set(fitting_handles.fig,'MenuBar','none'); 
% fitting_handles.save_figure = uimenu(fitting_handles.fig,'Label','Save Figure','Callback',@save_preprocess_fig_callback);

fitting_handles.h1 = uipanel(fitting_handles.fig,'Title','Voxel No','Units','normalized','Position',[15/x_max 190/y_max 70/x_max 270/y_max],...
    'BackgroundColor','white','FontSize',8, 'FontWeight','Bold');
% fitting_handles.voxel_list = uicontrol(fitting_handles.h1,'Style', 'listbox','Units','normalized','Position', [.05 .027 .9 .97],'String', vox_no,...
%     'Max', details_para.num_vox, 'Min', 1, 'Callback', @voxel_list_callback, 'Interruptible','off','BusyAction','cancel'); 
fitting_handles.voxel_list = uicontrol(fitting_handles.h1,'Style', 'listbox','Units','normalized','Position', [.05 .027 .9 .97],'String', vox_no,...
    'value',1, 'Callback', @voxel_list_callback, 'Interruptible','off','BusyAction','cancel'); 


fitting_handles.axis = axes('Units','Pixels','Units','normalized','Position',[120/x_max,200/y_max,420/x_max,250/y_max]);
% fitting_handles.axis_text = uicontrol(fitting_handles.fig,'Style','text','Units','normalized','Position',[220/x_max 460/y_max 130/x_max 15/y_max],...
%     'String','Signal Spectrum','BackgroundColor',[0.68,0.92,1],'FontSize',9,'FontWeight','Bold');

fitting_handles.TF_FD_display = uibuttongroup(fitting_handles.fig,'Units','normalized','Position',[0.1 0.92 .4 0.06],...
    'FontSize',9,'BackgroundColor',[0.68,0.92,1],'BorderType','none','SelectionChangeFcn',@TD_FD_display_callback);
fitting_handles.FD_disp_real = uicontrol(fitting_handles.TF_FD_display,'Style','radiobutton','Units','normalized','Position',[.01 .1 .3 .8],'String',...
    'FD signal real','FontWeight','Bold','BackgroundColor',[0.68,0.92,1]);
fitting_handles.FD_disp_abs = uicontrol(fitting_handles.TF_FD_display,'Style','radiobutton','Units','normalized','Position',[.35 .1 .3 .8],'String',...
    'FD signal abs','FontWeight','Bold','BackgroundColor',[0.68,0.92,1]);
fitting_handles.TD_disp = uicontrol(fitting_handles.TF_FD_display,'Style','radiobutton','Units','normalized','String','TD signal','Position',...
    [.71 .1 .3 .8],'FontWeight','Bold','BackgroundColor',[0.68,0.92,1]);


fitting_handles.display_panel = uipanel(fitting_handles.fig,'Title','Display Properties','Position',[.1 .8 .27 .1],'BackgroundColor','white','FontSize',9,'FontWeight','Bold');
fitting_handles.seg_min_text = uicontrol(fitting_handles.display_panel,'Style','text','Units','normalized','Position',[0.4,0.7,0.2,0.1*2.5],'String',...
    'Min','FontSize',8,'FontWeight','Bold','BackgroundColor','white');
fitting_handles.seg_max_text = uicontrol(fitting_handles.display_panel,'Style','text','Units','normalized','Position',[0.7,0.7,0.2,0.1*2.5],'String',...
    'Max','FontSize',8,'FontWeight','Bold','BackgroundColor','white');
fitting_handles.display_seg_mim = uicontrol(fitting_handles.display_panel,'Style','edit','Units','normalized','Position',[0.45,0.2,0.15,0.15*2],'String','0',...
    'Callback',{@display_seg_mim_Callback},'Interruptible','off','BusyAction','cancel'); 
fitting_handles.display_seg_max = uicontrol(fitting_handles.display_panel,'Style','edit','Units','normalized','Position',[0.75,0.2,0.15,0.15*2],'String','0',...
    'Callback',{@display_seg_max_Callback},'Interruptible','off','BusyAction','cancel'); 
% fitting_handles.zerofill_text = uicontrol(fitting_handles.fig,'Style','text','Units','normalized','Position',[557/x_max 460/y_max 70/x_max 15/y_max],...
%     'String','Zerofill','BackgroundColor',[0.68,0.92,1],'FontSize',9,'FontWeight','Bold');
% fitting_handles.zerofill_options = uicontrol(fitting_handles.fig,'Style','popup','String','No Fill|1 times|2 times|3 times','Units','normalized',...
%     'Position',[555/x_max 405/y_max 80/x_max 15/y_max],'Value',zerofill_value,'Interruptible','off','BusyAction','cancel','Callback', @zerofill_options_callback); 
fitting_handles.file_panel = uipanel(fitting_handles.fig,'Title','','Position',[.5 .9 .43 .1],'BackgroundColor','white','FontSize',5,'FontWeight','Bold');
fitting_handles.main_text = uicontrol(fitting_handles.file_panel,'Style','text','Units','normalized','Position',[0.005,0.65,0.05,0.1*2.5],'String',...
    'M','FontSize',6,'FontWeight','Bold','BackgroundColor','white');
fitting_handles.ref_text = uicontrol(fitting_handles.file_panel,'Style','text','Units','normalized','Position',[0.005,0.1,0.05,0.1*2.5],'String',...
    'R','FontSize',6,'FontWeight','Bold','BackgroundColor','white');
fitting_handles.main_val = uicontrol(fitting_handles.file_panel,'Style','edit','Units','normalized','Position',[0.05,0.65,0.95,0.15*2],'String',data_value.file_name{1},...
    'Callback',{@file_main_Callback},'Interruptible','off','BusyAction','cancel'); 
fitting_handles.ref_val = uicontrol(fitting_handles.file_panel,'Style','edit','Units','normalized','Position',[0.05,0.1,0.95,0.15*2],'String',data_value.file_name{2},...
    'Callback',{@file_ref_Callback},'Interruptible','off','BusyAction','cancel'); 

% apodization panel
fitting_handles.h2 = uipanel(fitting_handles.fig,'Title','Metabolite','Units','normalized','Position',[600/x_max 370/y_max 150/x_max 100/y_max],...
    'BackgroundColor','white','FontSize',9, 'FontWeight','Bold');
fitting_handles.metabolite = uicontrol(fitting_handles.h2,'Style','text','Units','normalized','String','Metabolite Type','Units','normalized',...
    'Position',[0.1 0.65 0.8 0.3],'BackgroundColor','white','FontSize',8,'FontWeight','Bold');
fitting_handles.metabolite_type = uicontrol(fitting_handles.h2,'Style','popup','String','NAA|Lactate|Creatine|Choline|Myo','Units','normalized','Position',...
    [0.1 0.19 0.8 0.3],'Value',apod_para.type(cur_vox), 'Callback',{@metabolite_type_Callback});


% phase correction panel
fitting_handles.h3 = uipanel(fitting_handles.fig,'Title','Fitting range','Units','normalized','Position',[20/x_max 30/y_max 520/x_max 150/y_max],...
    'BackgroundColor','white','FontSize',9, 'FontWeight','Bold');
fitting_handles.pha_corr0_text = uicontrol(fitting_handles.h3,'Style','text','Units','normalized','Position',[0.1 0.5 0.2 .2],'String',...
    'start ','FontSize',8,'FontWeight','Bold','BackgroundColor','white');
fitting_handles.pha_corr1_text = uicontrol(fitting_handles.h3,'Style','text','Units','normalized','Position',[0.1 0.2 0.2 0.2],'String',...
    'end','FontSize',8,'FontWeight','Bold','BackgroundColor','white');
fitting_handles.pha_corr_value_text = uicontrol(fitting_handles.h3,'Style','text','Units','normalized','Position',[0.31 0.71 0.2 .2],'String',...
    'Value (PPM)','FontSize',8,'FontWeight','Bold','BackgroundColor','white');
fitting_handles.pha_corr0_value = uicontrol(fitting_handles.h3,'Style','edit','Units','normalized','Position',[0.34 0.5 .15 .2],...
    'String',num2str(fit_para.pha_0(cur_vox)),'Callback',{@pha0_val_Callback});
fitting_handles.pha_corr1_value = uicontrol(fitting_handles.h3,'Style','edit','Units','normalized','Position',[0.34 0.2 .15 .2],...
    'String',num2str(fit_para.pha_1(cur_vox)),'Callback',{@pha1_val_Callback});
fitting_handles.order_text = uicontrol(fitting_handles.h3,'Style','text','Units','normalized','Position',[0.71 0.7 0.2 .2],'String',...
    'order: ','FontSize',8,'FontWeight','Bold','BackgroundColor','white');
fitting_handles.order_value = uicontrol(fitting_handles.h3,'Style','edit','Units','normalized','Position',[0.85 0.75 .1 .2],...
    'String',num2str(fit_para.order(cur_vox)),'Callback',{@order_val_Callback});


fitting_handles.close = uicontrol(fitting_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[680/x_max 10/y_max 100/x_max 25/y_max],...
'String', 'OK/Close','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@close_Callback);


fitting_handles.fit = uicontrol(fitting_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[300/x_max 130/y_max 100/x_max 25/y_max],...
'String', 'Fit','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@fitting_Callback);
fitting_handles.map = uicontrol(fitting_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[300/x_max 80/y_max 100/x_max 25/y_max],...
'String', 'Map','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@map_Callback);
fitting_handles.display_grid = uicontrol(fitting_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[300/x_max 40/y_max 100/x_max 25/y_max],...
'String', 'Display grid','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@display_grid_Callback);

% fitting_handles.wavelet = uicontrol(fitting_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[240/x_max 140/y_max 100/x_max 25/y_max],...
% 'String', 'Denoise','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@denoise_Callback);


set(fitting_handles.fig,'Visible','on');
set(fitting_handles.display_seg_mim, 'Enable', 'on','string',num2str(details_para.disp_ppm_seg(1)))
set(fitting_handles.display_seg_max, 'Enable', 'on','string',num2str(details_para.disp_ppm_seg(2)))
set(fitting_handles.pha_corr0_value, 'Enable', 'on','string',num2str(details_para.disp_fit(1)))
set(fitting_handles.pha_corr1_value, 'Enable', 'on','string',num2str(details_para.disp_fit(2)))
set(fitting_handles.map , 'Enable', 'off');
set(fitting_handles.display_grid , 'Enable', 'off');

set(fitting_handles.main_val, 'Enable', 'off')
set(fitting_handles.ref_val, 'Enable', 'off')

 
fitting_display = 1;

% Display the selected voxel 
display_data()
waitfor(fitting_handles.fig);

function voxel_list_callback(hObject, eventdata)
% Callback for 'voxel no.' get the selected voxel and update the plot according to the voxel selected
global details_para;
global fitting_handles;
global apod_para;
global fit_para;
global cur_vox;

% get the selected voxel
preprocess_sel_vox = get(hObject,'Value');

if(isempty(preprocess_sel_vox))
    preprocess_sel_vox = 1;
    set(fitting_handles.voxel_list,'value',1);
end

selected_vox = details_para.selected_voxels;
cur_vox = selected_vox(1,preprocess_sel_vox);
display_data()


% Update the GUI based on the selected voxel
temp = apod_para.type(cur_vox);
set(fitting_handles.metabolite_type,'value',temp);

set(fitting_handles.pha_corr0_value,'string',num2str(fit_para.pha_0(cur_vox)));
set(fitting_handles.pha_corr1_value,'string',num2str(fit_para.pha_1(cur_vox)));
% set(fitting_handles.display_seg_mim, 'Enable', 'on','string',num2str(details_para.disp_ppm_seg(1)))
% set(fitting_handles.display_seg_max, 'Enable', 'on','string',num2str(details_para.disp_ppm_seg(2)))



function metabolite_type_Callback(hObject, eventdata)
%Callback for 'Apodisation type' 
global fitting_handles;
global apod_para;
global ind_allvox;
global details_para;
global cur_vox;
global fit_para;
% global data_value;

% get the selected apodization type
temp = get(hObject,'Value');
apod_para.type(details_para.selected_voxels(1,:)) = temp;

% enable/disable based on the selected type
switch(temp)
  case 1 %NAA
        phi_0 = 1.9;
        fit_para.pha_0(details_para.selected_voxels(1,:)) = phi_0;
        phi_1 = 2.1;
        fit_para.pha_1(details_para.selected_voxels(1,:)) = phi_1;    
        set(fitting_handles.pha_corr0_value, 'Enable', 'on','string',num2str(phi_0));
        set(fitting_handles.pha_corr1_value, 'Enable', 'on','string',num2str(phi_1));
  case 2 %Lactate
         phi_0 = 1.1;
         fit_para.pha_0(details_para.selected_voxels(1,:)) = phi_0;    
         phi_1 = 1.4;
         fit_para.pha_1(details_para.selected_voxels(1,:)) = phi_1;    
        set(fitting_handles.pha_corr0_value, 'Enable', 'on','string',num2str(phi_0));
        set(fitting_handles.pha_corr1_value, 'Enable', 'on','string',num2str(phi_1));
  case 3 %Creatine
         phi_0 = 2.8;
         fit_para.pha_0(details_para.selected_voxels(1,:)) = phi_0;    
         phi_1 = 3.1;
         fit_para.pha_1(details_para.selected_voxels(1,:)) = phi_1;    
        set(fitting_handles.pha_corr0_value, 'Enable', 'on','string',num2str(phi_0));
        set(fitting_handles.pha_corr1_value, 'Enable', 'on','string',num2str(phi_1));
  case 4 %Choline
         phi_0 = 3.3;
         fit_para.pha_0(details_para.selected_voxels(1,:)) = phi_0;    
         phi_1 = 3.4;
         fit_para.pha_1(details_para.selected_voxels(1,:)) = phi_1;    
         set(fitting_handles.pha_corr0_value, 'Enable', 'on','string',num2str(phi_0));
         set(fitting_handles.pha_corr1_value, 'Enable', 'on','string',num2str(phi_1));
  case 5 %Myo inositol
         phi_0 = 3.4;
         fit_para.pha_0(details_para.selected_voxels(1,:)) = phi_0;    
         phi_1 = 3.65;
         fit_para.pha_1(details_para.selected_voxels(1,:)) = phi_1;    
         set(fitting_handles.pha_corr0_value, 'Enable', 'on','string',num2str(phi_0));
         set(fitting_handles.pha_corr1_value, 'Enable', 'on','string',num2str(phi_1));
end
 
% 
function pha0_val_Callback(hObject, eventdata)
%Callback for 'zero-order' edit box 
global fit_para;
global cur_vox;
global ind_allvox;
global fitting_handles;
global details_para;

% read and check the zero-order phase value from the edit box 
phi_0 = str2double(get(hObject,'string'));
if((phi_0 < -10)|| (phi_0 > 10) || (isnan(phi_0)))
    errordlg('Invalid value');
    set(hObject,'string',num2str(fit_para.pha_0(cur_vox)))
    return;
end
fit_para.pha_0(details_para.selected_voxels(1,:)) = phi_0;   
set(fitting_handles.pha_corr0_value,'string',num2str(fit_para.pha_0(cur_vox)));

% 
function pha1_val_Callback(hObject, eventdata)
%Callback for 'first-order' edit box 
global fit_para;
global cur_vox;
global ind_allvox;
global fitting_handles;
global details_para;

% read and check the first-order phase value from the edit box 
phi_1 = str2double(get(hObject,'string'));
if((phi_1 < -10)|| (phi_1 > 10) || (isnan(phi_1)))
    errordlg('Invalid value');
    set(hObject,'string',num2str(fit_para.pha_1(cur_vox)))
    return;
end
fit_para.pha_1(details_para.selected_voxels(1,:)) = phi_1; 
set(fitting_handles.pha_corr1_value,'string',num2str(fit_para.pha_1(cur_vox)));

function order_val_Callback(hObject, eventdata)
%Callback for 'first-order' edit box 
global fit_para;
global cur_vox;
global ind_allvox;
global fitting_handles;
global details_para;

% read and check the first-order phase value from the edit box 
ord_1 = str2double(get(hObject,'string'));
if((ord_1 < 0)|| (ord_1 > 100) || (isnan(ord_1)))
    errordlg('Invalid value');
    set(hObject,'string',num2str(fit_para.order(cur_vox)))
    return;
end
fit_para.order(details_para.selected_voxels(1,:)) = ord_1; 
set(fitting_handles.order_value,'string',num2str(fit_para.order(cur_vox)));



function TD_FD_display_callback(hObject, eventdata)
global fitting_display;

temp = get(get(hObject,'SelectedObject'),'String');

if strcmp(temp,'FD signal real')
    fitting_display = 1;
elseif strcmp(temp,'FD signal abs')
    fitting_display = 2;
else
    fitting_display = 0;
end
display_data()

function display_data()
global fitting_display;
global fitting_handles;
global cur_vox;
global details_para;
global data_value;
global apod_para;

seg = details_para.disp_seg;
ppm = details_para.ppm_referenced;

if (details_para.fit_flag==0)
    if (fitting_display == 1)
        cla(fitting_handles.axis);
        plot_var = real(data_value.FID_FD(seg(1):seg(2),cur_vox));
        plot(fitting_handles.axis,ppm(seg(1):seg(2))',plot_var);
        max_p = max(plot_var);
        min_p = min(plot_var);
        axis(fitting_handles.axis,[ppm(seg(1)), ppm(seg(2)),min_p - 0*max_p, max_p + 0.0*max_p])
        set(fitting_handles.axis,'XDir','reverse');
    elseif (fitting_display == 2)
        cla(fitting_handles.axis);
        plot_var = abs(data_value.FID_FD(seg(1):seg(2),cur_vox));
        plot(fitting_handles.axis,ppm(seg(1):seg(2))',plot_var);
        max_p = max(plot_var);
        min_p = min(plot_var);
        axis(fitting_handles.axis,[ppm(seg(1)), ppm(seg(2)),min_p - 0.5*max_p, max_p + 0.5*max_p])
        set(fitting_handles.axis,'XDir','reverse');
    else
        cla(fitting_handles.axis);
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
        plot(fitting_handles.axis,details_para.t',plot_var);
    %     hold on;
    %     plot(fitting_handles.axis,details_para.t',filter,'r');
    %     axis(fitting_handles.axis,[details_para.t(1), details_para.t(details_para.N),min_p , max_p])
    %     hold off
    end
else
     if (fitting_display == 1)
        cla(fitting_handles.axis);
        plot_var = real(data_value.FID_FD(seg(1):seg(2),cur_vox));
        plot_fit = real(data_value.fit_FD(seg(1):seg(2),cur_vox));
        plot(fitting_handles.axis,ppm(seg(1):seg(2))',plot_var);
        hold on 
        plot(fitting_handles.axis,ppm(seg(1):seg(2))',plot_fit);
        max_p = max(plot_var);
        min_p = min(plot_var);
        axis(fitting_handles.axis,[ppm(seg(1)), ppm(seg(2)),min_p - 0.0*max_p, max_p + 0.0*max_p])
        set(fitting_handles.axis,'XDir','reverse');
        legend(fitting_handles.axis,'Original','Fitted');
    elseif (fitting_display == 2)
        cla(fitting_handles.axis);
        plot_var = abs(data_value.FID_FD(seg(1):seg(2),cur_vox));
        plot_fit = abs(data_value.fit_FD(seg(1):seg(2),cur_vox));
        plot(fitting_handles.axis,ppm(seg(1):seg(2))',plot_var);
        hold on 
        plot(fitting_handles.axis,ppm(seg(1):seg(2))',plot_fit);
        max_p = max(plot_var);
        min_p = min(plot_var);
        axis(fitting_handles.axis,[ppm(seg(1)), ppm(seg(2)),min_p - 0.1*max_p, max_p + 0.1*max_p])
        set(fitting_handles.axis,'XDir','reverse');
        legend(fitting_handles.axis,'Original','Fitted');
    else
        cla(fitting_handles.axis);
        plot_var = real(data_value.FID_TD(:,cur_vox));
        plot_fit = real(data_value.fit_TD(:,cur_vox));
    %     plot_var= ifftshift(real(data_value.FID_TD(:,cur_vox))+flipud(real(data_value.FID_TD(:,cur_vox))));
    %     filter = apod_para.function_weight(:,cur_vox);
        max_p = max(abs(plot_var));
        min_p = min(abs(plot_var));
    %     if (max(filter)> max_p)
    %         max_p = max(filter);
    %         min_p = -max_p;
    %     else
    %         max_p = max(plot_var);
    %         min_p = -max_p;
    %     end   
        axis(fitting_handles.axis,[details_para.t(1), details_para.t(end),min_p - 0.5*max_p, max_p + 0.1*max_p])
        plot(fitting_handles.axis,details_para.t',plot_var);
        hold on 
        plot(fitting_handles.axis,details_para.t',plot_fit);
        legend(fitting_handles.axis,'Original','Fitted');
    %     hold on;
    %     plot(fitting_handles.axis,details_para.t',filter,'r');
    %     axis(fitting_handles.axis,[details_para.t(1), details_para.t(details_para.N),min_p , max_p])
    %     hold off
    end
end

function save_preprocess_fig_callback(hObject, eventdata)
%Callback for 'save figure'
global fitting_handles;
global log_file;
ini_path = log_file.saving_path;
[temp,s_path] = uiputfile(ini_path);
if(temp == 0)
    return;
else
    dot_pos = strfind(temp, '.');
    s_fln = temp(1:dot_pos-1);
    file_name = strcat(s_path,s_fln);
    saveas(fitting_handles.fig,file_name,'png');
    log_file.saving_path = s_path;
end



function close_Callback(hObject, eventdata)
%Callback for 'close' button
global fitting_handles;
global met_tick;
global data_value;
if (met_tick==1 && ~isempty(fitting_handles.met_fig))
    close(fitting_handles.fig)
    close(fitting_handles.met_fig);
else
     close(fitting_handles.fig);
end
    

function file_main_Callback(hObject, eventdata)
%Callback for 'first-order' edit box 
global pha_para;
global cur_vox;
global ind_allvox;
global fitting_handles;
global details_para;
 

function file_ref_Callback(hObject, eventdata)
%Callback for 'first-order' edit box 
global pha_para;
global cur_vox;
global ind_allvox;
global fitting_handles;
global details_para;
% 
%-----------------------------------------------------------------------------------------------------------------------------------
function fitting_Callback(hObject, eventdata)

global data_value;
global details_para;
global fitting_handles;
global preprocess_sel_vox;
global FID_FD_copy;
global fit_para;
global cur_vox;
global ind_allvox;

details_para.fit_flag=1;

set(fitting_handles.display_grid, 'Enable', 'off');
set(fitting_handles.map, 'Enable', 'off');

details_para.fit_freq(1) = ((fit_para.pha_0(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq(2) = ((fit_para.pha_1(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_NAA(1) = ((fit_para.pha_0_NAA(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_NAA(2) = ((fit_para.pha_1_NAA(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_Cr(1) = ((fit_para.pha_0_Cr(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_Cr(2) = ((fit_para.pha_1_Cr(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_Cho(1) = ((fit_para.pha_0_Cho(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_Cho(2) = ((fit_para.pha_1_Cho(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_Lac(1) = ((fit_para.pha_0_Lac(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_Lac(2) = ((fit_para.pha_1_Lac(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_mi(1) = ((fit_para.pha_0_mi(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_mi(2) = ((fit_para.pha_1_mi(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_disp_seg(1) = floor((details_para.fit_freq(1)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg(2) = ceil((details_para.fit_freq(2)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_NAA(1) = floor((details_para.fit_freq_NAA(1)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_NAA(2) = ceil((details_para.fit_freq_NAA(2)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_Cr(1) = floor((details_para.fit_freq_Cr(1)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_Cr(2) = ceil((details_para.fit_freq_Cr(2)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_Cho(1) = floor((details_para.fit_freq_Cho(1)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_Cho(2) = ceil((details_para.fit_freq_Cho(2)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_Lac(1) = floor((details_para.fit_freq_Lac(1)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_Lac(2) = ceil((details_para.fit_freq_Lac(2)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_mi(1) = floor((details_para.fit_freq_mi(1)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_mi(2) = ceil((details_para.fit_freq_mi(2)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
h = waitbar(0,'Please wait...');
if details_para.fit_method==0
    for i = 1:details_para.num_vox
        
%         load('data_baseline_19_3_Shifted.mat');
%         if i~=1
%             data_value.FID_FD(:,i)=signals_shifted(:,i-1);
%             data_value.FID_TD(:,i)=ifft(ifftshift( data_value.FID_FD(:,i)));
%         end
%         data_value.FID_FD(:,i)=data_value.FID_FD(:,i)+ abs(min(data_value.FID_FD(:,i)));
%         data_value.FID_TD(:,i)=ifft(ifftshift(data_value.FID_FD(:,i)));
%         sp1 = csaps(ppm,real(data_value.FID_FD(:,i)),30);
%         re_bl =fnval(sp1,ppm);
        re_bl = ssa(real(data_value.FID_FD(:,i)),40); % estimate the baseline for real part
%         sp2 = csaps(ppm,imag(data_value.FID_FD(:,i)),30);
%         im_bl = fnval(sp2,ppm);
        im_bl = ssa(imag(data_value.FID_FD(:,i)),40); % estimate the baseline for imaginary part
        data_value.FID_FD_base(:,i)=re_bl + 1i*im_bl;
        data_value.FID_FD_b(:,i)=data_value.FID_FD(:,i)-data_value.FID_FD_base(:,i);
        data_value.FID_TD_b(:,i)=ifft(ifftshift(data_value.FID_FD_b(:,i)));
        
%         details_para.fla=0;
%         details_para.seg = [1 length(data_value.FID_TD(:,i))];
%         if i==1
% %              [ fit_sig,amp,Ind_sig,freq,Z_sig,FWHM] = HSVDPK_new(data_value.FID_TD(:,i),details_para.fres,details_para.Fs,fit_para.order(cur_vox));
%             [No_pks,pk_loc,FWHM_res,signal_power_each] = find_peaks_new(data_value.FID_TD(:,3),details_para.fres);
%         end
%         freq=details_para.fres(pk_loc);
%         re_err1 = residual_calc(real(data_value.FID_FD(:,i)), 1024, 4,freq,FWHM_res,pk_loc);
%         im_err1 = residual_calc(imag(data_value.FID_FD(:,i)), 1024, 4,freq,FWHM_res,pk_loc);
%         err1 = re_err1 + 1i*im_err1;
%         data_value.FID_FD_base(:,i)=err1;
%         data_value.FID_FD_b(:,i)=data_value.FID_FD(:,i)-data_value.FID_FD_base(:,i);
%         data_value.FID_TD_b(:,i)=ifft(ifftshift(data_value.FID_FD_b(:,i)));
        
%         data_value.FID_TD(:,i)=data_value.FID_TD_b(:,i);
%         data_value.FID_FD(:,i)=fftshift(fft(data_value.FID_TD(:,i)));
        baseline(:,1)=data_value.FID_FD_base(:,i);
        baseline_t(:,1)=ifft(ifftshift(data_value.FID_FD_base(:,i)));
        [ fit_sig,amp,Ind_sig,freq,Z_sig,FWHM] = HSVDPK_new(data_value.FID_TD_b(:,i),details_para.fres,details_para.Fs,fit_para.order(cur_vox));
        signal_fitted(:,1)=fftshift(fft(fit_sig));
        count=2;
        while(1)
%                 sp1 = csaps(ppm,real(data_value.FID_FD(:,i)-signal_fitted(:,count-1)),30);
%                 re_bl =fnval(sp1,ppm);
                re_bl = ssa(real(data_value.FID_FD(:,i)-signal_fitted(:,count-1)),30); % estimate the baseline for real part
%                 sp2 = csaps(ppm,imag(data_value.FID_FD(:,i)-signal_fitted(:,count-1)),30);
%                 im_bl =fnval(sp2,ppm);
                im_bl = ssa(imag(data_value.FID_FD(:,i)-signal_fitted(:,count-1)),30); % estimate the baseline for imaginary part
                baseline(:,count)=re_bl + 1i*im_bl;
                baseline_t(:,count)=ifft(ifftshift(baseline(:,count)));
                [ fit_sig,amp,Ind_sig,freq,Z_sig,FWHM] = HSVDPK_new(data_value.FID_TD(:,i)-baseline_t(:,count),details_para.fres,details_para.Fs,fit_para.order(cur_vox));
                signal_fitted(:,count)=fftshift(fft(fit_sig));
                if (std(abs(baseline(:,count)-baseline(:,count-1)))<=1e-5 || count==5)
                    break;
                end
                count=count+1;
        end
%         [ fit_sig,amp,Ind_sig,freq,Z_sig,FWHM] = HSVDPK_new(data_value.FID_TD_b(:,i),details_para.fres,details_para.Fs,fit_para.order(cur_vox));
         model_select=1;
        Initx = [amp (FWHM)' freq'];
    %     lb= [(amp-(amp/4)) ((FWHM)-50)' (freq-20)'];
    %     ub= [(amp+(amp/4)) ((FWHM)+50)' (freq+20)'];
    %      options = optimset('MaxFunEvals', 500, 'MaxIter', 20, 'Algorithm','trust-region-reflective', 'Display', 'off'); % optimization parameters 
    %     [fit_param11,~,residual,~,~,~,J] = lsqnonlin(@(x)model_estimate_equation2(x,details_para.t,data_value.FID_TD(:,i),model_select),Initx,[],[],options); % non-linear least squares algorithm
    %     [fit_sig,Ind_sig,Z_sig] = model_estimate(fit_param11,details_para.t,model_select);
        ind=find((freq>=details_para.fit_freq(1)-20 & freq<=details_para.fit_freq(2)+20));
        ind_NAA=find((freq>=details_para.fit_freq_NAA(1)-20 & freq<=details_para.fit_freq_NAA(2)+20));
        ind_Cr=find((freq>=details_para.fit_freq_Cr(1)-20 & freq<=details_para.fit_freq_Cr(2)+20));
        ind_Cho=find((freq>=details_para.fit_freq_Cho(1)-20 & freq<=details_para.fit_freq_Cho(2)+20));
        ind_Myo=find((freq>=details_para.fit_freq_Cr(1)-20 & freq<=details_para.fit_freq_Cr(2)+20));
        ind_Lac=find((freq>=details_para.fit_freq_Cho(1)-20 & freq<=details_para.fit_freq_Cho(2)+20));
        half_sig2=sum(Ind_sig(:,ind),2);
        data_value.fit_TD(:,i)=half_sig2;
        data_value.fit_FD_NAA(:,i)=fftshift(fft(sum(Ind_sig(:,ind_NAA),2)));
        data_value.fit_FD_Cr(:,i)=fftshift(fft(sum(Ind_sig(:,ind_Cr),2)));
        data_value.fit_FD_Cho(:,i)=fftshift(fft(sum(Ind_sig(:,ind_Cho),2)));
        data_value.fit_FD_Myo(:,i)=fftshift(fft(sum(Ind_sig(:,ind_Myo),2)));
        data_value.fit_FD_Lac(:,i)=fftshift(fft(sum(Ind_sig(:,ind_Lac),2)));
        data_value.fit_FD(:,i)=fftshift(fft(data_value.fit_TD(:,i)));
        data_value.fit_FD_temp(:,i)=fftshift(fft(data_value.fit_TD(:,i)));
        data_value.inhom(:,i)=sum(FWHM(ind_NAA),2);
%         min_ref=min(real(data_value.fit_FD(1:30,i)));
%         data_value.fit_FD(:,i)=data_value.fit_FD(:,i);
%         data_value.peak_area(:,i)=trapz(details_para.ppm,real(data_value.fit_FD(:,i)));
        data_value.peak_area(:,i)=trapz(details_para.ppm(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2)),abs(real(data_value.fit_FD(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2),i))));
%         base_temp=real(data_value.FID_FD(:,i)) - real(fftshift(fft(baseline_t(:,count))));
%         data_value.peak_area(:,i)=trapz(details_para.ppm(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2)),abs(real(base_temp(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2)))));
%         data_value.peak_area(:,i)=trapz(details_para.ppm(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2)),abs(data_value.FID_FD(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2),i))-abs(baseline(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2),count-1)));%details_para.fit_disp_seg(1)
%         data_value.peak_area_NAA(:,i)=trapz(details_para.ppm,abs(real(data_value.fit_FD_NAA(:,i))));
%         data_value.peak_area_NAA(:,i)=trapz(details_para.ppm_referenced(details_para.fit_disp_seg_NAA(1):details_para.fit_disp_seg_NAA(2)),abs(real(data_value.fit_FD_NAA(details_para.fit_disp_seg_NAA(1):details_para.fit_disp_seg_NAA(2),i))));
%         data_value.peak_area_NAA(:,i)=max(real(data_value.fit_FD_NAA(:,i)));
        data_value.peak_area_NAA(:,i)=max(abs(data_value.FID_FD(details_para.fit_disp_seg_NAA(1):details_para.fit_disp_seg_NAA(2),i)));
        data_value.peak_area_Cr(:,i)=max(abs(data_value.FID_FD(details_para.fit_disp_seg_Cr(1):details_para.fit_disp_seg_Cr(2),i)));
        data_value.peak_area_Cho(:,i)=max(abs(data_value.FID_FD(details_para.fit_disp_seg_Cho(1):details_para.fit_disp_seg_Cho(2),i)));
        data_value.peak_area_Myo(:,i)=max(abs(data_value.FID_FD(details_para.fit_disp_seg_mi(1):details_para.fit_disp_seg_mi(2),i)));
        data_value.peak_area_Lac(:,i)=max(abs(data_value.FID_FD(details_para.fit_disp_seg_Lac(1):details_para.fit_disp_seg_Lac(2),i)));
%         data_value.peak_area_Cr(:,i)=trapz(details_para.ppm_referenced(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2)),real(data_value.fit_FD_Cr(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2),i)));
%         data_value.peak_area_Cho(:,i)=trapz(details_para.ppm_referenced(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2)),real(data_value.fit_FD_Cho(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2),i)));
%         data_value.peak_area_Myo(:,i)=trapz(details_para.ppm_referenced(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2)),real(data_value.fit_FD_Myo(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2),i)));
%         data_value.peak_area_Lac(:,i)=trapz(details_para.ppm_referenced(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2)),real(data_value.fit_FD_Lac(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2),i)));
        data_value.corr_baseline(:,i)=real(data_value.FID_FD(:,i)) - real(fftshift(fft(baseline_t(:,count))));
%         data_value.peak_area_NAA(:,i)=max(data_value.corr_baseline(:,i));
        data_value.noise(:,i)=std(data_value.corr_baseline(end-50:end,i));%std(real(data_value.FID_FD(:,i))-real(fftshift(fft(fit_sig))));
%         data_value.noise(:,i)=std(real(data_value.FID_FD(end-30:end,i)));
        data_value.snr(:,i)=data_value.peak_area_NAA(:,i)./data_value.noise(:,i);
        data_value.NAA_Cr(:,i)=data_value.peak_area_NAA(:,i)./data_value.peak_area_Cr(:,i);
        data_value.fit_temp(:,i)=data_value.fit_TD(:,i);
        data_value.Cho_Cr(:,i)=data_value.peak_area_Cho(:,i)./data_value.peak_area_Cr(:,i);
        data_value.Myo_Cr(:,i)=data_value.peak_area_Myo(:,i)./data_value.peak_area_Cr(:,i);
        data_value.Lac_Cr(:,i)=data_value.peak_area_Lac(:,i)./data_value.peak_area_Cr(:,i);
        data_value.fit_TD(:,i)=1*data_value.fit_TD(:,i) + 1*baseline_t(:,count);%baseline_t(:,count-1);%data_value.fit_TD(:,i) + ifft(ifftshift(data_value.FID_FD_base(:,i)));
        data_value.fit_FD(:,i)=fftshift(fft(data_value.fit_TD(:,i)));
        data_value.fit_baseline(:,i)= data_value.fit_FD(:,i)-fftshift(fft(data_value.fit_temp(:,i)));
        waitbar(i/details_para.num_vox);
    end
elseif (details_para.fit_method==1)
    for i = 1:details_para.num_vox
        details_para.fla=0;
        details_para.seg = [1 length(data_value.FID_TD(:,i))];
%         if i==1
%              [ fit_sig,amp,Ind_sig,freq,Z_sig,FWHM] = HSVDPK_new(data_value.FID_TD(:,i),details_para.fres,details_para.Fs,fit_para.order(cur_vox));
            [No_pks,pk_loc,FWHM,signal_power_each] = find_peaks_new(data_value.FID_TD(:,i),details_para.fres);
%         end
        freq=details_para.fres(pk_loc);
        FWHM=FWHM/1;
        [fit_sig,~,amp_est] = model_varpro(data_value.FID_TD(:,i),details_para.fres(pk_loc), FWHM, ones(1,length(pk_loc)),details_para.t);
        amp_est=amp_est;
        Initx = [amp_est' (FWHM/1)' freq'];
        lb= [(amp_est- data_value.FID_TD(1,i))' ((FWHM/1)- 10)' (freq-10)'];
        ub= [(amp_est + data_value.FID_TD(1,i))' ((FWHM/1) + 10)' (freq+10)'];
        options = optimset('MaxFunEvals', 5000, 'MaxIter', 5000, 'Algorithm','trust-region-reflective', 'Display', 'off'); % optimization parameters 
        [fit_param1,~,residual,~,~,~,J] = lsqnonlin(@(x)model_estimate_equation2(x,details_para.t,(data_value.FID_TD(:,i)),1),Initx,lb,ub,options); % non-linear least squares algorithm
        [sig_sum,C,Z_sig] = model_estimate(fit_param1,details_para.t,1);
        Ind_sig=C;
        ind=find((freq>=details_para.fit_freq(1)-20 & freq<=details_para.fit_freq(2)+20));
        ind_NAA=find((freq>=details_para.fit_freq_NAA(1)-20 & freq<=details_para.fit_freq_NAA(2)+20));
        ind_Cr=find((freq>=details_para.fit_freq_Cr(1)-20 & freq<=details_para.fit_freq_Cr(2)+20));
        ind_Cho=find((freq>=details_para.fit_freq_Cho(1)-20 & freq<=details_para.fit_freq_Cho(2)+20));
        half_sig2=sum(Ind_sig(:,ind),2);
        data_value.fit_TD(:,i)=half_sig2;
        data_value.fit_FD_NAA(:,i)=fftshift(fft(sum(Ind_sig(:,ind_NAA),2)));
        data_value.fit_FD_Cr(:,i)=fftshift(fft(sum(Ind_sig(:,ind_Cr),2)));
        data_value.fit_FD_Cho(:,i)=fftshift(fft(sum(Ind_sig(:,ind_Cho),2)));
        data_value.fit_FD(:,i)=fftshift(fft(data_value.fit_TD(:,i)));
        data_value.peak_area(:,i)=trapz(details_para.ppm_referenced,abs(data_value.fit_FD(:,i)));
        data_value.peak_area_NAA(:,i)=trapz(details_para.ppm_referenced,abs(data_value.fit_FD_NAA(:,i)));
        data_value.peak_area_Cr(:,i)=trapz(details_para.ppm_referenced,abs(data_value.fit_FD_Cr(:,i)));
        data_value.peak_area_Cho(:,i)=trapz(details_para.ppm_referenced,abs(data_value.fit_FD_Cho(:,i)));
        data_value.noise(:,i)=std(data_value.FID_FD(:,i)-fftshift(fft(fit_sig)));
        data_value.snr(:,i)=data_value.peak_area_NAA(:,i)./data_value.noise(:,i);
        data_value.NAA_Cr(:,i)=data_value.peak_area_NAA(:,i)./data_value.peak_area_Cr(:,i);
        data_value.Cho_Cr(:,i)=data_value.peak_area_Cho(:,i)./data_value.peak_area_Cr(:,i);
        waitbar(i/details_para.num_vox);
    end
else
    for i = 1:details_para.num_vox
        details_para.fla=0;
        data_value.FID_FD_base(:,i)=ssa(data_value.FID_FD(:,i),100);
        data_value.FID_FD_b(:,i)=data_value.FID_FD(:,i)-data_value.FID_FD_base(:,i);
        data_value.FID_TD_b(:,i)=ifft(ifftshift(data_value.FID_FD_b(:,i)));
        details_para.seg = [1 length(data_value.FID_TD(:,i))];
%         if i==1
            [No_pks,pk_loc,FWHM(1,:),signal_power(1,:)] = find_peaks_new(data_value.FID_TD_b(:,i),details_para.fres);
%         end
%         [ fit_sig,amp,Ind_sig,freq,Z_sig,FWHM] = HSVDPK_new(data_value.FID_TD(:,i),details_para.fres,details_para.Fs,fit_para.order(cur_vox));
        freq=details_para.fres(pk_loc);
        FWHM=1*FWHM*(1E6/details_para.Tf)./2;
        FWHM(1:end)=FWHM(1);
        for kk=1:length(FWHM)%length(pk_loc)
            Data2bFit=(fftshift(fft(data_value.FID_TD_b(:,i))))/1;
            PEAK_ppm=4.7 + freq(kk)*(1E6/details_para.Tf);%details_para.ppm_referenced(pk_loc(kk));
            z=abs(details_para.ppm_referenced-(PEAK_ppm+FWHM(kk)));
            ub=find(min(z)==z);
            z=abs(details_para.ppm_referenced-(PEAK_ppm-FWHM(kk)));
            lb=find(min(z)==z); 
            freqrange = details_para.ppm_referenced(lb:ub);
            SumSpec_half = (Data2bFit(lb:ub,:));
            max_peak=max(SumSpec_half);
            Init = [max_peak./2 FWHM(kk) PEAK_ppm 0 0 0 ];
            dataParams_tot(kk,:) = FitPeaks(freqrange, SumSpec_half, Init);
            dataParams_tot(kk,4:6)=0;
%             dataParams_tot(kk,2)=0;
            fitted_signal_tot(kk,:)=LorentzModel(dataParams_tot(kk,:), details_para.ppm_referenced);
            Ind_sig(:,kk)=fitted_signal_tot(kk,:);
            Ind_sig(:,kk)=ifft(ifftshift(Ind_sig(:,kk)));
        end
        fit_sig=sum(Ind_sig,2);
%         fit_sig=ifft(ifftshift(fit_sig));
        ind=find((freq>=details_para.fit_freq(1)-20 & freq<=details_para.fit_freq(2)+20));
        ind_NAA=find((freq>=details_para.fit_freq_NAA(1)-20 & freq<=details_para.fit_freq_NAA(2)+20));
        ind_Cr=find((freq>=details_para.fit_freq_Cr(1)-20 & freq<=details_para.fit_freq_Cr(2)+20));
        ind_Cho=find((freq>=details_para.fit_freq_Cho(1)-20 & freq<=details_para.fit_freq_Cho(2)+20));
        half_sig2=sum(Ind_sig(:,ind),2);
        data_value.fit_TD(:,i)=half_sig2;
        data_value.fit_FD_NAA(:,i)=fftshift(fft(sum(Ind_sig(:,ind_NAA),2)));
        data_value.fit_FD_Cr(:,i)=fftshift(fft(sum(Ind_sig(:,ind_Cr),2)));
        data_value.fit_FD_Cho(:,i)=fftshift(fft(sum(Ind_sig(:,ind_Cho),2)));
        data_value.fit_FD(:,i)=fftshift(fft(data_value.fit_TD(:,i)));
        data_value.peak_area(:,i)=trapz(details_para.ppm_referenced,abs(data_value.fit_FD(:,i)));
        data_value.peak_area_NAA(:,i)=trapz(details_para.ppm_referenced,abs(data_value.fit_FD_NAA(:,i)));
        data_value.peak_area_Cr(:,i)=trapz(details_para.ppm_referenced,abs(data_value.fit_FD_Cr(:,i)));
        data_value.peak_area_Cho(:,i)=trapz(details_para.ppm_referenced,abs(data_value.fit_FD_Cho(:,i)));
        data_value.noise(:,i)=std(data_value.FID_FD(:,i)-fftshift(fft(fit_sig)));
        data_value.snr(:,i)=data_value.peak_area_NAA(:,i)./data_value.noise(:,i);
        data_value.NAA_Cr(:,i)=data_value.peak_area_NAA(:,i)./data_value.peak_area_Cr(:,i);
        data_value.Cho_Cr(:,i)=data_value.peak_area_Cho(:,i)./data_value.peak_area_Cr(:,i);
        data_value.fit_TD(:,i)=data_value.fit_TD(:,i) + ifft(ifftshift(data_value.FID_FD_base(:,i)));
        data_value.fit_FD(:,i)=fftshift(fft(data_value.fit_TD(:,i)));
        waitbar(i/details_para.num_vox);
    end
end
close(h);
cur_data(:,1)=data_value.peak_area;
details_para.mean_peak_area=mean(cur_data(:,1));
details_para.std_peak_area=std(cur_data(:,1));
cur_data(:,2)=data_value.snr;
cur_data(:,3)=data_value.NAA_Cr;
cur_data(:,4)=data_value.Cho_Cr;
cur_data(:,5)=data_value.Myo_Cr;
cur_data(:,6)=data_value.Lac_Cr;
FID_FD_copy=data_value.FID_FD;
details_para.fit_flag=0;
x_max = 850;
y_max = 600;
col_name = {'Area','SNR','NAA/Cr','Cho/Cr','Myo/Cr','Lac/Cr'};
row_name = cell(1,details_para.num_vox);
for i=1:details_para.num_vox
    row_name{i} = strcat('voxel ', num2str(i));
end
diplay_handles.para_table = uitable(fitting_handles.fig,'Units','normalized','position',[550/x_max,100/y_max,260/x_max,220/y_max], 'data',cur_data,...
    'RowName', row_name, 'ColumnName', col_name,'ColumnWidth',{85});
% seg = details_para.disp_seg;
% ppm = details_para.ppm_referenced;
% cla(fitting_handles.axis);
% plot_var = real(data_value.FID_FD(seg(1):seg(2)));
% plot(fitting_handles.axis,ppm(seg(1):seg(2))',plot_var);
% max_p = max(plot_var);
% mim_p = min(plot_var);
% axis(fitting_handles.axis,[ppm(seg(1)), ppm(seg(2)),mim_p - 0.1*max_p, max_p + 0.1*max_p])
% set(fitting_handles.axis,'XDir','reverse')
details_para.fit_flag=1;
display_data()
set(fitting_handles.map , 'Enable', 'on');
% set(fitting_handles.display_grid , 'Enable', 'on');
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
%------------------------------------------------------------------------------------------------------------------------------------------
function map_Callback(hObject, eventdata)
global details_para;
global data_value;
global cur_vox;
global fitting_handles;

set(fitting_handles.fig,'Visible','off');
set(fitting_handles.display_grid, 'Enable', 'off');
set(fitting_handles.map, 'Enable', 'off');
set(fitting_handles.close, 'Enable', 'off');
PE=details_para.PE;
l=1;
for i=1:details_para.PE/1
    for j=1:details_para.RE
        details_para.map_csi(i,j)=data_value.peak_area(:,l);
        l=l+1;
    end
end
data_value.H0=figure();
set(data_value.H0,'Name','RAPID 1.0','NumberTitle', 'off');
imagesc(details_para.map_csi);
colorbar;
try
[x,y]=getpts();
% if (~isempty(x))
    row_no=round(x); 
    col_no= round(y);
    close;
    PE_len=1;
    RE_len=PE_len;
    for lll=1:PE_len
       grid(lll,:)= (PE*(lll-1)+ PE*(col_no-1) + row_no : PE*(lll-1)+PE*(col_no-1) +row_no+ RE_len -1)' ;
    end
    grid_all=unique(grid);
    selected_vox = details_para.selected_voxels;
    cur_vox = selected_vox(1,grid_all);
    set(fitting_handles.voxel_list,'value',cur_vox);
    display_data();
    set(fitting_handles.map, 'Enable', 'on');
    set(fitting_handles.close, 'Enable', 'on');
catch
    set(fitting_handles.fig,'Visible','on');
    set(fitting_handles.close, 'Enable', 'on');
    set(fitting_handles.map, 'Enable', 'on');
    return;
end
set(fitting_handles.display_grid, 'Enable', 'on');
set(fitting_handles.fig,'Visible','on');
%------------------------------------------------------------------------------------------------------------------------------------------
function display_grid_Callback(hObject, eventdata)
global details_para;
global data_value;
global fitting_display;
global fitting_handles;


set(fitting_handles.fig,'Visible','off');
set(fitting_handles.map , 'Enable', 'off');
set(fitting_handles.display_grid, 'Enable', 'off');
set(fitting_handles.close, 'Enable', 'off')

seg = details_para.disp_seg;
ppm = details_para.ppm_referenced;
PE=details_para.PE;
cnt=0;
try
while(1)
    x = inputdlg(['Enter grid size (for eg. enter 3 for 3x3 size): ']);
    if((str2double(x)<=0) || (str2double(x)>details_para.PE) || (isnan(str2double(x))))
        cnt=cnt+1;
        h1(cnt)=errordlg('Invalid input');
    else
        PE_len=str2double(x);
        for i=1:(cnt) 
            close(h1(i));
        end
        break;
    end
end
data_value.H1=figure();
set(data_value.H1,'Name','RAPID 1.0','NumberTitle', 'off');
imagesc(details_para.map_csi);
% set(fitting_handles.fig,'Visible','off');
[x,y]=getpts();
row_no=round(x); 
col_no= round(y);
close;
RE_len=PE_len;
for lll=1:PE_len
   grid(lll,:)= (PE*(lll-1)+ PE*(col_no-1) + row_no : PE*(lll-1)+PE*(col_no-1) +row_no+ RE_len -1)' ;
end
grid_all=unique(grid);
% set(fitting_handles.fig,'Visible','on');
data_value.H2=figure();
set(data_value.H2,'Name','RAPID 1.0','NumberTitle', 'off');
count=0;
for lll=1:PE_len
    for i=grid_all(1):grid_all(RE_len)%1:(PE/1);%1:(PE/1)*RE;%1:(PE/1)*RE;%1:(PE/1);%1:(PE/1)*RE ;%1:(PE/1)*RE %1:16
        subplot(PE_len,RE_len,count+1);%subplot(1,RE,i);%subplot(PE/4,RE,i);%subplot(1,RE,i);%subplot(PE/1,RE,i);%subplot(1,RE,i);%subplot(PE/1,RE,i);%subplot(1,RE,i);%subplot(PE/1,RE,i);%subplot(1,RE,i);%subplot(PE/1,RE,i);%subplot(1,RE,i);%subplot(PE/1,RE,i);
        if (fitting_display == 1)
            temp1=(real(data_value.FID_FD(seg(1):seg(2),i+PE*(lll-1))));%*rat;
        else
             temp1=(abs(data_value.FID_FD(seg(1):seg(2),i+PE*(lll-1))));%*rat;
        end
        plot(ppm(seg(1):seg(2))',temp1);
        max_p = max(temp1);
        min_p = min(temp1);
        axis([ppm(seg(1)), ppm(seg(2)),0, 15])
%         axis([ppm(seg(1)), ppm(seg(2)),min_p - 0.1*max_p, max_p + 0.1*max_p])
        ax=gca;
        li=get(ax,'children');%spessore
        set(li,'linewidth',1);
        axis off
    %     end
    %     legend('CSI data','RAPID data');
        set(gca,'XDir','reverse');
        if i==1
%         title('abs part'); 
        end
        count=count+1;
    end
end
waitfor(data_value.H2);
set(fitting_handles.display_grid, 'Enable', 'on');
set(fitting_handles.close, 'Enable', 'on');
catch
    set(fitting_handles.fig,'Visible','on');
    set(fitting_handles.close, 'Enable', 'on');
    set(fitting_handles.display_grid, 'Enable', 'on');
    return;
end
set(fitting_handles.map , 'Enable', 'on');
set(fitting_handles.fig,'Visible','on');
%------------------------------------------------------------------------------------------------------------------------------------------