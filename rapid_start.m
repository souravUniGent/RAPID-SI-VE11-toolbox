
function rapid_start()

% This function is used for starting the RAPID GUI

% Author: Sourav Bhaduri
% University of Ghent
% email: Sourav.Bhaduri@UGent.be
% Feb 2018; Last revision: Feb 2018

clear all
clear global all
close all
clc

global mainGUI_handles;
global details_para;
global globalpara;
global data_value;

globalpara=[];
x_max = 300;
y_max = 390;

data_value.scsa_flag=0;
data_value.sim_flag=1;
details_para.slice_no=1;
details_para.image_flag=0;
data_value.save_flag=0;

% main figure creation
mainGUI_handles.fig = figure('Visible','off','Position',[400,400,x_max,y_max],'Color',[0.68,0.92,1],'Name','RAPID 1.0','NumberTitle', 'off');
set(mainGUI_handles.fig,'MenuBar','none'); 

mainGUI_handles.load_image = uicontrol(mainGUI_handles.fig,'Style', 'pushbutton','Units','normalized','Position', [.04 .85 .4 .07],...
'String', 'Load image','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@load_image_Callback);
mainGUI_handles.load_wo_image = uicontrol(mainGUI_handles.fig,'Style', 'pushbutton','Units','normalized','Position', [.5 .85 .4 .07],...
'String', 'Load w/o image','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@load_wo_image_Callback);


% menu set-up
% set(mainGUI_handles.fig,'MenuBar','none'); 
% mainGUI_handles.project = uimenu(mainGUI_handles.fig,'Label','Project');
% mainGUI_handles.new_project = uimenu(mainGUI_handles.project,'Label','New Project','Callback',@new_project_callback);
% mainGUI_handles.load_project = uimenu(mainGUI_handles.project,'Label','Load Project','Callback',@load_project_callback);
% mainGUI_handles.save_project = uimenu(mainGUI_handles.project,'Label','Save Project','Callback',@save_project_callback);
% mainGUI_handles.save_figure = uimenu(mainGUI_handles.project,'Label','Save Figure','Callback',@save_main_fig_callback);
% 
% mainGUI_handles.load = uimenu(mainGUI_handles.fig,'Label','Load');
% mainGUI_handles.load_data = uimenu(mainGUI_handles.load,'Label','MRS Data','Callback',@load_data_callback);
% mainGUI_handles.load_mri = uimenu(mainGUI_handles.load,'Label','MRI Image','Callback',@MRI_load_callback);
% 
% mainGUI_handles.settings = uimenu(mainGUI_handles.fig,'Label','Settings');
% mainGUI_handles.save_setting = uimenu(mainGUI_handles.settings,'Label','Create and Save Settings','Callback',@save_settings_callback);
% mainGUI_handles.load_setting = uimenu(mainGUI_handles.settings,'Label','Load Settings');%,'Callback',@load_setting_callback);
% mainGUI_handles.load_default_setting = uimenu(mainGUI_handles.load_setting,'Label','Default','Callback',@load_default_setting_callback);
% mainGUI_handles.load_user_setting = uimenu(mainGUI_handles.load_setting,'Label','User Defined','Callback',@load_user_setting_callback);
% mainGUI_handles.Method = uimenu(mainGUI_handles.fig,'Label','Method');
% mainGUI_handles.automatic = uimenu(mainGUI_handles.Method ,'Label','Automatic','Callback',@auto_callback);
% mainGUI_handles.manual= uimenu(mainGUI_handles.Method ,'Label','Manual','Callback',@manual_callback);
% 
% mainGUI_handles.mathematics = uimenu(mainGUI_handles.fig,'Label','Mathematics');
% mainGUI_handles.save_setting = uimenu(mainGUI_handles.mathematics,'Label','Sum','Callback',@sum_callback);
% mainGUI_handles.load_setting = uimenu(mainGUI_handles.mathematics,'Label','Average','Callback',@average_callback);
% mainGUI_handles.load_default_setting = uimenu(mainGUI_handles.mathematics,'Label','Subtract','Callback',@subtract_callback);
% 
% mainGUI_handles.segmentation = uimenu(mainGUI_handles.fig,'Label','Segmentation');
% mainGUI_handles.segment = uimenu(mainGUI_handles.segmentation,'Label','Segment','Callback',@segment_callback);
% mainGUI_handles.segment_display = uimenu(mainGUI_handles.segmentation,'Label','Display','Callback',@segment_display_callback);
% mainGUI_handles.segment_bias = uimenu(mainGUI_handles.segmentation,'Label','Bias correct','Callback',@segment_bias_callback);


%axis and voxel list setup
mainGUI_handles.slice_selection_panel = uipanel('Title','','Units','normalized','Position',[.1 .02 .75 .8],'BackgroundColor',[0.68,0.92,1],...
    'FontSize',9,'FontWeight','Bold','TitlePosition','centertop');
% mainGUI_handles.slice_selection_list = uicontrol('Parent',mainGUI_handles.slice_selection_panel,'Style','listbox','Units','normalized','Position',...
%     [.05 .027 .9 .97],'Callback', @slice_selection_callback,'Interruptible','off','BusyAction','cancel');
mainGUI_handles.slice_load = uicontrol(mainGUI_handles.slice_selection_panel,'Style', 'pushbutton','Units','normalized','Position', [.09 .8 .8 .2],...
'String', 'Start','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@slice_load_Callback);
mainGUI_handles.slice_preprocess = uicontrol(mainGUI_handles.slice_selection_panel,'Style', 'pushbutton','Units','normalized','Position', [.09 .55 .8 .2],...
'String', 'Preprocess','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@slice_preprocess_Callback);
mainGUI_handles.slice_fit = uicontrol(mainGUI_handles.slice_selection_panel,'Style', 'pushbutton','Units','normalized','Position', [.09 .3 .8 .2],...
'String', 'Fit','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@slice_fit_Callback);
mainGUI_handles.slice_map = uicontrol(mainGUI_handles.slice_selection_panel,'Style', 'pushbutton','Units','normalized','Position', [.09 .05 .8 .2],...
'String', 'Metabolite map','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@slice_map_Callback);

mainGUI_handles.close = uicontrol(mainGUI_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[255/x_max 1/y_max 45/x_max 25/y_max],...
'String', 'Close','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@close_Callback);
mainGUI_handles.Help = uicontrol(mainGUI_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[5/x_max 370/y_max 45/x_max 18/y_max],...
'String', 'Help','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@help_Callback);
mainGUI_handles.save = uicontrol(mainGUI_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[255/x_max 35/y_max 45/x_max 25/y_max],...
'String', 'Save','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@save_Callback);

set(mainGUI_handles.fig,'Visible','on');
set(mainGUI_handles.slice_load, 'Enable', 'off');
set(mainGUI_handles.slice_preprocess, 'Enable', 'off');
set(mainGUI_handles.slice_fit, 'Enable', 'off');
set(mainGUI_handles.slice_map, 'Enable', 'off');
set(mainGUI_handles.load_image, 'Enable', 'on');
set(mainGUI_handles.save, 'Enable', 'off');

waitfor(mainGUI_handles.fig);

%------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function slice_load_Callback(hObject, eventdata)

global details_para;
global mainGUI_handles;
global data_value;
global preprocess_handles;
global preprocess_sel_vox;
global FID_FD_copy;
global met_tick;
global fitting_handles;
global preprocess_display;
global fitting_display;
global cur_vox;
global ind_allvox;
global zerofill_value;
global default_parameters;
global pha_cor_data;
global fit_para;
global apod_para;
global pha_para;

% set(mainGUI_handles.slice_preprocess, 'Enable', 'on');
set(mainGUI_handles.slice_load, 'Enable', 'off');
set(mainGUI_handles.slice_fit, 'Enable', 'off');
set(mainGUI_handles.slice_map, 'Enable', 'off');

try
    if (data_value.sense_flag==0)
        rapid_main_gui();
    else
        rapid_sense_gui();
    end
catch
    set(mainGUI_handles.slice_preprocess, 'Enable', 'off');
    return;
end
try
if (data_value.data_error==0)
    set(mainGUI_handles.slice_preprocess, 'Enable', 'on');
else
    set(mainGUI_handles.slice_preprocess, 'Enable', 'off');
    set(mainGUI_handles.slice_fit, 'Enable', 'off');
    set(mainGUI_handles.slice_map, 'Enable', 'off');
end
if (data_value.file_flag==1)
    if (data_value.scsa_flag==1)
        load(data_value.temp_filename,'details_para','data_value');
        set(mainGUI_handles.slice_preprocess, 'Enable', 'on');
    elseif (data_value.sim_flag==1)
         rapid_sim_gui();
         set(mainGUI_handles.slice_preprocess, 'Enable', 'on');
    else
        load(data_value.temp_filename,'details_para','data_value','preprocess_handles','FID_FD_copy','met_tick','preprocess_sel_vox','fitting_handles','apod_para','fit_para','ind_allvox','pha_cor_data','default_parameters','zerofill_value','cur_vox','fitting_display','preprocess_display','pha_para');
        set(mainGUI_handles.slice_preprocess, 'Enable', 'on');
    end
end
set(mainGUI_handles.slice_load, 'Enable', 'on');
catch
    return;
end

%------------------------------------------------------------------------------------------------------------------------------------------------------

function slice_preprocess_Callback(hObject, eventdata)

global details_para;
global mainGUI_handles;
global data_value;
global preprocess_handles;

set(mainGUI_handles.fig,'Visible','off');
set(mainGUI_handles.slice_fit, 'Enable', 'off');
set(mainGUI_handles.slice_map, 'Enable', 'off');
set(mainGUI_handles.slice_preprocess, 'Enable', 'off');

if (isempty(preprocess_handles))
    preprocessing_new();
else
    set(mainGUI_handles.slice_preprocess, 'Enable', 'off');
    preprocessing_new();
end
set(mainGUI_handles.slice_preprocess, 'Enable', 'on');   
set(mainGUI_handles.slice_fit, 'Enable', 'on');
set(mainGUI_handles.fig,'Visible','on');
% set(mainGUI_handles.load_image, 'Enable', 'on');

%------------------------------------------------------------------------------------------------------------------------------------------------------

function slice_fit_Callback(hObject, eventdata)

global details_para;
global mainGUI_handles;
global data_value;
global fitting_handles;

set(mainGUI_handles.fig,'Visible','off');
set(mainGUI_handles.slice_fit, 'Enable', 'off');
set(mainGUI_handles.slice_preprocess, 'Enable', 'off');
set(mainGUI_handles.slice_map, 'Enable', 'off');

if (isempty(fitting_handles))
    fitting_new();
else
     set(mainGUI_handles.slice_fit, 'Enable', 'off');
     fitting_new();
end
try
    set(mainGUI_handles.slice_fit, 'Enable', 'on');
catch
    return;
end
set(mainGUI_handles.slice_preprocess, 'Enable', 'on');
set(mainGUI_handles.fig,'Visible','on');
set(mainGUI_handles.slice_map, 'Enable', 'on');
set(mainGUI_handles.save, 'Enable', 'on');
%------------------------------------------------------------------------------------------------------------------------------------------------------

function slice_map_Callback(hObject, eventdata)

global details_para;
global mainGUI_handles;
global data_value;
global preprocess_handles;
global preprocess_sel_vox;
global FID_FD_copy;
global met_tick;

set(mainGUI_handles.close, 'Enable', 'off');
set(mainGUI_handles.slice_map, 'Enable', 'off');

set(mainGUI_handles.fig,'Visible','off');

% for i = 1:details_para.num_vox
%     [ fit_sig,amp,Ind_sig,freq] = HSVDPK_new(data_value.FID_TD(:,i),details_para.fres,details_para.Fs,10);
%     ind=find((freq>=-80 & freq<=80));
%     half_sig2=sum(Ind_sig(:,ind),2);
%     data_value.FID_TD(:,i)=data_value.FID_TD(:,i)- half_sig2;
%     data_value.FID_FD(:,i)=fftshift(fft(data_value.FID_TD(:,i)));
% end
PE=details_para.PE;
RE=details_para.RE;
seg_NAA=details_para.fit_disp_seg_NAA(1):details_para.fit_disp_seg_NAA(2);
seg_Lac=details_para.fit_disp_seg_Lac(1):details_para.fit_disp_seg_Lac(2);%210:220;
seg_Cho=details_para.fit_disp_seg_Cho(1):details_para.fit_disp_seg_Cho(2);%205:210;%200:207;%205:210;%205:210;%205:210;%200:205;198:210;
seg_Cr=details_para.fit_disp_seg_Cr(1):details_para.fit_disp_seg_Cr(2);%200:205;%200:205;198:210;
seg_mi=details_para.fit_disp_seg_mi(1):details_para.fit_disp_seg_mi(2);%200:205;%200:205;198:210;
% if round(details_para.Fs)==2000
%     seg_NAA=160:180;
%     seg_Lac=130:159;%210:220;
%     seg_Cho=209:218;%205:210;%200:207;%205:210;%205:210;%205:210;%200:205;198:210;
%     seg_Cr=196:205;%200:205;%200:205;198:210;
% elseif round(details_para.Fs)==1600
%     seg_NAA=135:165;
%     seg_Lac=110:130;
%     seg_Cho=185:200;
% elseif round(details_para.Fs)==1000
%     seg_NAA=60:95;
%     seg_Lac=30:50;
%     seg_Cho=145:165;
%     seg_Cr=145:165;
% elseif round(details_para.Fs)==2500
%     seg_NAA=186:200;
%     seg_Lac=165:185;
%     seg_Cho=210:225;
% end
for i=1:(PE/1)*RE
    temp1=abs(data_value.FID_FD(:,i));
%     if i==1
%     details_para.fla=0;
%     details_para.seg = [1 length(data_value.FID_TD(:,1))];
%     [No_pks,pk_loc,FWHM,signal_power] = find_peaks_new(data_value.FID_TD(:,424),details_para.fres);
%     details_para.No_pks=No_pks;
%     end
% %     else
%         details_para.fla=1;
%         [No_pks,pk_loc,FWHM,signal_power] = find_peaks_new(data_value.FID_TD(:,i),details_para.fres);
% %     end
    NAA_CSI(i)=(max(temp1(seg_NAA))/1) - (min(temp1(seg_NAA))/1);
    Lac_CSI(i)=(max(temp1(seg_Lac))/1) - ((min(temp1(seg_Lac))/1));
    Cho_CSI(i)=(max(temp1(seg_Cho))/1) - (min(temp1(seg_Cho))/1);
    Cr_CSI(i)=(max(temp1(seg_Cr))) - (min(temp1(seg_Cr)));
    MI_CSI(i)=(max(temp1(seg_mi))) - ((min(temp1(seg_mi))));

end
l=1;
for i=1:PE/1
    for j=1:RE
        NAA_map_csi(i,j)=NAA_CSI(l).*1;%5;
        Lac_map_csi(i,j)=Lac_CSI(l).*1;%5;
        Cho_map_csi(i,j)=Cho_CSI(l).*1;%5;
        Cr_map_csi(i,j)=Cr_CSI(l).*1;%5;
        MI_map_csi(i,j)=MI_CSI(l).*1;%5;
        l=l+1;
    end
end
met_tick=1;
preprocess_handles.met_fig=figure();
set(preprocess_handles.met_fig,'Name','RAPID 1.0','NumberTitle', 'off');
subplot(2,3,1)
imagesc(Cr_map_csi);%./Cr_map_csi);
title('Cr map');
colormap('hot')
colorbar
subplot(2,3,2)
imagesc(NAA_map_csi);%./Cr_map_csi);
title('NAA map');
colormap('hot')
colorbar
subplot(2,3,3)
imagesc(Cho_map_csi);%./Cr_map_csi);
title('Cho map');
colormap('hot')
colorbar
subplot(2,3,4)
imagesc(Lac_map_csi);%./Cr_map_csi);
title('Lac map');
colormap('hot')
colorbar
subplot(2,3,5)
imagesc(MI_map_csi);%./Cr_map_csi);
title('MI map');
colormap('hot')
colorbar
waitfor(preprocess_handles.met_fig);
set(mainGUI_handles.fig,'Visible','on');
set(mainGUI_handles.slice_map, 'Enable', 'on');
set(mainGUI_handles.close, 'Enable', 'on');


%------------------------------------------------------------------------------------------------------------------------------------------------------
function load_image_Callback(hObject, eventdata)

global details_para;
global mainGUI_handles;
global data_value;
global preprocess_handles;
global preprocess_sel_vox;
global FID_FD_copy;
global met_tick;

set(mainGUI_handles.slice_load, 'Enable', 'on');

try
load_image_new();
if (data_value.image_error==0)
    set(mainGUI_handles.slice_load, 'Enable', 'on');
else
    set(mainGUI_handles.slice_load, 'Enable', 'off');
end
details_para.image_flag=1;
catch
    return;
end





%--------------------------------------------------------------------------------------------------------------------------------------------------------
function load_wo_image_Callback(hObject, eventdata)

global details_para;
global mainGUI_handles;
global data_value;
global preprocess_handles;
global preprocess_sel_vox;
global FID_FD_copy;
global met_tick;

data_value.sense_flag=0;
set(mainGUI_handles.slice_load, 'Enable', 'on');


function close_Callback(hObject, eventdata)
%Callback for 'close' button

try
    close all;
catch 
    return;
end

function help_Callback(hObject, eventdata)
%Callback for 'close' button

file=strcat(pwd,'\HELP.docx');

try
    open(file);
catch 
    return;
end

function save_Callback(hObject, eventdata)
%Callback for 'close' button
global details_para;
global mainGUI_handles;
global data_value;
global preprocess_handles;
global preprocess_sel_vox;
global FID_FD_copy;
global met_tick;
global fitting_handles;
global preprocess_display;
global fitting_display;
global cur_vox;
global ind_allvox;
global zerofill_value;
global default_parameters;
global pha_cor_data;
global fit_para;
global apod_para;
global pha_para;

set(mainGUI_handles.save, 'Enable', 'off');

file_exten_ind=strfind(data_value.spect_name,'.');
file_exten=data_value.spect_name(file_exten_ind-4:1:file_exten_ind-1);
file=strcat(pwd,'\',file_exten,data_value.slice_profile,'.mat');
data_value.save_flag=1;
h = waitbar(0,'Please wait...');
try
    save(file,'details_para','mainGUI_handles','data_value','preprocess_handles','FID_FD_copy','met_tick','preprocess_sel_vox','fitting_handles','apod_para','fit_para','ind_allvox','pha_cor_data','default_parameters','zerofill_value','cur_vox','fitting_display','preprocess_display','pha_para');
catch 
    set(mainGUI_handles.save, 'Enable', 'on');
    return;
end
waitbar(1);
close(h);
set(mainGUI_handles.save, 'Enable', 'on');




%--------------------------------------------------------------------------------------------------------------------------------------------------------