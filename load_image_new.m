function[]= load_image_new()

% This function is used for loading the MRI localizer image and initiate
% parallel imaging 

% Author: Sourav Bhaduri
% University of Ghent
% email: Sourav.Bhaduri@UGent.be
% Feb 2018; Last revision: Feb 2018

global data_value;
global image_handles;
global details_para;
global apod_para;
global pha_para;
global ind_allvox;
global pha_cor_data;
global default_parameters;
global zerofill_value;
global info_struct;
global cur_vox;
global met_tick;

data_value.image_name=cell(0);
data_value.sense_flag=0;
data_value.image_error=0;
met_tick=0;
try
    twix_obj = mapVBVD();
    close(data_value.h);
catch
    close(data_value.h);
    return;
end
filename=twix_obj.image.filename;
data_value.image_name{1}=filename(length(pwd)+1:end);
data=twix_obj.image;
global_rescale=1e5;
pointsBeforeEcho=data.freeParam(1);
fid=fopen(filename);
line=fgets(fid);
index=findstr(line,'sProtConsistencyInfo.flNominalB0');
equals_index=findstr(line,'= ');
while isempty(index) || isempty(equals_index)
    line=fgets(fid);
    index=findstr(line,'sProtConsistencyInfo.flNominalB0');
    equals_index=findstr(line,'= ');
end
Bo=line(equals_index+1:end);
Bo=str2double(Bo);
fclose(fid);
has_header = isfield(twix_obj,'hdr');
if has_header
    dwelltime = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1,1} * 1e-9;
    spectralwidth=1/dwelltime;
    PE=twix_obj.hdr.Dicom.lPhaseEncodingLines;
    RE=twix_obj.hdr.Dicom.lBaseResolution;
else
    fid=fopen(filename);
    line=fgets(fid);
    index=findstr(line,'sRXSPEC.alDwellTime[0]');
    equals_index=findstr(line,'= ');
    while isempty(index) || isempty(equals_index)
        line=fgets(fid);
        index=findstr(line,'sRXSPEC.alDwellTime[0]');
        equals_index=findstr(line,'= ');
    end
    dwelltime=line(equals_index+1:end);
    dwelltime=str2double(dwelltime)*1e-9;
    spectralwidth=1/dwelltime;
    fclose(fid);
end
fid=fopen(filename);
line=fgets(fid);
index=findstr(line,'sTXSPEC.asNucleusInfo[0].lFrequency');
equals_index=findstr(line,'= ');
while isempty(index) || isempty(equals_index)
    line=fgets(fid);
    index=findstr(line,'sTXSPEC.asNucleusInfo[0].lFrequency');
    equals_index=findstr(line,'= ');
end
txfrq=line(equals_index+1:end);
txfrq=str2double(txfrq);
fclose(fid);
struct.sw=spectralwidth;
struct.LarmorFreq = Bo*42.577;          
struct.nrows = twix_obj.image.NAcq;
rc_xres = double(twix_obj.image.NCol);
rc_yres = double(twix_obj.image.NAcq);
nreceivers = double(twix_obj.image.NCha);
FullData=permute(reshape(double(twix_obj.image()),[twix_obj.image.NCol twix_obj.image.NCha twix_obj.image.NLin*twix_obj.image.NPhs twix_obj.image.NSet twix_obj.image.NSli twix_obj.image.NAve]),[4 2 1 3 5 6]);
FullData=(squeeze(sum(FullData(1,:,:,:,:,:),1)));
FullData=sum(FullData,5)./twix_obj.image.NAve;
FullData = FullData*global_rescale;
data_value.FullData=FullData;
for ll=1:nreceivers
    if twix_obj.image.NSli>1
        data_value.FullData1(:,:,:,ll) = conj(squeeze(sum(FullData(ll,:,:,:),1)));
    else
%         data_value.FullData1(:,:,:,ll) = conj(squeeze(sum(FullData(ll,:,:),1)));
            data_value.image_error=1;
            return;
    end
end
x_max = 600;
y_max = 500;
data_value.FullData = conj(squeeze(sum(data_value.FullData,1)));
image_no = cell(1,size(data_value.FullData,3));
details_para.selected_image=1:size(data_value.FullData,3);
for i = 1:size(data_value.FullData,3)
image_no{i} =  num2str(details_para.selected_image(1,i));
end
cur_vox = details_para.selected_image(1,1);
data_value.final_image=cell(0);
h = waitbar(0,'Please wait...');
for iii=1:size(data_value.FullData,3)
    img_kspace=data_value.FullData(:,:,iii);
    d_kspace = max(max(max(abs(img_kspace))));
    img_image=ifft2c(img_kspace);
    d_image = max(max(abs(img_image)));
    img_kspace1=data_value.FullData1(:,:,iii,:);
    img_kspace1=squeeze(permute((img_kspace1),[3 1 2 4]));
    d_kspace1 = max(max(max(abs(img_kspace))));
    for ll=1:size(img_kspace1,3)
        img_image1(:,:,ll)=ifft2c(img_kspace1(:,:,ll));
    %     img_image2(:,:,ll)=(imresize(img_image1(:,:,ll),[8 8]));
    end
    d_image1 = max(max(max(abs(img_image1))));
    data_value.im_map_mod=(abs(img_image1)/d_image1);
    for i=1:size(data_value.im_map_mod,3)
         temp_img_mod=(data_value.im_map_mod(:,:,i));
         temp_img_mod=imresize(temp_img_mod,([900 600]));
         data_value.im_map_final_mod(:,:,i)=imrotate(temp_img_mod,270*1);
    end
    data_value.final_image{iii}=sum(data_value.im_map_final_mod.^1,3);
    waitbar(iii/size(data_value.FullData,3));
end
close(h);
% Cerate GUI
image_handles.fig = figure('Visible','off','Position',[100,100,x_max,y_max],'Color',[0.68,0.92,1],'Name','RAPID 1.0','NumberTitle', 'off');
set(image_handles.fig,'MenuBar','none'); 
% image_handles.save_figure = uimenu(image_handles.fig,'Label','Save Figure','Callback',@save_preprocess_fig_callback);

image_handles.h1 = uipanel(image_handles.fig,'Title','Slice No','Units','normalized','Position',[15/x_max 190/y_max 70/x_max 270/y_max],...
    'BackgroundColor','white','FontSize',8, 'FontWeight','Bold');
% image_handles.voxel_list = uicontrol(image_handles.h1,'Style', 'listbox','Units','normalized','Position', [.05 .027 .9 .97],'String', vox_no,...
%     'Max', details_para.num_vox, 'Min', 1, 'Callback', @voxel_list_callback, 'Interruptible','off','BusyAction','cancel'); 
image_handles.voxel_list = uicontrol(image_handles.h1,'Style', 'listbox','Units','normalized','Position', [.05 .027 .9 .97],'String',image_no ,...
    'value',1, 'Callback', @voxel_list_callback, 'Interruptible','off','BusyAction','cancel'); 

image_handles.axis = axes('Units','Pixels','Units','normalized','Position',[120/x_max,200/y_max,420/x_max,250/y_max]);
% image_handles.axis_text = uicontrol(image_handles.fig,'Style','text','Units','normalized','Position',[220/x_max 460/y_max 130/x_max 15/y_max],...
%     'String','Signal Spectrum','BackgroundColor',[0.68,0.92,1],'FontSize',9,'FontWeight','Bold');

image_handles.file_panel = uipanel(image_handles.fig,'Title','','Position',[.3 .94 .48 .05],'BackgroundColor','white','FontSize',5,'FontWeight','Bold');
image_handles.main_text = uicontrol(image_handles.file_panel,'Style','text','Units','normalized','Position',[0.005,0.001,0.05,0.1*5],'String',...
    'I','FontSize',6,'FontWeight','Bold','BackgroundColor','white');
image_handles.main_val = uicontrol(image_handles.file_panel,'Style','edit','Units','normalized','Position',[0.05,0.001,0.95,0.15*5],'String',data_value.image_name{1},...
    'Callback',{@file_main_Callback},'Interruptible','off','BusyAction','cancel'); 

% image_handles.close = uicontrol(image_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[480/x_max 8/y_max 100/x_max 25/y_max],...
% 'String', 'Close','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@close_Callback);
% set(image_handles.fig,'Visible','on');

image_handles.sense_map = uicontrol(image_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[480/x_max 40/y_max 100/x_max 45/y_max],...
'String', 'Sensitivity map','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@sense_Callback);
set(image_handles.fig,'Visible','on');

set(image_handles.main_val, 'Enable', 'off')

display_image()
waitfor(image_handles.fig);

function voxel_list_callback(hObject, eventdata)
% Callback for 'voxel no.' get the selected voxel and update the plot according to the voxel selected
global details_para;
global image_handles;
global cur_vox;

% get the selected voxel
preprocess_sel_vox = get(hObject,'Value');

if(isempty(preprocess_sel_vox))
    preprocess_sel_vox = 1;
    set(image_handles.voxel_list,'value',1);
end

selected_vox = details_para.selected_image;
cur_vox = selected_vox(1,preprocess_sel_vox);
display_image()

function display_image()

global image_handles;
global cur_vox;
global details_para;
global data_value;

% img_kspace=data_value.FullData(:,:,cur_vox);
% d_kspace = max(max(max(abs(img_kspace))));
% img_image=ifft2c(img_kspace);
% d_image = max(max(abs(img_image)));
% img_kspace1=data_value.FullData1(:,:,cur_vox,:);
% img_kspace1=squeeze(permute((img_kspace1),[3 1 2 4]));
% d_kspace1 = max(max(max(abs(img_kspace))));
% for ll=1:size(img_kspace1,3)
%     img_image1(:,:,ll)=ifft2c(img_kspace1(:,:,ll));
% %     img_image2(:,:,ll)=(imresize(img_image1(:,:,ll),[8 8]));
% end
% d_image1 = max(max(max(abs(img_image1))));
% data_value.im_map_mod=(abs(img_image1)/d_image1);
% for i=1:size(data_value.im_map_mod,3)
%      temp_img_mod=(data_value.im_map_mod(:,:,i));
%      temp_img_mod=imresize(temp_img_mod,([900 600]));
%      data_value.im_map_final_mod(:,:,i)=imrotate(temp_img_mod,270*1);
% end
% data_value.final_image=sum(data_value.im_map_final_mod.^1,3);
imshow(data_value.final_image{cur_vox},'parent',image_handles.axis);

function close_Callback(hObject, eventdata)
%Callback for 'close' button
global image_handles;
global met_tick;
if (met_tick==1 && ~isempty(image_handles.met_fig))
    close(image_handles.fig)
    close(image_handles.met_fig);
else
     close(image_handles.fig);
end

function file_main_Callback(hObject, eventdata)
%Callback for 'first-order' edit box 
global pha_para;
global cur_vox;
global ind_allvox;
global image_handles;
global details_para;

function sense_Callback(hObject, eventdata)
%Callback for 'close' button
global image_handles;
global data_value;
global details_para;

cnt=0;
try
while(1)
    x = inputdlg(['Which slice do you want to use for the maps?']);
    waitfor(x);
    if((str2double(x)<=0) || (str2double(x)>size(data_value.FullData,3)) || (isnan(str2double(x))))
        cnt=cnt+1;
        h1(cnt)=errordlg('Invalid input');
    else
        slice_no=str2double(x);
        for i=1:(cnt) 
            close(h1(i));
        end
        break;
    end
end
set(image_handles.fig,'Visible','off');
img_kspace=data_value.FullData(:,:,slice_no);
d_kspace = max(max(max(abs(img_kspace))));
img_image=ifft2c(img_kspace);
d_image = max(max(abs(img_image)));
img_kspace1=data_value.FullData1(:,:,slice_no,:);
img_kspace1=squeeze(permute((img_kspace1),[3 1 2 4]));
d_kspace1 = max(max(max(abs(img_kspace))));
for ll=1:size(img_kspace1,3)
    img_image1(:,:,ll)=ifft2c(img_kspace1(:,:,ll));
%     img_image2(:,:,ll)=(imresize(img_image1(:,:,ll),[8 8]));
end
d_image1 = max(max(max(abs(img_image1))));
im_map_mod=(abs(img_image1)/d_image1);
im_map_real=(real(img_image1)/d_image1);%/ max(max(max(real(img_image1)))));
im_map_imag=(imag(img_image1)/d_image1);%/ max(max(max(imag(img_image1)))));
im_map_phase=(angle(img_image1));
for i=1:size(im_map_mod,3)
     temp_img_mod=(im_map_mod(:,:,i));
     temp_img_real=(im_map_real(:,:,i));
     temp_img_imag=(im_map_imag(:,:,i));
     temp_img_phase=(im_map_phase(:,:,i));
     im_map_final_mod(:,:,i)=imrotate(temp_img_mod,270*1);
     im_map_final_real(:,:,i)=imrotate(temp_img_real,270*1);
     im_map_final_imag(:,:,i)=imrotate(temp_img_imag,270*1);
     im_map_final_phase(:,:,i)=imrotate(temp_img_phase,270*1);
end
for i=1:size(im_map_final_mod,3)
    h=figure();
    set(h,'Name','RAPID 1.0','NumberTitle', 'off');
    if i==1
        [img2,rect]=imcrop((im_map_final_mod(:,:,i)));
    else
       [img2]=imcrop((im_map_final_mod(:,:,i)),rect);
    end
    if isempty(rect)
        data_value.image_error=1;
        data_value.sense_flag=0;
        set(image_handles.fig,'Visible','on');
        return;
    end
    [img2_real]=imcrop((im_map_final_real(:,:,i)),rect);
    [img2_imag]=imcrop((im_map_final_imag(:,:,i)),rect);
    [img2_phase]=imcrop((im_map_final_phase(:,:,i)),rect);
    img_sense(:,:,i)=(imresize(img2,[16/1 16/1]));
    img_sense_real(:,:,i)=(imresize(img2_real,[16/1 16/1]));
    img_sense_imag(:,:,i)=(imresize(img2_imag,[16/1 16/1]));
    img_sense_phase(:,:,i)=(imresize(img2_phase,[16/1 16/1]));
    img_sense_8(:,:,i)=(imresize(img2,[8/1 8/1]));
    img_sense_real_8(:,:,i)=(imresize(img2_real,[8/1 8/1]));
    img_sense_imag_8(:,:,i)=(imresize(img2_imag,[8/1 8/1]));
    img_sense_phase_8(:,:,i)=(imresize(img2_phase,[8/1 8/1]));
    close(h);
end
img_image2_16=img_sense;
img_image2_real_16=img_sense_real;
img_image2_imag_16=img_sense_imag;
img_image2_phase_16=img_sense_phase;
img_image2_8=img_sense_8;
img_image2_real_8=img_sense_real_8;
img_image2_imag_8=img_sense_imag_8;
img_image2_phase_8=img_sense_phase_8;
save('sense_map8_invivo.mat','img_image2_8','img_image2_real_8','img_image2_imag_8','img_image2_phase_8');
save('sense_map16_invivo.mat','img_image2_16','img_image2_real_16','img_image2_imag_16','img_image2_phase_16');
data_value.cropped_image=sum(img_image2_16,3);
% clear  data_value.final_image;
data_value.final_image{slice_no}=data_value.cropped_image;
data_value.sense_flag=1;
display_image();
catch
    data_value.image_error=1;
    data_value.sense_flag=0;
    set(image_handles.fig,'Visible','on');
    return;
end
set(image_handles.fig,'Visible','on');
%-----------------------------------------------------------------------------------------------------------------------------------
