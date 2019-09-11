function rapid_main_gui()

% This function is used for loading and preparing RAPID data with
% preprocessing

% Author: Sourav Bhaduri
% University of Ghent
% email: Sourav.Bhaduri@UGent.be
% Feb 2018; Last revision: Feb 2018

% addpath (strcat(pwd,'/regtools/'));
global details_para;
global data_value;
sig=cell(0);
rev_mode=0;
data_value.file_flag=0;
data_value.file_name=cell(0);
data_value.data_error=0;
for i=1:2
    try
        twix_obj = mapVBVD();
        close(data_value.h);
    catch
        data_value.data_error=1;
        if ((data_value.file_flag==1))
            return;
        end
        close(data_value.h);
        return;
    end
    filename=twix_obj{2}.image.filename;
    file_exten_ind=strfind(filename,'.');
    file_exten=filename(file_exten_ind:end);
    if (~strcmp(file_exten,'.dat') && ~strcmp(file_exten,'.mat'))
        data_value.data_error=1;
        return;
    end
    data_value.file_name{i}=filename(length(pwd)+1:end);
    data_value.spect_name=data_value.file_name{i};
    data=twix_obj{2}.image;
    try
        pos=twix_obj{2}.hdr.Phoenix.sAdjData.sAdjVolume;
    catch
        data_value.data_error=1;
        return;
    end
    FOV_RL=pos.dReadoutFOV;
    FOV_PE=pos.dPhaseFOV;
    FOV_HF=pos.dThickness;
    if isfield(pos,'sPosition')
        if isfield(pos.sPosition,'dSag')
            RL_off=pos.sPosition.dSag;
        else
            RL_off=0;
        end
        if isfield(pos.sPosition,'dCor')
            PE_off=pos.sPosition.dCor;
        else
            PE_off=0;
        end
        if isfield(pos.sPosition,'dTra')
            HF_off=pos.sPosition.dTra;
        else
             HF_off=0;
        end
    else
        RL_off=0;
        PE_off=0;
         HF_off=0;
    end
        
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
    has_header = isfield(twix_obj{2},'hdr');
    if has_header
        dwelltime = twix_obj{2}.hdr.MeasYaps.sRXSPEC.alDwellTime{1,1} * 1e-9;
        spectralwidth=1/dwelltime;
        if i==1
        PE=twix_obj{2}.hdr.Dicom.lPhaseEncodingLines;
        RE=twix_obj{2}.hdr.Dicom.lBaseResolution;
        end
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
    struct.nrows = twix_obj{2}.image.NAcq;
    rc_xres = double(twix_obj{2}.image.NCol);
    rc_yres = double(twix_obj{2}.image.NAcq);
    nreceivers = double(twix_obj{2}.image.NCha);
    FullData=permute(reshape(double(twix_obj{2}.image()),[twix_obj{2}.image.NCol twix_obj{2}.image.NCha twix_obj{2}.image.NPar twix_obj{2}.image.NAve twix_obj{2}.image.NSeg*twix_obj{2}.image.NPhs]),[2 1 5 4 3]);
    FullData=permute(FullData,[1,2,3,5,4]);
    FullData=sum(FullData,5)./twix_obj{2}.image.NAve;
    firstpoint=mean(conj(FullData(:,1,:)),3);
    channels_scale=squeeze(sqrt(sum(firstpoint.*conj(firstpoint))));
    firstpoint=repmat(firstpoint, [1 twix_obj{2}.image.NCol twix_obj{2}.image.NLin*twix_obj{2}.image.NPhs])/channels_scale;
    FullData = FullData*global_rescale;
    for ll=1:nreceivers
        if twix_obj{2}.image.NPar>1
            FullData1(:,:,:,ll) = conj(squeeze(sum(FullData(ll,:,:,:),1)));
        else
            FullData1(:,:,:,ll) = conj(squeeze(sum(FullData(ll,:,:),1)));
        end
    end
%     sig{i}= FullData1;
    sig{i}= (sum(FullData1,4));
    clear FullData1;
end
fit_param1=[];
Z_sig_total=[];D_seg3=[]; D_seg3_interp=[];
signal_tempo=sig{1};
cnt=0;
while(1)
    x = inputdlg(['Total slices: ',num2str(size(signal_tempo,3)), '   Which slice do you want to analyse?']);
    waitfor(x);
    if((str2double(x)<=0) || (str2double(x)>size(signal_tempo,3)) || (isnan(str2double(x))))
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
data_value.slice_profile=['_slice_',num2str(slice_no)];
cnt=0;
while(1)
    x = inputdlg(['Enter zerofill value: ']);
    waitfor(x);
    if((str2double(x)<=0) || (str2double(x)>10) || (isnan(str2double(x))))
        cnt=cnt+1;
        h1(cnt)=errordlg('Invalid input');
    else
        zerofill_value=str2double(x);
        for i=1:(cnt) 
            close(h1(i));
        end
        break;
    end
end
% PE=1;
h = waitbar(0,'Please wait...');
for ll=1:nreceivers./nreceivers
    clear fit_param1;clear sub_data_coeff_re;clear sub_data_coeff_im;
    clear Z_sig_total; clear sub_data_coeff; clear phase_diff; clear ref_param;
    clear D_seg3; clear D_seg3_interp; clear phase_param;clear phase_ref_param;
    signal=sig{1};signal=signal(:,:,slice_no,ll);
    signal_temp=signal;
    if size(signal,1)>size(signal,2)
        signal=(signal)';
    end
    % signal=randn(32,256);
%     RE=1;
    Rx=1;Ry=1;rp=2;
    AS=Rx.*size(signal,1)./RE;
    % AS=1;
    % RE=1;
    PE=PE./AS;
    PE=PE*AS;
    sig_temp=[];
    sig_temp1=signal;
    for i=1:RE*AS/Rx
        if AS==1 && RE==1
             sig_temp=[sig_temp,(signal(i,:))];
        else
            sig_temp=[sig_temp,conj(signal(i,:))];
        end
    end
    signal=sig_temp;   
    reference=sig{2};
    reference=reference(:,slice_no,1,ll);
    reference=(conj(reference'));
    fac2=Ry*Rx*size(signal,2)./(RE*AS)./(PE./AS); 
    blip_pt=2;
    % rev_mode=0;
    % load('test3.mat');
    % signal=(sub_signal');reference=(scale_half);
    % blip_pt=0;
    % fac2=size(sig_temp1,2)./AS;
    for ii=0:(RE/Rx)-1
        for kk=0:AS-1
            for l=1:(PE/AS/Ry)
                 sub_data(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry)) =  signal(1,1 + 1*fac2*(l-1) + (kk*size(sig_temp1,2))+ (ii*size(sig_temp1,2)*AS):1*fac2*l+(kk*size(sig_temp1,2))+(ii*size(sig_temp1,2)*AS));
                 temp=sub_data(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry));
                 temp(1:blip_pt)=[];
                 sub_data_part(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry))=ssa(temp,1);
            end
        end
    end
    % if rev_mode==1
    %     for i=1:RE
    %         for ii=1:PE
    %             sub_data_part_rev(:,ii+(i-1)*PE)= sub_data_part(:,(i)*PE- (ii+(i-1)*PE)+1 + (i-1)*PE);
    %         end
    %     end
    %     sub_data_part=sub_data_part_rev;
    % end
    pointsBeforeEcho=1;
    scale_half=reference(1, pointsBeforeEcho:end);
%     load('ref_signal_RAPID.mat');
%     scale_half=water_signal_CSI';
%     scale_half(1:blip_pt)=0;
%     sub_data_part(:,PE/2 +1+1*(RE*PE*0.5./Rx./Ry))=max(scale_half)*sub_data_part(:,PE/2 +1+1*(RE*PE*0.5./Rx./Ry))./max(sub_data_part(:,PE/2 +1+1*(RE*PE*0.5./Rx./Ry)));
%     scale_half(blip_pt+1:length(sub_data_part(:,1))+blip_pt)=sub_data_part(:,PE/2 +1+1*(RE*PE*0.5));
    water_signal_RAPID=scale_half';
    flag=round(pointsBeforeEcho./size(sub_data_part,1));
    water_supress=2;
    if water_supress==1
        half_temp_real=real(fftshift(fft(scale_half)));half_temp_imag=imag(fftshift(fft(scale_half)));
        half_temp_real_abs=abs(half_temp_real);half_temp_imag_abs=abs(half_temp_imag);
        half_temp_real_diff=half_temp_real_abs-half_temp_real;half_temp_imag_diff=half_temp_imag_abs-half_temp_imag;
        [half_real_recon,half_real_recon_vec,psi] = diff_shrodinger(1,half_temp_real_abs,30);
        half_real_recon1=(4/psi)*half_real_recon_vec(:,1);
        [half_imag_recon,half_imag_recon_vec,psi] = diff_shrodinger(1,half_temp_imag_abs,20);
        half_imag_recon1=(4/psi)*half_imag_recon_vec(:,1);
        scale_half_temp=scale_half;
        scale_half=ifft(ifftshift((half_temp_real-half_real_recon1')+ 1i*(half_temp_imag-half_imag_recon1')));
        for i=1:PE*size(signal,1)
            half_temp_real=real(fftshift(fft(sub_data_part(:,i)')));half_temp_imag=imag(fftshift(fft(sub_data_part(:,i)')));
            half_temp_real_abs=abs(half_temp_real);half_temp_imag_abs=abs(half_temp_imag);
            half_temp_real_diff=half_temp_real_abs-half_temp_real;half_temp_imag_diff=half_temp_imag_abs-half_temp_imag;
            [half_real_recon,half_real_recon_vec,psi] = diff_shrodinger(1,half_temp_real_abs,2);
            half_real_recon1=(4/psi)*half_real_recon_vec(:,1);
            [half_imag_recon,half_imag_recon_vec,psi] = diff_shrodinger(1,half_temp_imag_abs,2);
            half_imag_recon1=(4/psi)*half_imag_recon_vec(:,1);
            sub_data_part_temp(:,i)=sub_data_part(:,i);
            sub_data_part(:,i)=conj(ifft(ifftshift((half_temp_real-half_real_recon1')+ 1i*(half_temp_imag-half_imag_recon1'))));
         end      
    elseif water_supress==2
         model_order=15;
         Fs= spectralwidth;
         fres=(-( spectralwidth)/2 + (( spectralwidth)/(length(scale_half)*1))*(0:((length(scale_half)*1-1))));
         [fit_sig,amp,Ind_sig,freq] = HSVDPK_new(scale_half',fres,Fs,model_order) ;
         ind=find((freq>=-80 & freq<=80));
         half_sig2=sum(Ind_sig(:,ind),2);
         ind2=find((freq>=-800 & freq<=-350));
         half_sig3=sum(Ind_sig(:,ind2),2);
         scale_half_temp=scale_half;
         scale_half=scale_half- 1*half_sig2' - 0.0*half_sig3';
%          scale_fft=fftshift(fft(scale_half'));
%          baseline=ssa(scale_fft,20);
%          scale_half=scale_half-ifft(ifftshift(baseline'));
%          scale_half=ssa(scale_half',7);
%          scale_half=scale_half';
         for i=1:PE*RE./Ry./Rx
             model_order=10;
             fres=(-( spectralwidth)/2 + (( spectralwidth)/(length(sub_data_part(:,i))*1))*(0:((length(sub_data_part(:,i))*1-1))));
            [fit_sig,amp,Ind_sig,freq] = HSVDPK_new(sub_data_part(:,i),fres,Fs,model_order) ;
            ind=find((freq>=-80 & freq<=80));
            half_sig2=sum(Ind_sig(:,ind),2);
            ind2=find((freq>=-800 & freq<=-350));
            half_sig3=sum(Ind_sig(:,ind2),2);
            sub_data_part_temp(:,i)=sub_data_part(:,i);
            sub_data_part(:,i)=sub_data_part(:,i)- 1*half_sig2 - 0.0*half_sig3;
%             sub_data_part(:,i)=ssa(sub_data_part(:,i),4);
%             a= (pi*5);
%             sub_data_part(:,i)=exp((-a)*((1:length(sub_data_part(:,i)))'./Fs)).*sub_data_part(:,i);
            if max(abs(fftshift(fft(half_sig2))))<max(abs(fftshift(fft(sub_data_part(:,i)))))/(PE*RE)
               [fit_sig,amp,Ind_sig,freq] = HSVDPK_new(sub_data_part(:,i),fres,Fs,2) ;
                ind=find((freq>=-80 & freq<=80));
                half_sig2=sum(Ind_sig(:,ind),2);
                ind2=find((freq>=-800 & freq<=-350));
                half_sig3=sum(Ind_sig(:,ind2),2);
                sub_data_part(:,i)=sub_data_part(:,i)- 1*half_sig2 - 0.0*half_sig3;
%                 sub_data_part(:,i)=ssa(sub_data_part(:,i),4);
%                 sub_data_part(:,i)=exp((-a)*((1:length(sub_data_part(:,i)))'./Fs)).*sub_data_part(:,i);
%                 sub_data_part(:,i)=0;%sub_data_part(:,i)- 1*half_sig2;
%                 flag_ind(i)=0;
            end
         end      
    end
    org_half=scale_half;
    if size(signal,1)>size(reference,1)
         fl=0;
     else
         fl=1;
    end
    for i=1:PE/Ry
        if i==1
            shift_fac(i)=0;
        else
            shift_fac(i)=0;
        end
    end
    % scale_half=scale_half.*exp(-0*((1:512)/2000));
    % scale_half=scale_half_temp;
    fres=(-(Fs)/2 + ((Fs)/(length(scale_half)*1))*(0:((length(scale_half)*1-1))));
    scale_half_temp=scale_half;
    details_para.fla=0;
    time_sample=(1:length(scale_half))./Fs;
    details_para.seg = [1 length(scale_half)];
    Tf = 127776603;
    ref = 4.7;
    Acq_time=time_sample(end);
    ppm = (fres)*(1E6/Tf);
    details_para.ppm =  ppm;
    details_para.ppm_referenced = details_para.ppm + ref ;
    a= (pi*0);
    w = exp((-a)*time_sample);
    scale_half=scale_half.*w;
    if AS<=4
    scale_half=2*(PE/8)*(scale_half)./AS;
    end
%     for lll=1:size(sub_data_part,1)
%         temp_img_sub=reshape(sub_data_part(lll,:),PE/Ry,RE/Rx);
%         sub_data_part(lll,:)=reshape(temp_img_sub,PE*RE,1);
%     end
    % [No_pks,pk_loc,FWHM,signal_power_each] = find_peaks_new(scale_half',fres);
    % FWHM_ref=FWHM;
    % [fit_sig,~,amp_est] = model_varpro(scale_half',fres(pk_loc), FWHM, ones(1,length(pk_loc)),time_sample);
    % Initx = [amp_est' (pi*FWHM/2)' fres(pk_loc)'];
    % lb= [(amp_est-(amp_est/4))' ((pi*FWHM/2)-50)' (fres(pk_loc)-20)'];
    % ub= [(amp_est+(amp_est/4))' ((pi*FWHM/2)+50)' (fres(pk_loc)+20)'];
    % options = optimset('MaxFunEvals', 5000, 'MaxIter', 5000, 'Algorithm','trust-region-reflective', 'Display', 'off'); % optimization parameters 
    % [fit_param1,~,residual,~,~,~,J] = lsqnonlin(@(x)model_estimate_equation2(x,time_sample,(scale_half')),Initx,lb,ub,options); % non-linear least squares algorithm
    % %          fit_param(:,3)=fres1(pk_loc)';
    % [sig_sum,C,Z_sig] = model_estimate(fit_param1,time_sample);
    % FD_sig_wbc=fftshift(fft(scale_half'));
    % baseline_para=25;
    % blc_FD(:,1) = FD_sig_wbc; % for the first time baseline corrected signal (blc)= original FD signal
    % bs_lin(:,1) = FD_sig_wbc - FD_sig_wbc;
    % [err_val(1),~] = calculate_error(FD_sig_wbc,sig_sum,bs_lin(:,1),details_para.seg);
    % Max_ite = 20;
    % wh = waitbar(0,'Please wait');
    % lop = 2; % iteration starts from 2 because first iteration is done above 
    % while(1) % iterate between baseline correction and fitting
    % cmp_fit_sig = fftshift(fft(sig_sum)); % get the fitted signal by summing the individual signal 
    % [blc_TD,blc_FD(:,lop),bs_lin(:,lop)] = baseline_correction(FD_sig_wbc,cmp_fit_sig,details_para.ppm_referenced,baseline_para); %perform baseline correction
    % [fit_param1,~,residual,~,~,~,J] = lsqnonlin(@(x)model_estimate_equation2(x,time_sample,blc_TD),Initx,lb,ub,options); % non-linear least squares algorithm
    % [sig_sum,C,Z_sig] = model_estimate(fit_param1,time_sample);
    % [err_val(lop),~] = calculate_error(FD_sig_wbc,fftshift(fft(sig_sum)),bs_lin(:,lop),details_para.seg); %calculate error 
    %     if((abs(err_val(lop-1)-err_val(lop))< 1e-8)||(lop == Max_ite))
    %         no_ite = lop;
    %         break
    %     end
    %     waitbar(lop/Max_ite);
    %     lop = lop+1;
    % end
    % close(wh);
    % fit_param2=fit_param1;
%     a= (pi*0);
%     scale_half=exp((-a)*time_sample).*scale_half;
%     scale_rat=abs(max(sub_data_part(:,PE/2 +1+1*(RE*PE*0.5))))./abs(max(scale_half));
%     scale_half=scale_rat.*scale_half;
    model_order=PE-PE+4;
    Fs= spectralwidth;
    fres=(-( spectralwidth)/2 + (( spectralwidth)/(length(scale_half)*1))*(0:((length(scale_half)*1-1))));
    [fit_sig,amp,Ind_sig,freq,Z_sig,FWHM] = HSVDPK_new(scale_half',fres,Fs,model_order) ;
    ind=find((freq>=-700 & freq<=100));
    model_order=length(ind);
    Ind_sig=Ind_sig(:,ind);Z_sig=Z_sig(:,ind);
    freq=freq(ind); FWHM=FWHM(ind); amp=amp(ind);
    if PE==4
        F_rat=8;
    elseif PE==8
        F_rat=4;
    else
        F_rat=2;
    end
    % FWHM(2:3)=6*FWHM(2:3)./F_rat;
    % FWHM(4)=8*FWHM(4)./F_rat;
    % FWHM(5:end)=4*FWHM(5:end)./F_rat;
    %  FWHM(5:7)=2*FWHM(5:7);
    %  FWHM(8)=5*FWHM(8);
    %  FWHM(9:end)=7*FWHM(9:end);
    % details_para.fla=0;
    % [No_pks,pk_loc,FWHM_ref,signal_power_each] = find_peaks_new(sum(Ind_sig(:,2:8),2),fres);
    % FWHM(2:3)=FWHM_ref(3:4);
    % FWHM(4)=FWHM_ref(2);
    % FWHM(8)=FWHM_ref(1);
    % FWHM(6)=FWHM_ref(1);
    % FWHM(3)=3*FWHM(3);
    % FWHM(4)=FWHM(4)*4;
    % FWHM(1:2)=FWHM_ref(1:2);
    % FWHM(2)=FWHM_ref(4);%FWHM(3)=FWHM_ref(3);
    % FWHM(2)=FWHM(2)*2;
    % FWHM(4)=FWHM(4)*6;
    % FWHM(6)=FWHM(6)*2;
    % FWHM(4)=FWHM_ref(2);
    % FWHM(6)=FWHM_ref(1);
    % FWHM(6)=3*FWHM(6);
%     model_select=1;
%     Initx = [amp (FWHM)' freq'];
% %     lb= [(amp-(amp/4)) ((FWHM)-50)' (freq-20)'];
% %     ub= [(amp+(amp/4)) ((FWHM)+50)' (freq+20)'];
%      options = optimset('MaxFunEvals', 5000, 'MaxIter', 5000, 'Algorithm','trust-region-reflective', 'Display', 'off'); % optimization parameters 
%     [fit_param11,~,residual,~,~,~,J] = lsqnonlin(@(x)model_estimate_equation2(x,time_sample,scale_half',model_select),Initx,[],[],options); % non-linear least squares algorithm
%     [fit_sig,Ind_sig,Z_sig] = model_estimate(fit_param11,time_sample,model_select);
    Z_sig_total=Z_sig;
%     freso=(-( spectralwidth)/2 + (( spectralwidth)/(length(sub_data_part(:,i))*1))*(0:((length(sub_data_part(:,i))*1-1))));
%     [fit_sig_tempo,amp_tempo,Ind_sig_tempo,freq_tempo] = HSVDPK_new(conj(sub_data_part(:,PE/2 +1+1*(RE*PE*0.5))),freso,Fs,model_order) ;
%     freq_tempo2=freq;
%     freq=freq_tempo;
    for i=1:model_order
%         Z_sig(:,i)=exp((1i*2*pi*(freq(i))).*time_sample);
        Z_sig(:,i)=exp((-FWHM(i)/(4*pi/pi) + 1i*2*pi*(freq(i))).*time_sample);
%           Z_sig_total(:,i)=exp((-FWHM(i)*0/(8*pi) + 1i*2*pi*(freq(i))).*time_sample);
%         Z_sig(:,i)=Z_sig_total(:,i);
%         Z_sig(:,i)=(1/2)*exp((-FWHM(i)*2.5 + 1i*2*pi*(freq(i))).*time_sample);
%         if i==2 || i==3
%             a=(pi*2);
%             Z_sig(:,i)=exp((-a)*time_sample').* Z_sig(:,i);
%         end
%         Z_sig(:,i)=(1./pi)*Z_sig(:,i);
    end
%     Z_sig_total=Z_sig;
    fit_sig=sum(Ind_sig(:,1:end),2);
%     mean_scale(:,ll)=sum(Ind_sig(:,1:8),2)./(PE/1);
    noise_sig=scale_half'-fit_sig;
    nonbaseline_sig=scale_half'-Ind_sig(:,1);
    [fit_sig2,amp2,Ind_sig2,freq2,Z_sig2,FWHM2] = HSVDPK_new(nonbaseline_sig,fres,Fs,model_order) ;
    basis_rest=scale_half'-fit_sig;
    basis_rest=basis_rest./max(basis_rest);
    amp_basis=amp./1;
    T2_basis=1./(1*FWHM);
    noise_level=std(noise_sig(:));
    noise_sig=noise_level*randn(size(scale_half'));
    freq_met=freq;
    [amp_basis,phas_ref]=phas_correct(amp_basis,model_order);
    ref_sig=scale_half';
    ref_sig_total(:,:,ll)=scale_half';
    met_ind=[1:model_order];
    ind=~ismember(met_ind,[3]);
    ind_met=[3];
    count_space=0;
%     Z_sig=Ind_sig;
    for l = 1:PE/AS
                sub_ref(:,l)=ref_sig((1+blip_pt + fac2*(l-1-(flag/1)):fac2*(l-(flag/1))));
%                 modu_ref(:,l)=((sub_ref(:,1))./(sub_ref(:,l)));
%                 modu_ref(:,l)=ssa(modu_ref(:,l),20);
    end
    if PE==16
        amp_rat=.5*2;
    else
        amp_rat=.5*2;
%         amp_basis(ind_met-1)=.5*amp_basis(ind_met-1);
    end
%     amp_basis(1:16)=amp_rat*amp_basis(1:16);
%     lac_ind=find(max(abs(amp_basis(2:3))));
%     lac_ind=lac_ind+1;
%     Lac_ind_min=[2:3];
%     lac_ind_min=~ismember(Lac_ind_min,[lac_ind]);
%     lac_ind_min=Lac_ind_min(lac_ind_min==1);
%     amp_basis(ind_met)=0;
    % met_ind=met_ind(ind==1);
    % amp_basis=amp_basis./2;
    % for i=1:PE*RE
    %          model_order=20;
    %          fres=(-( spectralwidth)/2 + (( spectralwidth)/(length(sub_data_part(:,i))*1))*(0:((length(sub_data_part(:,i))*1-1))));
    %         [fit_sig,amp,Ind_sig,freq] = HSVDPK_new(sub_data_part(:,i),fres,Fs,model_order) ;
    %         ind=find((freq>=freq_met(2)-30 & freq<=freq_met(2)+30));
    %         half_sig_met_sig=(sum(Ind_sig(:,ind),2));
    %         half_sig_met(i)=max(abs(fftshift(fft(half_sig_met_sig))));
    % end      
    for i=1:size(T2_basis,2)
        for l=1:PE/Ry
            if (T2_basis(i)>(1*Acq_time/PE/Ry))
                    T2_flag(i,l)=0;
            else
                T2_flag(i,l)=1;
            end   
        end
    end
%     for i=1:size(sub_data_part,1)
%         sub_temp=conj(sub_data_part(i,:))';
%         sub_temp=fftshift(fft(sub_temp));
%         sub_temp2=ssa(sub_temp,20);
%         sub_data_part(i,:)=ifft(ifftshift(sub_temp2));
%     end
    for ii=0:(RE/Rx)-1
        for kk=0: AS-1  
            count_space=0;
            for l = 1:((PE/Ry)/AS)
                time_sample_seg=time_sample(1+blip_pt + fac2*(l-1-(flag/1)):fac2*(l-(flag/1)));
%                 n_var=estimatenoise((sub_data_part(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry))));
                if (mod(count_space,2)==1)
%                     sub_data_part(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry))=exp(1i*-50*pi/180).*(sub_data_part(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry)));
                end
                D_seg=Z_sig(((1+blip_pt + fac2*(l-1-(flag/1)):fac2*(l-(flag/1)))),[met_ind]);
                D_seg2=D_seg;
                count_space=count_space+1;
                if l==1
                    D_seg_temp=D_seg;
                end
               sub_sig_interp_store(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry))=(sub_data_part(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry)));
               for iiii=1:length(met_ind)
                    D_seg3(:,iiii)=D_seg(:,iiii);
               end
%                D_seg3=D_seg3.*(repmat(conj(amp_basis(met_ind))',size(D_seg3,1),1));
               for iiii=1:length(met_ind)
                     D_seg3_interp(:,iiii)=interp(((D_seg3(:,iiii))),1);
               end
               sub_sig_interp_store_interp(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry))=interp((sub_data_part(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry))),1);
               sub_ref(:,l)=ref_sig((1+blip_pt + fac2*(l-1-(flag/1)):fac2*(l-(flag/1))));
%                 rat_ref(:,l)=(angle(sub_ref(:,l)))- angle(abs(sub_ref(:,1)));
               rat_ref(:,l)=(((sub_ref(:,1))))./((sub_ref(:,l)));
%                pha_ref(:,l)=(angle(sub_ref(1,1))) - (angle(sub_ref(1,l)));
%                pha_ref(:,l)=exp(-pha_ref(:,l));
%                 rat_ref(:,l)=ssa(rat_ref(:,l),10);
                sub_ref_interp(:,l)=interp(sub_ref(:,l),1);
                basis_rest_seg=basis_rest(((1+blip_pt + fac2*(l-1-(flag/1)):fac2*(l-(flag/1)))));
                P = pseudoinverse(D_seg3_interp,[],'lsqr','tikhonov', {@(x,r) r*normest(D_seg3_interp)*x, 1e-18});
%                 sub_data_coeff(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry))=P*conj(sub_sig_interp_store_interp(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry)));
                sub_data_coeff_re(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry))=P*real(conj(sub_sig_interp_store_interp(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry))));
                sub_data_coeff_im(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry))=P*imag(conj(sub_sig_interp_store_interp(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry))));
%                 sub_data_coeff(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry))=D_seg3_interp\conj(sub_sig_interp_store_interp(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry)));
                sub_data_coeff(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry))=sub_data_coeff_re(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry))+(1i*sub_data_coeff_im(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry)));
%                 sub_data_coeff2(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry))=D_seg_temp\conj(sub_sig_interp_store_interp(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry)));
    %             con_mat(l+(kk*(PE/AS))+(ii*PE))=cond(D_seg,2);
    %             sig_reg(:,1)=lsqr(D_seg,conj(sub_data_part(:,l+(kk*(PE/AS))+(ii*PE))));
                lambda = .1;
                hh=0;
                sig_reg=[];error_regu=[];
    %             sub_data_coeff(:,l+(kk*(PE/AS))+(ii*PE))=lsqr(D_seg,conj(sub_data_part(:,l+(kk*(PE/AS))+(ii*PE))));
                w=[];
                for i=1:size(D_seg,2)
                  w=[w norm(D_seg(:,i))];
                end
                w=diag(w);
                w=inv(w);
                D_seg1=D_seg*w;
%                 ref_coeff=D_seg3_interp\sub_ref_interp(:,l);
                ref_coeff_re=P*real((sub_ref_interp(:,l)));
                ref_coeff_im=P*imag((sub_ref_interp(:,l)));
                ref_coeff=ref_coeff_re + (1i*ref_coeff_im);
%                 ref_coeff2=D_seg_temp\sub_ref_interp(:,l);
                fit_param1(:,1)=(sub_data_coeff(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry)));
                [~,phas_ref]=phas_correct(ref_coeff,size(fit_param1,1));
                [ref_param(:,1),phas_ref]=phas_correct(ref_coeff,size(fit_param1,1));
%                 [ref_param2(:,l),phas_ref2]=phas_correct(ref_coeff2,size(fit_param1,1));
%                 fit_param2(:,l)=(sub_data_coeff2(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry)));
%                 [fit_param2(:,l),phas2]=phas_correct(fit_param2(:,l),size(fit_param2,1));
                if l==1
                    phas_ref0=phas_ref;
                end
                [fit_param1(:,1),phas]=phas_correct(fit_param1(:,1),size(fit_param1,1));
                phase_diff(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry))=(phas-phas_ref)';
                phase_param(l+(kk*(PE/AS/Ry))+(ii*PE/Ry),:)=phas;
                phase_ref_param(l+(kk*(PE/AS/Ry))+(ii*PE/Ry),:)=phas_ref;
%                 ref_param_rat2(:,l)=(ref_param2(:,l))./(ref_param2(:,1));
                fit_param1(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry))=fit_param1(:,1);
%                 ssignal(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry))=D_seg3_interp(:,[met_ind])*fit_param1(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry));
%                 model_select=1;   
%                 Initx = [(fit_param1(:,1)) FWHM' freq'];
            %     lb= [(amp-(amp/4)) ((FWHM)-50)' (freq-20)'];
            %     ub= [(amp+(amp/4)) ((FWHM)+50)' (freq+20)'];
%                  options = optimset('MaxFunEvals', 50, 'MaxIter', 50, 'Algorithm','trust-region-reflective', 'Display', 'off'); % optimization parameters 
%                 [fit_param_temp,~,residual,~,~,~,J] = lsqnonlin(@(x)sub_estimate_equation(x,D_seg3_interp(:,[met_ind]),conj(sub_sig_interp_store_interp(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry)))),Initx,[],[],options); % non-linear least squares algorithm
%                 [fit_param_temp,~,residual,~,~,~,J] = lsqnonlin(@(x)sub_estimate_equation(x,time_sample,met_ind,blip_pt,fac2,flag,model_order,l,conj(sub_sig_interp_store_interp(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry)))),Initx,[],[],options); % non-linear least squares algorithm
%                 ssignal2(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry))=D_seg3_interp(:,[met_ind])*fit_param_temp(:,1);
%                 fit_param1(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry))=fit_param_temp(:,1);
            end
        end
    end
    for i=1:model_order
        phas_mean=mean(phase_param(:,i));
%         phase_param(:,i)=ssa(phase_param(:,i),2);%(phase_param(:,i));
%         p = polyfit([1:1:size(phase_param,1)]',phase_param(:,i),6);
%         y1 = polyval(p,[1:1:size(phase_param,1)]');
%         phase_param(:,i)=y1;
%         phase_param(:,i)= ssa(phase_param(:,i),4);%-mean(phase_param(:,i));
%         phase_param(1:3,i)=0;
%         phase_param(end:end-3,i)=0;
          phase_param(:,i)= (unwrap((phase_param(:,i))));%
          phase_ref_param(:,i)= (unwrap((phase_ref_param(:,i))));%
          val_par=conj(fit_param1(i,:)');
%           fit_param1_real(i,:)=ssa(real(val_par),1);
%           fit_param1_imag(i,:)=ssa(imag(val_par),1);
%           fit_param1(i,:)=fit_param1_real(i,:)+(1i*fit_param1_imag(i,:));
%           phase_param(:,i)=unwrap(angle(fit_param1(i,:)));
%           phase_param(:,i)=ssa(phase_param(:,i),5);
%           vector = phase_param(:,i);
%           percntiles = prctile(vector,[5 95]); %5th and 95th percentile
%           outlierIndex = vector < percntiles(1) | vector > percntiles(2);
%           %remove outlier values
%           vector(outlierIndex) = 0;
%           phase_param(:,i)=vector;
    end
    for ii=0:(RE/Rx)-1
        for kk=0: AS-1  
            count_space=0;
            for l = 1:((PE/Ry)/AS)
                 fit_param1(:,1)=abs(fit_param1(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry))).*exp(1i*(phase_param(l+(kk*(PE/AS/Ry))+(ii*PE/Ry),:) - 0*phase_ref_param(l+(kk*(PE/AS/Ry))+(ii*PE/Ry),:))');%.*abs(amp_basis(met_ind))./abs(fit_param1(:,1));
%                 fit_param1(:,1)=abs(amp_basis).*exp(1i*(phase_param(l+(kk*(PE/AS/Ry))+(ii*PE/Ry),:)- 1*phas_ref)');
                 for iii=1:model_order
                     if(abs(fit_param1(iii,1))>=1.5*abs(amp_basis(iii)))% && abs(fit_param1(iii,1))<=noise_level)% || T2_flag(iii,l)==1)
%                          fit_param1(iii,1)=0;
                     end
                 end
%                  for iii=1:length(fit_param2(:,1))
%                      if(abs(fit_param2(iii,l))>=1.5*abs(amp_basis(iii)))% || T2_flag(iii,l)==1)% && abs(fit_param1(iii,1))<=noise_level)% || T2_flag(iii,l)==1)
%                          fit_param2(iii,1)=0;
%                      end
%                  end
                 sig_sum_temp(:,1)=1*Z_sig_total(:,[met_ind])*fit_param1(:,1);
                 sig_sum=sum(sig_sum_temp,2);
                 sig_sum3=((conj(sub_data_part(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry)))));%.*modu_ref(:,l);
                 sig_sum2=zeros(size(sig_sum));%sig_sum;
%                  sig_sum2=scale_half';
                 sig_sum4=zeros(size(sig_sum));%sig_sum;
%                  sig_sum2=(scale_half)'.*exp(1i*(max(angle(conj(sub_data_part(:,l+(kk*(PE/AS))+(ii*PE)))))-max(angle((scale_half)'))));
%                  sig_sum(((1+blip_pt + fac2*(l-1-(flag/1)):fac2*(l-(flag/1)))))=(conj(sub_data_part(:,l+(kk*(PE/AS))+(ii*PE))));
                 sig_sum2(((1+blip_pt + fac2*(l-1-(flag/1)):fac2*(l-(flag/1)))))=conj(sub_data_part(:,l+(kk*(PE/AS))+(ii*PE)));
                 sig_sum4(1+blip_pt :fac2)=sig_sum3;%.*rat_ref(:,l);%.*rat_ref(:,l);%./ref_param_rat2(1,l);
                 extrapolated_signal(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry))=sig_sum;
                 extrapolated_signal2(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry))=sig_sum2;
                 if l+(kk*(PE/AS))+(ii*PE)==PE/2 +1+1*(RE*PE*0.5) 
                         extrapolated_signal(:,l+(kk*(PE/AS))+(ii*PE))=  (scale_half)';
                         extrapolated_signal2(:,l+(kk*(PE/AS))+(ii*PE))=  (scale_half)';
                 end
%                  [ fit_sig,amp,Ind_sig,freq] = HSVDPK_new(extrapolated_signal(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry)),fres,Fs,2) ;
%                  ind=find((freq>=-80 & freq<=80));
%                  half_sig2=sum(Ind_sig(:,ind),2);
%                  extrapolated_signal(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry))=extrapolated_signal(:,l+(kk*(PE/AS/Ry))+(ii*PE/Ry))-half_sig2;
            end
        end
    end
%     extrapolated_signal(:,1:2:end)=0;
    net_blipped_sig_org1=zeros(PE*RE./Ry./Rx,length(scale_half_temp));
    net_blipped_sig_org1_2=net_blipped_sig_org1;
    extrapolated_signal3=extrapolated_signal(shift_fac(1)+1:end,:);
    net_blipped_sig_org1(:,1:length(scale_half_temp))=conj(extrapolated_signal)';
    net_blipped_sig_org1_2(:,1:length(scale_half_temp))=conj(extrapolated_signal2)';
    net_blipped_sig_org2=zeros(PE*RE*zerofill_value.^2,length(scale_half_temp));
    net_blipped_sig_org2_2=zeros(PE*RE,length(scale_half_temp));
    half_fourier=0;
%----------------------------HALF FOURIER-----------------------------------------------------------------% 
if half_fourier==1
    net_blipped_sig_org1(1:PE,:)=0;
%     for kkk=2:(RE/2+1)
%         for l = 1:PE/2 -1
%           part_re=(real((((net_blipped_sig_org1(((RE-kkk+2)*PE)-l+1,:))))));
%           part_im=(imag((((net_blipped_sig_org1(((RE-kkk+2)*PE)-l+1,:))))));
%           part_im=-(part_im);
%           net_blipped_sig_org1(l+((kkk-1)*PE)+1,:)= part_re + (1i*part_im);
%         end
%     end
    for kkk=2:(RE/2)
        for l = 1:PE/2 -1
          part_re=(real((((net_blipped_sig_org1(((RE-kkk+2)*PE)-l+1,:))))));
          part_im=(imag((((net_blipped_sig_org1(((RE-kkk+2)*PE)-l+1,:))))));
          part_im=-(part_im);
          net_blipped_sig_org1(l+((kkk-1)*PE)+1,:)= part_re + (1i*part_im);
        end
    end
    for kkk=2:(RE/2)
        for l = 1:PE/2-1 
          part_re=(real((((net_blipped_sig_org1(((RE-kkk+2)*PE)-l+1-(PE/2),:))))));
          part_im=(imag((((net_blipped_sig_org1(((RE-kkk+2)*PE)-l+1-(PE/2),:))))));
          part_im=-(part_im);
          net_blipped_sig_org1(l+((kkk-1)*PE)+1+(PE/2),:)= part_re + (1i*part_im);
        end
    end
    for kkk=(RE/2+2):RE
        for l =  PE/2+1  
          part_re=(real((((net_blipped_sig_org1(((RE-kkk+2)*PE)-l+2,:))))));
          part_im=(imag((((net_blipped_sig_org1(((RE-kkk+2)*PE)-l+2,:))))));
          part_im=-(part_im);
          net_blipped_sig_org1(((RE-(RE-kkk+2)+2)*PE)-l+2,:)= part_re + (1i*part_im);
        end
    end
    for kkk=(RE/2+2):RE+1
        for l = PE/2+1 
          net_blipped_sig_org1(((RE-kkk+2)*PE)-l+1-(PE/2)+1,:)= zeros(1,1*length(scale_half));
        end
    end
end
%------------------------------------------------------------------------------------------------------------%
    for i=1:((RE*PE)./Rx./Ry)
        net_blipped_sig_org1(i,:)=(fftshift(fft(net_blipped_sig_org1(i,:))))./((RE*PE)/(Rx*Ry));
        net_blipped_sig_org1(i,:)=conj(net_blipped_sig_org1(i,:));
        net_blipped_sig_org1_2(i,:)=(fftshift(fft(net_blipped_sig_org1_2(i,:))))./((RE*PE)/(Rx*Ry));
        net_blipped_sig_org1_2(i,:)=conj(net_blipped_sig_org1_2(i,:));
        SNR_each(i)=max(abs(fftshift(fft(conj(sub_data_part(:,i))))));
    end
    % SNR_each(9)=[];
    SNR=mean(SNR_each);
    vox_RL=round(RE*RL_off./FOV_RL);
    vox_PE=round(PE*PE_off./FOV_PE);
    for kk = 1:length(scale_half)*1    %fac
                temp1=conj(net_blipped_sig_org1(:,kk))';%half_signal(kk,:);%net_blipped_sig_org1(:,kk);%half_signal(kk,:);%net_blipped_sig_org1(:,kk);%half_signal(kk,:);%net_blipped_sig_org1(:,kk);%half_signal(kk,:);%net_blipped_sig_org1(:,kk);%half_signal(kk,:);%net_blipped_sig_org1(:,kk);%half_signal(kk,:);%net_blipped_sig_org1(:,kk);%half_signal(kk,:);%net_blipped_sig_org1(:,kk);%half_signal(kk,:);%net_blipped_sig_org1(:,kk);% half_signal(kk,:);%net_blipped_sig_org1(:,kk);
                 temp1_2=conj(net_blipped_sig_org1_2(:,kk))';
                 if size(signal_temp,1)>size(signal_temp,2)
    %                 temp1=conj(temp1)';
                 else
    %                 temp1=conj(temp1)';
                 end
                temp1_for=temp1((PE/2)+1:end);temp1_for=temp1_for';
                temp1_back=temp1(1:(PE/2));temp1_back=temp1_back';
    %             [temp1_for]=auto_regression_arma(temp1_for,zeros(length(zeros(1,((PE*size(signal,1)*(zerofill_value))-length(temp1))./2))+length(temp1_for),1),(PE/2)-1,temp1_for);
    %             [temp1_back]=auto_regression_arma(flipud(temp1_back),zeros(length(zeros(1,((PE*size(signal,1)*(zerofill_value))-length(temp1))./2))+length(temp1_back),1),(PE/2)-1,flipud(temp1_back));
    %             temp1_back=fliplr(temp1_back);
    %             temp1=conj([temp1_back,temp1_for]);
                temp1=[zeros(1,(((PE/Ry)*(RE/Rx)*(1.^2))-length(temp1))./2),temp1,zeros(1,(((PE/Ry)*(RE/Rx)*(1.^2))-length(temp1))./2)];
                temp1_2=[zeros(1,(((PE/Ry)*(RE/Rx)*(1.^2))-length(temp1))./2),temp1_2,zeros(1,(((PE/Ry)*(RE/Rx)*(1.^2))-length(temp1))./2)];
                temp=reshape(temp1,(PE/Ry)*1,(RE/Rx)*1);
                tmp=temp;
%                 temp_2=reshape(temp1_2,(PE/Ry)*zerofill_value,(RE/Rx)*zerofill_value);
%                 tmp_2=temp_2;
                for ii=1:RE/Rx
%                     tmp(:,ii)=(tmp(:,ii)).*hamming(length(tmp(:,ii)));
%                     tmp(:,ii)=ssa(tmp(:,ii),2.*PE./8);
%                     tmp_real(:,ii)=ssa(real(tmp(:,ii)),1);%2.*PE./8);
%                     tmp_imag(:,ii)=ssa(imag(tmp(:,ii)),1);%2.*PE./8);
%                     tmp_real(:,ii)=image_smooth((real(tmp(:,ii))),0.1);%2.*PE./8);
%                     tmp_imag(:,ii)=image_smooth((imag(tmp(:,ii))),0.1);%2.*PE./8);
%                     tmp(:,ii)=tmp_real(:,ii)+(1i*tmp_imag(:,ii));
                      tmp(:,ii)=image_smooth(((tmp(:,ii))),0.6);%2.*PE./8);
%                     tmp_2(:,ii)=(tmp_2(:,ii)).*hamming(length(tmp(:,ii)));
%                     tmp_2(:,ii)=ssa(tmp_2(:,ii),2.*PE./8);
%                     tmp_zero(:,ii)=[zeros(1,(((PE/Ry)*(zerofill_value))-length(tmp(:,ii)))./2)';tmp(:,ii);zeros(1,(((PE/Ry)*(zerofill_value))-length(tmp(:,ii)))./2)'];
                end
                for ii=1:PE/Ry
                    tmp(ii,:)=(tmp(ii,:)).*hamming(length(tmp(ii,:)))';
%                     tmp(ii,:)=ssa(tmp(ii,:)',2.*PE./8);
                    tmp_real(ii,:)=ssa(real(tmp(ii,:)'),1);%(2./Rx).*PE./8);
                    tmp_imag(ii,:)=ssa(imag(tmp(ii,:)'),1);%(2./Rx).*PE./8);
                    tmp(ii,:)= tmp_real(ii,:) + (1i* tmp_imag(ii,:));
                    tmp(ii,:)=conj(tmp(ii,:));
%                     tmp_2(ii,:)=(tmp_2(ii,:)).*hamming(length(tmp(ii,:)))';
%                     tmp_2(ii,:)=ssa(tmp_2(ii,:)',2.*PE./8);
%                     tmp_2(ii,:)=conj(tmp_2(ii,:));
                end
                tmp_zero=zeros(length(tmp(:,ii))*zerofill_value,length(tmp(:,ii))*zerofill_value);
                tmp_zero((length(zeros(1,(((PE/Ry)*(zerofill_value))-length(tmp(:,ii)))./2)')+1):(length(zeros(1,(((PE/Ry)*(zerofill_value))-length(tmp(:,ii)))./2)')+(PE/Ry)),(length(zeros(1,(((RE/Rx)*(zerofill_value))-length(tmp(:,ii)))./2)')+1):(length(zeros(1,(((RE/Rx)*(zerofill_value))-length(tmp(:,ii)))./2)')+(RE/Rx)))=tmp;
%                 tmp=image_smooth(tmp,.45);
    %             for ii=1:RE      
    %                         tmp_y(:,ii)=temp1(:,(ii-1)*PE+1 : (ii*PE));
    %                         tmp_y(:,ii)=tmp_y(:,ii).*hamming(length(tmp_y(:,ii)));
    %                         tmp_y(:,ii) =(fftshift(fft2(((tmp_y(:,ii))))));
    %                         if PE==4
    %                             tmp_y(:,ii) =circshift(tmp_y(:,ii) ,[0 0]);
    %                         elseif PE==8
    %                             tmp_y(:,ii) =circshift(tmp_y(:,ii) ,[-1 0]);
    %                         elseif PE==16
    %                             tmp_y(:,ii) =circshift(tmp_y(:,ii) ,[-2 0]);
    %                         end
    %             end 
    %             for ii=1:PE 
    %                 tmp_x(:,ii)=(((((tmp_y(ii,:))))));
    %                 tmp_x(:,ii) =(fftshift(fft2(((tmp_x(:,ii))))));
    %             end
    %             tmp=tmp_x;
    %             for ii=1:RE       
    %                     tmp_y(:,ii)=tmp(ii,:);
    % %                     tmp_y(:,ii)=tmp_y(:,ii).*hamming(length(tmp(ii,:)));
    %                     tmp_y(:,ii) =((fft(((tmp_y(:,ii))))));
    %                     if PE==4
    %                         tmp_y(:,ii) =circshift(tmp_y(:,ii) ,[0 0]);
    %                     elseif PE==8
    %                         tmp_y(:,ii) =circshift(tmp_y(:,ii) ,[0 0]);
    %                     elseif PE==16
    %                         tmp_y(:,ii) =circshift(tmp_y(:,ii) ,[-2 0]);
    %                     end
    %             end
    %             for ii=1:PE       
    %                     tmp_x(:,ii)=tmp_y(ii,:);
    %                     tmp_x(:,ii) =(fftshift(fft(((tmp_x(:,ii))))));
    %             end
    %             tmp=tmp_x;
    %             tmp=tmp.*hamming(length(tmp));
    %            tmp(1:RE,1:2)=0;
               tmp_full=zeros(PE*zerofill_value,RE*zerofill_value);
               tmp_full(1:Ry:PE*zerofill_value,1:Rx:RE*zerofill_value)=tmp_zero;
%                tmp_full_2=zeros(PE,RE);
%                tmp_full_2(1:Ry:PE,1:Rx:RE)=tmp_2;
%                tmp_full(1:1:PE,2:2:RE)=0;
               tmp_zero=tmp_full;
               tmp_zero1 = tmp_zero;%image_smooth((tmp_zero),0.5);
               tmp_zero =(fftshift(fft2(((tmp_zero1)))));
%                tmp_zero = image_smooth((tmp_zero),0.06);
%                tmp_2=tmp_full_2;
%                tmp_2 =(fftshift(fft2(((tmp_2)))));
               tmp_interp=tmp_zero;
%                tmp_interp_2=tmp_2;
               tmp_zero=tmp_zero(1:PE*zerofill_value,1:RE*zerofill_value);
%                tmp_2=tmp_2(1:zerofill_value:PE*zerofill_value,1:RE);
    %             if size(signal_temp,1)>size(signal_temp,2)
    %                 tmp=circshift(tmp,[1 0]);
    %                 tmp_interp=circshift(tmp_interp,[0 0]);
    %             end
                if rev_mode==1
                    sh=1; si=1;
                else
                    sh=-2;si=-1;
                end
                tmp_zero=circshift(tmp_zero,[0 0]);
                tmp_zero=circshift(tmp_zero,[vox_PE*zerofill_value -vox_RL*zerofill_value]);
                tmp_interp=circshift(tmp_interp,[-1 -1]);
                for ii=1:RE*zerofill_value
                    if rev_mode==1
    %                 tmp(:,ii)=flipud(tmp(:,ii));
                    else
                        tmp_zero(:,ii)=flipud(tmp_zero(:,ii));
%                         tmp_2(:,ii)=flipud(tmp_2(:,ii));
                    end
                end
                for ii=1:RE
%                     tmp(:,ii)=fliplr(tmp(ii,:));
%                     tmp(ii,:)=fliplr(tmp(ii,:));
                end
                for ii=1:PE
    %                 tmp(ii,:)=(tmp(:,ii));
                end
                temp=(tmp_zero);
                temp_interp=tmp_interp;
%                 temp_2=(tmp_2);
%                 temp_interp_2=tmp_interp_2;
                net_blipped_sig_org2(:,kk) = reshape(temp,PE*RE*zerofill_value*zerofill_value,1); %circshift(tmp,-1);
%                 net_blipped_sig_org2_2(:,kk) = reshape(temp_2,PE*RE,1); %circshift(tmp,-1);
    %             net_blipped_sig_org_interp(:,kk) = reshape(temp_interp,PE*RE*zerofill_value,1);
                tmp = 0;
%                 tmp_2 = 0;
    end
    PE=PE*zerofill_value;RE=RE*zerofill_value;
    details_para.PE=PE;details_para.RE=RE;
    for i=1:PE*RE
        X_temp=((fftshift(fft(scale_half'))))./(PE*RE./Rx./Ry);%-abs((fftshift(fft(fit_sig))));
        noise_X=std(X_temp(end-50:end));
        Y_temp=abs(net_blipped_sig_org2(i,:));
        if max(abs(Y_temp(1:230)))<=1.5*noise_X
    %         net_blipped_sig_org1(i,:)=max(abs(((net_blipped_sig_org1(i,:))))).*randn(size(net_blipped_sig_org1(i,:)))./(PE*RE);
        end
%         recon_sig(:,i)=ifft(ifftshift((net_blipped_sig_org2(i,:))));
        recon_sig(:,i)=((net_blipped_sig_org2(i,:)'));
        recon_sig(:,i)=recon_sig(:,i);
        recon_sig_phase(:,i)=angle((net_blipped_sig_org2(i,:)));
%         recon_sig_2(:,i)=abs((net_blipped_sig_org2_2(i,:)));
%         recon_sig_2(:,i)=recon_sig_2(:,i);
%         recon_sig_phase_2(:,i)=angle((net_blipped_sig_org2_2(i,:)));
    end
    recon_sig_tempo=recon_sig;
    recon_sig_tempo_phase=recon_sig_phase;
%     recon_sig_tempo_2=recon_sig_2;
%     recon_sig_tempo_phase_2=recon_sig_phase_2;
    for ii=0:PE-1
            for l=1:(RE)
                 recon_sig(:,l+(ii*RE)) =  recon_sig_tempo(:,(l-1)*RE +(ii)+ 1);
                 recon_sig_phase(:,l+(ii*RE)) =  recon_sig_tempo_phase(:,(l-1)*RE +(ii)+ 1);
%                  recon_sig_2(:,l+(ii*RE)) =  recon_sig_tempo_2(:,(l-1)*RE +(ii)+ 1);
%                  recon_sig_phase_2(:,l+(ii*RE)) =  recon_sig_tempo_phase_2(:,(l-1)*RE +(ii)+ 1);
            end
    end
    % for i=1:PE*RE*zerofill_value
    %     recon_sig_interp(:,i)=ifft(ifftshift((net_blipped_sig_org_interp(i,:))));
    % end
    fres=(-( spectralwidth)/2 + (( spectralwidth)/(length(scale_half)*1))*(0:((length(scale_half)*1-1))));
    details_para.seg = [1 length(scale_half)];
    Tf = 127776603;
    ref = 4.7;
    ppm = (fres)*(1E6/Tf);
    details_para.ppm =  ppm;
    details_para.ppm_referenced = details_para.ppm + ref ;
    seg=round(length(scale_half)/4)-80:length(scale_half)/2 + 100;
    recon_sig_total(:,:,ll)=recon_sig;
    recon_sig_total_phase(:,:,ll)=recon_sig_phase;
    recon_sig_map(:,:,ll)=squeeze(reshape(recon_sig_total(1,:,ll),RE,PE));
%     recon_sig_total_2(:,:,ll)=recon_sig_2;
%     recon_sig_total_phase_2(:,:,ll)=recon_sig_phase_2;
%     recon_sig_map_2(:,:,ll)=squeeze(reshape(recon_sig_total_2(1,:,ll),RE,PE));
    for lll=1:length(scale_half)
        recon_sig_ima(:,:,ll,lll)=squeeze(reshape(recon_sig_total(lll,:,ll),RE,PE));
        recon_sig_ima_phase(:,:,ll,lll)=squeeze(reshape(recon_sig_total_phase(lll,:,ll),RE,PE));
%         recon_sig_ima_2(:,:,ll,lll)=squeeze(reshape(recon_sig_total_2(lll,:,ll),RE,PE));
%         recon_sig_ima_phase_2(:,:,ll,lll)=squeeze(reshape(recon_sig_total_phase_2(lll,:,ll),RE,PE));
%         recon_sig_map(:,:,ll,lll)=squeeze(reshape(recon_sig_total(lll,:,ll),RE,PE));
    end
    waitbar(1);
end
close(h);
arg=1;
% load('sense_map.mat','img_image2');
% recon_sig_map=img_image2;
% sosimg = sqrt(sum(abs(recon_sig_map.^2),3));
% recon_sig_map= (recon_sig_map./repmat(sosimg,[1 1 nreceivers]));
dir_flag=1;
% [recon_sig_map]=sensitivity_map([RE ,PE],nreceivers);
% load('PI_map.mat','recon_sig_map');
% load('sense_map.mat','recon_sig_map');
%-------------------------------------------------------------------------------------------
% load('sense_map8_invivo.mat','img_image2','img_image2_real','img_image2_imag','img_image2_phase');
% recon_sig_map=img_image2;
% recon_sig_map_phase=img_image2_phase;
% recon_sig_map_2=img_image2;
% recon_sig_map_phase_2=img_image2_phase;
% for lll=1:length(scale_half)
%     for ll=1:nreceivers
%         temp_img=(recon_sig_ima(:,:,ll,lll));
%         temp_img_phase=(recon_sig_ima_phase(:,:,ll,lll));
%         temp_img4=(recon_sig_map(:,:,ll,1));
%         temp_img_2=(recon_sig_ima_2(:,:,ll,lll));
%         temp_img_phase_2=(recon_sig_ima_phase_2(:,:,ll,lll));
%         temp_img4_2=(recon_sig_map_2(:,:,ll,1));
%         for ii=1:PE
%             temp_img1(ii,:)=temp_img(:,ii);
%             temp_img5(ii,:)=temp_img4(:,ii);
%             temp_img1_phase(ii,:)=temp_img_phase(:,ii);
%             temp_img1_2(ii,:)=temp_img_2(:,ii);
%             temp_img5_2(ii,:)=temp_img4_2(:,ii);
%             temp_img1_phase_2(ii,:)=temp_img_phase_2(:,ii);
%         end
%         for ii=1:RE
%             temp_img1(:,ii)=temp_img(ii,:);
%             temp_img5(:,ii)=temp_img4(ii,:);
%             temp_img1_phase(:,ii)=temp_img_phase(ii,:);
%             temp_img1_2(:,ii)=temp_img_2(ii,:);
%             temp_img5_2(:,ii)=temp_img4_2(ii,:);
%             temp_img1_phase_2(:,ii)=temp_img_phase_2(ii,:);
%         end
%     recon_sig_ima(:,:,ll,lll)=temp_img1;
%     recon_sig_ima_phase(:,:,ll,lll)=temp_img1_phase;
%     recon_sig_map2(:,:,ll,lll)=temp_img5;
%     recon_sig_ima_2(:,:,ll,lll)=temp_img1_2;
%     recon_sig_ima_phase_2(:,:,ll,lll)=temp_img1_phase_2;
%     recon_sig_map2_2(:,:,ll,lll)=temp_img5_2;
%     end
% end
% %---------------------------------------------------------------------------------------------------
% for lll=1:length(scale_half)
%     im_full(:,:,lll) = sense((recon_sig_ima(:,:,:,lll)),(recon_sig_map(:,:,:,1)), 1, Rx*1);
%     im_full_phase(:,:,lll) = sense((recon_sig_ima_phase(:,:,:,lll)),(recon_sig_map_phase(:,:,:,1)), 1, Rx*1 );
%     im_full(:,:,lll)=im_full(:,:,lll).*exp(1i*im_full_phase(:,:,lll));
%     im_full_2(:,:,lll) = sense((recon_sig_ima_2(:,:,:,lll)),(recon_sig_map_2(:,:,:,1)), 1, Rx*1);
%     im_full_phase_2(:,:,lll) = sense((recon_sig_ima_phase_2(:,:,:,lll)),(recon_sig_map_phase_2(:,:,:,1)), 1, Rx*1 );
%     im_full_2(:,:,lll)=im_full_2(:,:,lll).*exp(1i*im_full_phase_2(:,:,lll));
% %     im_full(:,:,lll) = SENSE_recon((recon_sig_map),(recon_sig_ima(:,1:rp/1:RE,:,lll)),rp/1,arg,dir_flag);
% end
% for lll=1:length(scale_half)
%         temp_img=(im_full(:,:,lll));
%         temp_img_2=(im_full_2(:,:,lll));
%         for ii=1:PE
%             temp_img1(ii,:)=temp_img(:,ii);
%             temp_img1_2(ii,:)=temp_img_2(:,ii);
%         end
%         for ii=1:RE
%             temp_img1(:,ii)=temp_img(ii,:);
%             temp_img1_2(:,ii)=temp_img_2(ii,:);
%         end
%        im_full(:,:,lll)=temp_img1;
%        im_full_2(:,:,lll)=temp_img1_2;
% end
% recon_sig = reshape(im_full,PE*RE,length(scale_half));
% recon_sig=conj(recon_sig)';%./(Rx*Ry);
% recon_sig_2 = reshape(im_full_2,PE*RE,length(scale_half));
% recon_sig_2=conj(recon_sig_2)';%./(Rx*Ry);
if AS<=4
recon_sig=recon_sig./(2*(PE/8)./AS);
% recon_sig_2=recon_sig_2./(2*(PE/8)./AS);
end
for lll=1:length(scale_half)
%     recon_sig3(lll,:)=recon_sig(lll,:)./max(((recon_sig(lll,:))));
%     recon_sig(lll,:)=recon_sig3(lll,:).*recon_sig(lll,:);
end
% for lll=1:size(recon_sig,2)
%     recon_sig(:,lll)=ifft(ifftshift(recon_sig(:,lll)));
%     P = pseudoinverse(Z_sig_total(:,[met_ind]),[],'lsqr','tikhonov', {@(x,r) r*normest(Z_sig_total(:,[met_ind]))*x, 1e-8});
%     fit_param_total(:,lll)=P*(recon_sig(:,lll));
%     [fit_param_total(:,lll),phas]=phas_correct(fit_param_total(:,lll),size(fit_param_total(:,lll),1));
%     recon_sig(:,lll)=Z_sig_total(:,[met_ind])*fit_param_total(:,lll);
%     try
%         [fit_sig,amp,Ind_sig,freq] = HSVDPK_new(recon_sig(:,lll),fres,Fs,20) ;
%         f2=1;
%     catch 
%         f2=0;
%     end
%     if f2==1
%         ind=find((freq>=-50 & freq<=50));
%         half_sig2=sum(Ind_sig(:,ind),2);
%     else
%         half_sig2=0;
%     end
%     recon_sig(:,lll)=recon_sig(:,lll)-half_sig2;
%     recon_sig(:,lll)=fftshift(fft(recon_sig(:,lll)));
% end
% recon_sig=recon_sig./max(max(max(abs(recon_sig))));
% recon_sig_2=recon_sig_2.*(PE/AS);
% recon_sig=recon_sig.*recon_sig_2;
% recon_sig=recon_sig_2;
% recon_sig=(sum((recon_sig_total),3));
scale_half=(sum(ref_sig_total,3))';
mean_scale=sum(recon_sig,2)./(PE/2);
mean_scale_thresh=max(abs(fftshift(fft((sum(recon_sig,2)./((PE)*RE))))));
mean_scale_TD=mean_scale;%sum(mean_scale,2)./nreceivers;
mean_scale_TD(2:end)=0;
% recon_sig=recon_sig - (repmat(abs(recon_sig(:,1)),1,size(recon_sig,2)));
% recon_sig=(recon_sig) - (repmat(mean_scale_TD,1,size(recon_sig,2)));
recon_shift=mean_scale_TD(1);
% recon_sig=recon_sig + (mean_scale_TD(1));
% recon_sig=(recon_sig) - repmat(abs(mean_scale_TD),1,size(recon_sig,2));
for i=1:RE*PE
    recon_sig_fft=(fftshift(fft(recon_sig(:,i))));
    if round(Fs)==2000
        seg_NAA=160:180;
        seg_Lac=130:159;
        seg_Cho=198:210;
    elseif round(Fs)==1600
        seg_NAA=135:165;
        seg_Lac=110:130;
        seg_Cho=185:200;
    elseif round(Fs)==1000
        seg_NAA=60:95;
        seg_Lac=30:50;
        seg_Cho=145:165;
    else
        %     elseif round(Fs)==2500
        seg_NAA=186:200;
        seg_Lac=165:185;
        seg_Cho=210:225;
    end
    ind_base(i)=max(abs(recon_sig_fft(seg_Lac)));
end
ind_base_sorted=sort(ind_base);
div=mean(ind_base_sorted);
div_max=max(ind_base_sorted);
tmp_div = abs(ind_base_sorted-div);
[idx idx] = min(tmp_div) ;
ind_base_index=find(ind_base==ind_base_sorted(idx));
temp_sig=scale_half';%recon_sig(:,ind_base_index(1));
for i=1:RE*PE
%     recon_sig_fft=(fftshift(fft(recon_sig(:,i))));
    recon_sig_fft=(((recon_sig(:,i))));
    recon_sig_noise=(fftshift(fft(temp_sig)))./(PE*RE);
%     if max(abs(recon_sig_fft))<=3*max(abs(recon_sig_noise))
%         recon_sig_fft=abs(recon_sig_fft)-(((abs(recon_sig_noise))));%*max(abs(recon_sig_fft))./max(abs(recon_sig_noise)));
%     else
%         recon_sig_fft=abs(recon_sig_fft)-(abs(recon_sig_noise)./3);
%     end
    recon_sig(:,i)=ifft(ifftshift(recon_sig_fft));
%     ind_base=find(max(recon_sig_fft)<=6*mean_scale_thresh);
%     if length(ind_base)>=1
%         recon_sig(:,i)=0;
%     end
end    
% recon_sig(1:5,:)=0;
mean_scale_fft=abs(fftshift(fft(mean_scale)));
% figure();
for i=1:PE*RE
%     subplot(PE,RE,i)
    temp=(abs(fftshift(fft((recon_sig(:,i))))));
    recon_sig_FD(:,i)=((fftshift(fft((recon_sig(:,i))))));
%     plot(details_para.ppm_referenced(seg) ,temp(seg),'r');
%     axis([0.5 3.7 0 5]);
%     set(gca,'XDir','reverse');
%     ax=gca;
%     li=get(ax,'children');%spessore
%     set(li,'linewidth',1.5);
%     axis off;
end
% figure();
% for i=1:PE*RE*zerofill_value
%     subplot(PE*zerofill_value,RE,i)
%     temp=(abs(fftshift(fft( recon_sig_interp(:,i)))));
%     plot(temp(:));
%     set(gca,'XDir','reverse');
% end
data_value.FID_FD=recon_sig_FD;
data_value.FID_TD=recon_sig;
details_para.fres=fres;
details_para.Fs=Fs;
details_para.t=time_sample;
details_para.PE=PE;
details_para.RE=RE;
details_para.Tf = 127776603;
details_para.ref = 4.7;
