%% RTME-FMRI-PROCESSING
 
%% DEFINITIONS AND VARIABLE DECLARATIONS
 
% Initialize other path definitions
% cd '/Users/jheunis/Documents/MATLAB/';
data_dir = '/Users/sssssong/Documents/TUe/Internship/multi-echo-fmri-rtme';
spm_dir = '/Users/sssssong/Documents/TUe/Internship/spm12';
% mni_fn = '/Users/jheunis/Documents/MATLAB/ME/Volunteer_data/single_subj_T1.nii';
% mni_lowres_fn = '/Users/jheunis/Documents/MATLAB/ME/Volunteer_data/rsingle_subj_T1.nii';
 
% Initialize subject file names
 
    T1_img = cell(5,1);
    T2_img = cell(5,1);
    T3_img = cell(5,1);
    T1_img_s = cell(5,1);
    T2_img_s = cell(5,1);
    T3_img_s = cell(5,1);
    T1_img_s_HMP = cell(5,1);
    T2_img_s_HMP = cell(5,1);
    T3_img_s_HMP = cell(5,1);
    T4_img = cell(5,1);
    T5_img = cell(5,1);
    T6_img = cell(5,1);
    T7_img = cell(5,1);
    T8_img = cell(5,1);
    T9_img = cell(5,1);
 
    S_combined_img_rt= cell(5,1);
    S_combined_img = cell(5,1);
    S_combined_img_s = cell(5,1);
    S_combined_img_s_HMP = cell(5,1);
    S_combined_img2 = cell(5,1);
    S_combined_img_us_HMP = cell(5,1);
    S_combined_img_SNRweighted = cell(5,1);
 
    % smooth the combined echo
    S_combined_fdyn_img = cell(5,1);
    S_combined_img_3D = cell(5,1);
    S_combined_dyn_smoothed = cell(5,1);
    
    run = 3;
% for run = 1:5
    subj = 1;
    subj_dir = ['subj_' num2str(subj) '_' num2str(run+2)];
    s_fn = [data_dir filesep subj_dir filesep subj_dir '_T1W.nii'];
    f_me1_fn = [data_dir filesep subj_dir filesep '-e001-d0208.nii'];
    f_me2_fn = [data_dir filesep subj_dir filesep '-e002-d0208.nii'];
    f_me3_fn = [data_dir filesep subj_dir filesep '-e003-d0208.nii'];
    f_me_fn = {f_me1_fn, f_me2_fn, f_me3_fn};
    f_me_fnparts = {'-e001-d0208.nii','-e002-d0208.nii','-e003-d0208.nii'};
 
    % functional data variables
    Nt = 208; % number of volumes
    Ne = 3; % number of echos
    TE = [12 35 58]; % Echo times in ms
 
 
    %% PREPROCESSING
    disp('STEP 1: PREPROCESSING')
    % Preprocess structural and f0 images
    [d, f, e] = fileparts(s_fn);
    if exist([d filesep 'rc1' f e], 'file')
        disp('preproc already done, writing variable names')
        preproc_data = struct;
        preproc_data.forward_transformation = [d filesep 'y_' f e];
        preproc_data.inverse_transformation = [d filesep 'iy_' f e];
        preproc_data.gm_fn = [d filesep 'c1' f e];
        preproc_data.wm_fn = [d filesep 'c2' f e];
        preproc_data.csf_fn = [d filesep 'c3' f e];
        preproc_data.bone_fn = [d filesep 'c4' f e];
        preproc_data.soft_fn = [d filesep 'c5' f e];
        preproc_data.air_fn = [d filesep 'c6' f e];
        preproc_data.rstructural_fn = [d filesep 'r' f e];
        preproc_data.rgm_fn = [d filesep 'rc1' f e];
        preproc_data.rwm_fn = [d filesep 'rc2' f e];
        preproc_data.rcsf_fn = [d filesep 'rc3' f e];
        preproc_data.rbone_fn = [d filesep 'rc4' f e];
        preproc_data.rsoft_fn = [d filesep 'rc5' f e];
        preproc_data.rair_fn = [d filesep 'rc6' f e];
    else
        disp('doing preprocessing...')
        preproc_data = preRtPreProcME(f_me_fn, s_fn, spm_dir);
    end
    T1_spm = spm_vol(f_me1_fn);
    T2_spm = spm_vol(f_me2_fn);
    T3_spm = spm_vol(f_me3_fn);
    T1_img_test = spm_read_vols(T1_spm);
    T2_img_test = spm_read_vols(T2_spm);
    T3_img_test = spm_read_vols(T3_spm);
    T_img = {T1_img_test, T2_img_test, T3_img_test};
 
    %% Construct GM, WM and CSF masks (separately and together)
 
    % [GM_img_bin, WM_img_bin, CSF_img_bin] = create_binary_segments(gm_fn, wm_fn, csf_fn, threshold)
    [GM_img_bin, WM_img_bin, CSF_img_bin] = createBinarySegments(preproc_data.rgm_fn, preproc_data.rwm_fn, preproc_data.rcsf_fn, 0.1);
    I_GM = find(GM_img_bin);
    I_WM = find(WM_img_bin);
    I_CSF = find(CSF_img_bin);
    I_total_test = numel(I_GM)+numel(I_WM)+numel(I_CSF);
    mask_reshaped = GM_img_bin | WM_img_bin | CSF_img_bin;
    I_mask = find(mask_reshaped);
    test = 1:110592;
    I_rest = setdiff(test,I_mask);
    Nmaskvox = numel(I_mask); % total number of voxels in brain mask (GM+WM+CSF)
    Nvox = numel(GM_img_bin); % total number of voxels in whole image
    [Ni, Nj, Nk] = size(GM_img_bin);
 
    %% Data preparation/initialization
 
    % Predefine some matrices/structures
    ref_fn = [f_me_fn{2} ',1'];
    fdyn_fn = cell(Ne,1);
    currentVol = cell(Ne,1);
    F_dyn_img = cell(Ne,1);
    F_dyn = cell(Ne,1);
    F_dyn_resliced = cell(Ne,1);
    F_denoised = cell(Ne,1);
    F_dyn_denoised = cell(Ne,1);
    F_dyn_smoothed_HMP = cell(Ne,1);
    F_dyn_unsmoothed_HMP = cell(Ne,1);
    realign_params = cell(Ne,1);
    F_dyn_resliced_masked = cell(Ne,1);
    F_dyn_resliced_masked_img = cell(Ne,1);
    F_dyn_smoothed = cell(Ne,1);
    F_dyn_smoothed_masked = cell(Ne,1);
    F_dyn_smoothed_masked_img = cell(Ne,1);
    SNR = zeros(Nt,Ne);
    for e = 1:Ne
        F_dyn_resliced{e} = zeros(Ni*Nj*Nk, Nt);
        F_dyn_denoised{e} = zeros(Ni*Nj*Nk, Nt);
        F_dyn_smoothed_HMP{e} = zeros(Ni*Nj*Nk, Nt);
        F_dyn_unsmoothed_HMP{e} = zeros(Ni*Nj*Nk, Nt);
        F_dyn_resliced_masked{e} = zeros(Ni*Nj*Nk, Nt);
        F_dyn_resliced_masked_img{e} = zeros(Ni,Nj,Nk, Nt);
        F_dyn_smoothed{e} = zeros(Ni*Nj*Nk, Nt);
        F_dyn_smoothed_masked{e} = zeros(Ni*Nj*Nk, Nt);
        F_dyn_smoothed_masked_img{e} = zeros(Ni,Nj,Nk, Nt);
    end
 
    % Run SPM 1st level anlysis in order to get design matrix regressors
    sess_params = struct;
    sess_params.cond_name = 'name';
    sess_params.timing_units = 'scans';
    sess_params.timing_RT = 1.98;
    sess_params.cond_onset = [17;49;81;113;145;177];
    sess_params.cond_duration = [16;16;16;16;16;16];
    if ~exist([data_dir filesep subj_dir filesep 'SPM.mat'], 'file')
        spm_specify1stlevel_jsh([data_dir filesep subj_dir], f_me1_fn, '', sess_params)
    end
    load([data_dir filesep subj_dir filesep 'SPM.mat']);
    convolved_task_design = SPM.xX.X(:,1); % convolved task time course regressor
    drift_regressors = SPM.xX.K.X0; % cosine basis set for drift regressors
    X_design = [convolved_task_design drift_regressors ones(Nt,1)]; % design matrix, including task + drift + constat regressors
 
    % Load "resting state" HMP for real time calculation
    for k =1:Ne
        filename = ['/Users/sssssong/Documents/TUe/Internship/ME/subj_1_' num2str(run+2) '/rp_-e00' num2str(k) '-d0208.txt'];
        formatSpec = '%16f%16f%16f%16f%16f%f%[^\n\r]';
        fileID = fopen(filename,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string',  'ReturnOnError', false);
        fclose(fileID);
        rp_e = table(dataArray{1:end-1}, 'VariableNames', {'x','y','z','pitch','roll','yaw'});
        clearvars filename formatSpec fileID dataArray ans;
        % Data preparation 
        headmotion_regressors = table2cell(rp_e);
        headmotion_regressors = cell2mat(headmotion_regressors);
        X_design_HMP{k,1} = [convolved_task_design headmotion_regressors ones(Nt,1)];
    end
 
    %% Setup parameters for real-time realignment
 
    flagsSpmRealign = struct('quality',.9,'fwhm',5,'sep',4,...
        'interp',4,'wrap',[0 0 0],'rtm',0,'PW','','lkp',1:6);
    flagsSpmReslice = struct('quality',.9,'fwhm',5,'sep',4,...
        'interp',4,'wrap',[0 0 0],'mask',1,'mean',0,'which', 2);
    infoVolTempl = spm_vol(ref_fn);
    imgVolTempl  = spm_read_vols(infoVolTempl);
    dimTemplMotCorr     = infoVolTempl.dim;
    matTemplMotCorr     = infoVolTempl.mat;
    dicomInfoVox   = sqrt(sum(matTemplMotCorr(1:3,1:3).^2));
    nrSkipVol = 0;
    for e = 1:Ne
        realign_params{e}.A0=[];realign_params{e}.x1=[];realign_params{e}.x2=[];
        realign_params{e}.x3=[];realign_params{e}.wt=[];realign_params{e}.deg=[];realign_params{e}.b=[];
        realign_params{e}.R(1,1).mat = matTemplMotCorr;
        realign_params{e}.R(1,1).dim = dimTemplMotCorr;
        realign_params{e}.R(1,1).Vol = imgVolTempl;
    end
 
    %% SNR precalculation 
 
    for e = 1:Ne
        T_img_reshaped{e} = reshape(T_img{e}, Ni*Nj*Nk, Nt);
        for k = 1:Nt
            SNR(k,e) = mean(T_img_reshaped{e}(I_mask,k))/std(T_img_reshaped{e}(I_rest,k));
            
            % 这里改一下std的参数！！！！！！！！！！！
            % 这里改一下std的参数！！！！！！！！！！！
            % 这里改一下std的参数！！！！！！！！！！！
            % 这里改一下std的参数！！！！！！！！！！！
            % 这里改一下std的参数！！！！！！！！！！！
            % 这里改一下std的参数！！！！！！！！！！！
            % 这里改一下std的参数！！！！！！！！！！！
            % 原本是std没有0，2
            % mean也改了个,2    
        end
    end
 
    %% Initialize T2star and S0 estimation parameters
    X=[ones(3,1) -TE(:)]; %This means the estimated beta = [ln(S0); R2star];
    S0 = zeros(Nvox,Nt);
    T2star = zeros(Nvox,Nt);
    S0_u = zeros(Nvox,Nt);
    T2star_u = zeros(Nvox,Nt);
    S0_s = zeros(Nvox,Nt);
    T2star_s = zeros(Nvox,Nt);
    S0_u_HMP = zeros(Nvox,Nt);
    T2star_u_HMP = zeros(Nvox,Nt);
    S0_s_HMP = zeros(Nvox,Nt);
    T2star_s_HMP = zeros(Nvox,Nt);
    T2star_corrected = T2star;
    T2star_s_img = zeros(Ni, Nj, Nk, Nt);
    S0_img = zeros(Ni, Nj, Nk, Nt);
    S_combined_img{run} = zeros(Ni, Nj, Nk, Nt);
    S_combined_img2{run} = zeros(Ni, Nj, Nk, Nt);
    T = zeros(Nt,1);
 
    
    F_tsnr{run} = cell(14,1);
    F_tsnr_img{run} = cell(14,1);
    
    % Load the precalculated T2star 3D image
    load('T2star_PREmean_img');
    load('T2star_PREmean_img2');
    load('T2star_PRE_img');
 
 
 
    %% Create figure to display real-time images
 
    figure;
    slice = 10;
    ax1 = subplot(2,3,1); imagesc(squeeze(T2star_s_img(:,:,slice,1))); % just use T2star to create the subplots, this will be replaced in real-time by the appropriate figures
    ax2 = subplot(2,3,2); imagesc(squeeze(T2star_s_img(:,:,slice,1)));
    ax3 = subplot(2,3,3); imagesc(squeeze(T2star_s_img(:,:,slice,1)));
    ax4 = subplot(2,3,4); imagesc(squeeze(T2star_s_img(:,:,slice,1)));
    ax5 = subplot(2,3,5); imagesc(squeeze(T2star_s_img(:,:,slice,1)));
    ax6 = subplot(2,3,6); imagesc(squeeze(T2star_s_img(:,:,slice,1)));
 
 
    %% Real-Time Simulation
 
    for i = 1:Nt
        tic;
        % 0: Load multi-echo volumes for current iteration
        for e = 1:Ne
            fdyn_fn{e} = [f_me_fn{e} ',' num2str(i)]; % filename of dynamic functional image (fully unprocessed)
            currentVol{e} = spm_vol(fdyn_fn{e});
            F_dyn_img{e} = spm_read_vols(currentVol{e}); % this is the unprocessed image
            F_dyn{e}(:,i) = F_dyn_img{e}(:); % 2D
        end
 
        % 1: Realign, reslice and smooth each volume to the reference volume
        for e = 1:Ne
            realign_params{e}.R(2,1).mat = currentVol{e}.mat;
            realign_params{e}.R(2,1).dim = currentVol{e}.dim;
            realign_params{e}.R(2,1).Vol = F_dyn_img{e};
 
            % realign (FROM OPENNFT: preprVol.m)
            [realign_params{e}.R, realign_params{e}.A0, realign_params{e}.x1, realign_params{e}.x2, realign_params{e}.x3, realign_params{e}.wt, realign_params{e}.deg, realign_params{e}.b, realign_params{e}.nrIter] = spm_realign_rt(realign_params{e}.R, flagsSpmRealign, i, nrSkipVol + 1, realign_params{e}.A0, realign_params{e}.x1, realign_params{e}.x2, realign_params{e}.x3, realign_params{e}.wt, realign_params{e}.deg, realign_params{e}.b);
 
            % MC params (FROM OPENNFT: preprVol.m). STEPHAN NOTE: I don't understand this part, but it runs fine
            tmpMCParam = spm_imatrix(realign_params{e}.R(2,1).mat / realign_params{e}.R(1,1).mat);
            if (i == nrSkipVol + 1)
                realign_params{e}.offsetMCParam = tmpMCParam(1:6);
            end
            realign_params{e}.motCorrParam(i,:) = tmpMCParam(1:6)-realign_params{e}.offsetMCParam; % STEPHAN NOTE: I changed indVolNorm to indVol due to error, not sure if this okay or wrong?
            realign_params{e}.MP(i,:) = realign_params{e}.motCorrParam(i,:);
 
            % Reslice (FROM OPENNFT: preprVol.m)
            realign_params{e}.reslVol = spm_reslice_rt(realign_params{e}.R, flagsSpmReslice);
            F_dyn_resliced{e}(:,i) = realign_params{e}.reslVol(:);
 
            % Smooth (FROM OPENNFT: preprVol.m)
            s_fdyn_img = zeros(Ni, Nj, Nk); % this empty matrix is created for each volume, it gets updated by the smoothing operation
            gKernel = [8 8 8] ./ dicomInfoVox;
            spm_smooth(realign_params{e}.reslVol, s_fdyn_img, gKernel);
            F_dyn_smoothed{e}(:,i) = s_fdyn_img(:); % this is the realigned, resliced, and smoothed volume
 
            F_dyn_resliced_masked{e}(I_mask,i) = F_dyn_resliced{e}(I_mask,i);
            F_dyn_smoothed_masked{e}(I_mask,i) = F_dyn_smoothed{e}(I_mask,i);
        end
 
 
        % 2: GLM denoising
        for e = 1:Ne
            x_design = X_design(1:i, :); % only use the design matrix with the same amount of timepoints as are available until the current iteration
            beta_func = x_design\F_dyn_smoothed{e}(I_mask,1:i)'; % func = X*beta + e ==> beta = X\func ==> func_detrended = mp - X(i)*beta(i)
            F_denoised{e} = F_dyn_smoothed{e}(I_mask,1:i)' - x_design(:, 2:(end-1))*beta_func(2:(end-1), :); % remove effects of all regressors except constant and task
            F_denoised{e} = F_denoised{e}';
            F_dyn_denoised{e}(I_mask,i) = F_denoised{e}(:,i);
        end
 
        % Head motion correction
        for e = 1:Ne
            x_design_HMP{e} = X_design_HMP{e}(1:i, :);
            beta_func2{e} = x_design_HMP{e}\F_dyn_smoothed{e}(I_mask,1:i)';
            beta_func3{e} = x_design_HMP{e}\F_dyn_resliced{e}(I_mask,1:i)';
            F_smoothed_HMP{e} = F_dyn_smoothed{e}(I_mask,1:i)' - x_design_HMP{e}(:, 2:(end-1))*beta_func2{e}(2:(end-1), :); % remove effects of all regressors except constant and task
            F_smoothed_HMP{e} = F_smoothed_HMP{e}';
            F_unsmoothed_HMP{e} = F_dyn_resliced{e}(I_mask,1:i)' - x_design_HMP{e}(:, 2:(end-1))*beta_func3{e}(2:(end-1), :); % remove effects of all regressors except constant and task
            F_unsmoothed_HMP{e} = F_unsmoothed_HMP{e}';
            F_dyn_smoothed_HMP{e}(I_mask,i) = F_smoothed_HMP{e}(:,i);
            F_dyn_unsmoothed_HMP{e}(I_mask,i) = F_unsmoothed_HMP{e}(:,i);
        end
 
 
        % 3: Estimate T2star and S0 parameters
        % S = S0*exp(-t/T2star) = S0*exp(-t*R2star)
        % ln(S) = ln(S0*exp(-t/T2star)) = ln(S0) - t*R2star
        % [ln(S(TE1))]        = [ln(S0) - TE1*R2star]
        % [ln(S(TE2))]          [ln(S0) - TE2*R2star]
        % [ln(S(TE3))]          [ln(S0) - TE3*R2star]
        %                     = [1 - TE1]     [ln(S0)]
        %                       [1 - TE2]  *  [R2star]
        %                       [1 - TE3]
        % A = X*b ==> b = X\A;
        %â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”
        % unsmooth + no HMP T2* & S0 estimation
        clear S b x_design_HMP;
        S = [F_dyn_resliced{1}(I_mask,i)'; F_dyn_resliced{2}(I_mask,i)'; F_dyn_resliced{3}(I_mask,i)'];
        S = max(S,1e-11); %negative or zero signal values not possible
        b = X\log(S);
        S0_u(I_mask,i)=exp(b(1,:)); % 2D S0 estimation
        T2star_u(I_mask,i)=1./b(2,:); % 2D T2star estimation
        
        T2star_u_processed = zeros(Nvox,Nt);
        spm_smooth(T2star_u, T2star_u_processed, gKernel);
%         š„regressorçš„å€¼ï¼Ÿ
    %     switch S(I_mask,i)
    %         case F_dyn_resliced{1}(I_mask,i)
    %             x_design_HMP = X_design_HMP(1, i);
    %             beta = x_design_HMP\T2star_u_processed(I_mask,1:i)';
    %             T2star_u_processed = T2star_u_processed(I_mask,1:i)' - x_design_HMP(:, 2:(end-1))*beta(2:(end-1), :);
    %             T2star_u_processed = T2star_u_processed';
    %             T2star_u_processed(I_mask,i) = T2star_u_processed(:,i);
    %         case F_dyn_resliced{2}(I_mask,i)
    %             x_design_HMP = X_design_HMP(2, i);
    %             beta = x_design_HMP\T2star_u_processed(I_mask,1:i)';
    %             T2star_u_processed = T2star_u_processed(I_mask,1:i)' - x_design_HMP(:, 2:(end-1))*beta(2:(end-1), :);
    %             T2star_u_processed = T2star_u_processed';
    %             T2star_u_processed(I_mask,i) = T2star_u_processed(:,i);
    %         case F_dyn_resliced{3}(I_mask,i)
    %             x_design_HMP = X_design_HMP(3, i);
    %             beta = x_design_HMP\T2star_u_processed(I_mask,1:i)';
    %             T2star_u_processed = T2star_u_processed(I_mask,1:i)' - x_design_HMP(:, 2:(end-1))*beta(2:(end-1), :);
    %             T2star_u_processed = T2star_u_processed';
    %             T2star_u_processed(I_mask,i) = T2star_u_processed(:,i);
    %     end
 
        % smooth + no HMP T2* & S0 estimation
        clear S b;
        S = [F_dyn_smoothed_masked{1}(I_mask,i)'; F_dyn_smoothed_masked{2}(I_mask,i)'; F_dyn_smoothed_masked{3}(I_mask,i)'];
        S = max(S,1e-11); %negative or zero signal values not possible
        b = X\log(S);
        S0_s(I_mask,i)=exp(b(1,:)); % 2D S0 estimation
        T2star_s(I_mask,i)=1./b(2,:); % 2D T2star estimation
        T2star_s_img = reshape(T2star_s, Ni, Nj, Nk, Nt);
 
        % unsmooth + HMP T2* & S0 estimation
        clear S b;
        S = [F_dyn_unsmoothed_HMP{1}(I_mask,i)'; F_dyn_unsmoothed_HMP{2}(I_mask,i)'; F_dyn_unsmoothed_HMP{3}(I_mask,i)'];
        S = max(S,1e-11); %negative or zero signal values not possible
        b = X\log(S);
        S0_u_HMP(I_mask,i)=exp(b(1,:)); % 2D S0 estimation
        T2star_u_HMP(I_mask,i)=1./b(2,:); % 2D T2star estimation
        T2star_u_HMP_processed = zeros(Nvox,Nt);
        spm_smooth(T2star_u_HMP, T2star_u_HMP_processed, gKernel);
 
        % smooth + HMP T2* & S0 estimation
        clear S b;
        S = [F_dyn_smoothed_HMP{1}(I_mask,i)'; F_dyn_smoothed_HMP{2}(I_mask,i)'; F_dyn_smoothed_HMP{3}(I_mask,i)'];
        S = max(S,1e-11); %negative or zero signal values not possible
        b = X\log(S);
        S0_s_HMP(I_mask,i)=exp(b(1,:)); % 2D S0 estimation
        T2star_s_HMP(I_mask,i)=1./b(2,:); % 2D T2star estimation
        %â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”
 
        % check if any values are Not-A-Number
        if isnan(b)
            disp(['isnan: i = ' i ])
        end
        % Threshold T2star values 
        T2s_max = 100;
        T2star_corrected(I_mask,i) = T2star_s(I_mask,i);
        T2star_corrected((T2star_corrected(:,i)<0), i) = 0;
        T2star_corrected((T2star_corrected(:,i)>T2s_max), i) = 0;
 
        % Reshape estimated parameters to 3D volumes
        T2star_s_img(:,:,:,i) = reshape(T2star_corrected(:,i), Ni, Nj, Nk);
        S0_img(:,:,:,i) = reshape(S0(:,i), Ni, Nj, Nk);
 
 
        % 4: Multi-echo combination (we can change which images we would like to combine and display in the figures)
        T1_img{run} = reshape(F_dyn_resliced_masked{1}(:,i), Ni, Nj, Nk); 
        T2_img{run} = reshape(F_dyn_resliced_masked{2}(:,i), Ni, Nj, Nk);
        T3_img{run} = reshape(F_dyn_resliced_masked{3}(:,i), Ni, Nj, Nk);
        T1_img_s{run} = reshape(F_dyn_smoothed_masked{1}(:,i), Ni, Nj, Nk);
        T2_img_s{run} = reshape(F_dyn_smoothed_masked{2}(:,i), Ni, Nj, Nk);
        T3_img_s{run} = reshape(F_dyn_smoothed_masked{3}(:,i), Ni, Nj, Nk);
        T1_img_s_HMP{run} = reshape(F_dyn_smoothed_HMP{1}(:,i), Ni, Nj, Nk);
        T2_img_s_HMP{run} = reshape(F_dyn_smoothed_HMP{2}(:,i), Ni, Nj, Nk);
        T3_img_s_HMP{run} = reshape(F_dyn_smoothed_HMP{3}(:,i), Ni, Nj, Nk);
        T4_img{run} = reshape(F_dyn_denoised{1}(:,i), Ni, Nj, Nk);
        T5_img{run} = reshape(F_dyn_denoised{2}(:,i), Ni, Nj, Nk);
        T6_img{run} = reshape(F_dyn_denoised{3}(:,i), Ni, Nj, Nk);
        T7_img{run} = reshape(F_dyn_unsmoothed_HMP{1}(:,i), Ni, Nj, Nk);
        T8_img{run} = reshape(F_dyn_unsmoothed_HMP{2}(:,i), Ni, Nj, Nk);
        T9_img{run} = reshape(F_dyn_unsmoothed_HMP{3}(:,i), Ni, Nj, Nk);
 
        % Use real-time T2star for combination:
        S_combined_img_rt{run}(:,:,:,i) = combine_echoes(TE, T2star_s_img(:,:,:,i), I_mask, 1, T1_img{run}, T2_img{run}, T3_img{run});
    %     S_combined_img2(:,:,:,i) = combine_echoes(TE, T2star_img(:,:,:,i), I_mask, 1, T4_img, T5_img, T6_img);
 
        % Use pre-calculated T2star for combination
        S_combined_img{run}(:,:,:,i) = combine_echoes(TE, T2star_PREmean_img2, I_mask, 1, T1_img{run}, T2_img{run}, T3_img{run}); 
    %     S_combined_img_test(:,:,:,i) = combine_echoes(TE, T2star_PRE_img, I_mask, 3, T1_img_test, T2_img_test, T3_img_test); 
        S_combined_img_s{run}(:,:,:,i) = combine_echoes(TE, T2star_PREmean_img2, I_mask, 1, T1_img_s{run}, T2_img_s{run}, T3_img_s{run}); 
        S_combined_img_s_HMP{run}(:,:,:,i) = combine_echoes(TE, T2star_PREmean_img2, I_mask, 1, T1_img_s_HMP{run}, T2_img_s_HMP{run}, T3_img_s_HMP{run}); 
        S_combined_img2{run}(:,:,:,i) = combine_echoes(TE, T2star_PREmean_img2, I_mask, 1, T4_img{run}, T5_img{run}, T6_img{run});
        S_combined_img_us_HMP{run}(:,:,:,i) = combine_echoes(TE, T2star_PREmean_img2, I_mask, 1, T7_img{run}, T8_img{run}, T9_img{run});
        S_combined_img_SNRweighted{run}(:,:,:,i) = combine_echoes(TE, SNR(i,:), I_mask, 4, T1_img{run}, T2_img{run}, T3_img{run});
 
        % smooth the combined echo
        S_combined_fdyn_img{run} = zeros(Ni, Nj, Nk);
        S_combined_img_3D{run} = reshape(S_combined_img{run}(:,:,:,i),Ni, Nj, Nk);
        spm_smooth(S_combined_img_3D{run}, S_combined_fdyn_img{run}, gKernel);
        S_combined_dyn_smoothed{run}(:,i) = S_combined_fdyn_img{run}(:); % this is the realigned, resliced, and smoothed volume
 
 
        % 5: Update figure
        % 2 x 3 subplots in order: Echo 1, Echo 2, Echo 3, combined 1, T2star,
        % combined 2
        imagesc(ax1, squeeze(T1_img{run}(:,:,slice))); colormap(ax1, 'bone'); colorbar(ax1); ax1.CLim = [0 15e5]; 
        imagesc(ax2, squeeze(T2_img{run}(:,:,slice))); colormap(ax2, 'bone'); colorbar(ax2); ax2.CLim = [0 15e5]; 
        imagesc(ax3, squeeze(T3_img{run}(:,:,slice))); colormap(ax3, 'bone'); colorbar(ax3); ax3.CLim = [0 15e5]; 
        imagesc(ax4, squeeze(S_combined_img{run}(:,:,slice,i))); colormap(ax4, 'bone'); colorbar(ax4); ax4.CLim = [0 15e5]; 
        imagesc(ax5, squeeze(T2star_s_img(:,:,slice,i))); colormap(ax5, 'hot'); colorbar(ax5); ax5.CLim = [0 300]; 
        imagesc(ax6, squeeze(S_combined_img2{run}(:,:,slice,i))); colormap(ax6, 'bone'); colorbar(ax6); ax6.CLim = [0 15e5]; 
        drawnow;
        T(i) = toc;
        disp(['i=' num2str(i) ': ' num2str(T(i))]); % display iteration time
 
    end
 
 
    %% generate and plot tSNR images (THIS SHOULD CHANGE BASED ON WHAT WE WANT TO TEST) 
    F_tsnr{run} = cell(Ne+3,1);
    F_tsnr_img{run} = cell(Ne+3,1);
 
    % for e = 1:Ne
    %     clear m stddev;
    %     m = mean(F_dyn_smoothed_masked{e}(I_mask,:), 2);
    %     stddev = std(F_dyn_smoothed_masked{e}(I_mask,:), 0, 2);
    %     F_tsnr{e} = zeros(Nvox,1);
    %     F_tsnr{e}(I_mask) = m./stddev;
    %     F_tsnr{e}(isnan(F_tsnr{e}))=0;
    %     F_tsnr_img{e} = reshape(F_tsnr{e}, Ni, Nj, Nk);
    %     
    % end
    %â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”
    % middle echo vs. combined echo(with different combination method)

    % tSNR - middle echo - F_tsnr_img{1}
    clear m stddev;
    m = mean(F_dyn_resliced_masked{2}(I_mask,:), 2);
    stddev = std(F_dyn_resliced_masked{2}(I_mask,:), 0, 2);
    F_tsnr{run}{1} = zeros(Nvox,1);
    F_tsnr{run}{1}(I_mask) = m./stddev;
    F_tsnr{run}{1}(isnan(F_tsnr{run}{1}))=0;
    F_tsnr_img{run}{1} = reshape(F_tsnr{run}{1}, Ni, Nj, Nk);
 
    % tSNR - combined echo (pre-calculated T2*) - F_tsnr_img{2}
    clear m stddev combined;
    combined = reshape(S_combined_img{run}, Nvox, Nt);
    combined(isnan(combined))=0;
    m = mean(combined(I_mask,:), 2);
    stddev = std(combined(I_mask,:), 0, 2);
    F_tsnr{run}{2} = zeros(Nvox,1);
    F_tsnr{run}{2}(I_mask) = m./stddev;
    F_tsnr{run}{2}(isnan(F_tsnr{run}{2}))=0;
    F_tsnr_img{run}{2} = reshape(F_tsnr{run}{2}, Ni, Nj, Nk);
 
    % tSNR - combined echo (real time T2*) - F_tsnr_img{3}
    clear m stddev combined;
    combined = reshape(S_combined_img_rt{run}, Nvox, Nt);
    combined(isnan(combined))=0;
    m = mean(combined(I_mask,:), 2);
    stddev = std(combined(I_mask,:), 0, 2);
    F_tsnr{run}{3} = zeros(Nvox,1);
    F_tsnr{run}{3}(I_mask) = m./stddev;
    F_tsnr{run}{3}(isnan(F_tsnr{run}{3}))=0;
    F_tsnr_img{run}{3} = reshape(F_tsnr{run}{3}, Ni, Nj, Nk);
 
    % tSNR - combined echo (pre-calculated SNR) - F_tsnr_img{4}
    clear m stddev combined;
    combined = reshape(S_combined_img_SNRweighted{run}, Nvox, Nt);
    combined(isnan(combined))=0;
    m = mean(combined(I_mask,:), 2);
    stddev = std(combined(I_mask,:), 0, 2);
    F_tsnr{run}{4} = zeros(Nvox,1);
    F_tsnr{run}{4}(I_mask) = m./stddev;
    F_tsnr{run}{4}(isnan(F_tsnr{run}{4}))=0;
    F_tsnr_img{run}{4} = reshape(F_tsnr{run}{4}, Ni, Nj, Nk);
    % %â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”
    % Influence of order of smoothing & HMP operations on tSNR of T2* &
    % combined echo (estimated in real time)
 
    % _____combined echo_____
    % unsmoothed + no HMP - F_tsnr_img{5}
    clear m stddev combined;
    combined = reshape(S_combined_img{run}, Nvox, Nt);
    combined(isnan(combined))=0;
    m = mean(combined(I_mask,:), 2);
    stddev = std(combined(I_mask,:), 0, 2);
    F_tsnr{run}{5} = zeros(Nvox,1);
    F_tsnr{run}{5}(I_mask) = m./stddev;
    F_tsnr{run}{5}(isnan(F_tsnr{run}{5}))=0;
    F_tsnr_img{run}{5} = reshape(F_tsnr{run}{5}, Ni, Nj, Nk);
 
 
    % smoothed + no HMP - F_tsnr_img{6}
    clear m stddev combined;
    combined = reshape(S_combined_img_s{run}, Nvox, Nt);
    m = mean(combined(I_mask,:), 2);
    stddev = std(combined(I_mask,:), 0, 2);
    F_tsnr{run}{6} = zeros(Nvox,1);
    F_tsnr{run}{6}(I_mask) = m./stddev;
    F_tsnr{run}{6}(isnan(F_tsnr{run}{6}))=0;
    F_tsnr_img{run}{6} = reshape(F_tsnr{run}{6}, Ni, Nj, Nk);
 
 
    % unsmoothed + HMP - F_tsnr_img{7}
    clear m stddev combined;
    combined = reshape(S_combined_img_us_HMP{run}, Nvox, Nt);
    combined(isnan(combined))=0;
    m = mean(combined(I_mask,:), 2);
    stddev = std(combined(I_mask,:), 0, 2);
    F_tsnr{run}{7} = zeros(Nvox,1);
    F_tsnr{run}{7}(I_mask) = m./stddev;
    F_tsnr{run}{7}(isnan(F_tsnr{run}{7}))=0;
    F_tsnr_img{run}{7} = reshape(F_tsnr{run}{7}, Ni, Nj, Nk);
 
    % smoothed + HMP - F_tsnr_img{8}
    clear m stddev combined;
    combined = reshape(S_combined_img_s_HMP{run}, Nvox, Nt);
    combined(isnan(combined))=0;
    m = mean(combined(I_mask,:), 2);
    stddev = std(combined(I_mask,:), 0, 2);
    F_tsnr{run}{8} = zeros(Nvox,1);
    F_tsnr{run}{8}(I_mask) = m./stddev;
    F_tsnr{run}{8}(isnan(F_tsnr{run}{8}))=0;
    F_tsnr_img{run}{8} = reshape(F_tsnr{run}{8}, Ni, Nj, Nk);
 
    % _____T2*_____
    % unsmoothed + no HMP - F_tsnr_img{9}
    clear m stddev combined;
    m = mean(T2star_u_processed(I_mask,:), 2);
    stddev = std(T2star_u_processed(I_mask,:), 0, 2);
    F_tsnr{run}{9} = zeros(Nvox,1);
    F_tsnr{run}{9}(I_mask) = m./stddev;
    F_tsnr{run}{9}(isnan(F_tsnr{run}{9}))=0;
    F_tsnr_img{run}{9} = reshape(F_tsnr{run}{9}, Ni, Nj, Nk);
 
    % smoothed + no HMP - F_tsnr_img{10}
    clear m stddev combined;
    % m = mean(T2star_s_processed(I_mask,:), 2);
    % stddev = std(T2star_s_processed(I_mask,:), 0, 2);
    m = mean(T2star_s(I_mask,:), 2);
    stddev = std(T2star_s(I_mask,:), 0, 2);
    F_tsnr{run}{10} = zeros(Nvox,1);
    F_tsnr{run}{10}(I_mask) = m./stddev;
    F_tsnr{run}{10}(isnan(F_tsnr{run}{10}))=0;
    F_tsnr_img{run}{10} = reshape(F_tsnr{run}{10}, Ni, Nj, Nk);
 
    % unsmoothed + HMP - F_tsnr_img{11}
    clear m stddev combined;
    m = mean(T2star_u_HMP_processed(I_mask,:), 2);
    stddev = std(T2star_u_HMP_processed(I_mask,:), 0, 2);
    F_tsnr{run}{11} = zeros(Nvox,1);
    F_tsnr{run}{11}(I_mask) = m./stddev;
    F_tsnr{run}{11}(isnan(F_tsnr{run}{11}))=0;
    F_tsnr_img{run}{11} = reshape(F_tsnr{run}{11}, Ni, Nj, Nk);
 
    % smoothed + HMP - F_tsnr_img{12}
    clear m stddev combined;
    m = mean(T2star_s_HMP(I_mask,:), 2);
    stddev = std(T2star_s_HMP(I_mask,:), 0, 2);
    F_tsnr{run}{12} = zeros(Nvox,1);
    F_tsnr{run}{12}(I_mask) = m./stddev;
    F_tsnr{run}{12}(isnan(F_tsnr{run}{12}))=0;
    F_tsnr_img{run}{12} = reshape(F_tsnr{run}{12}, Ni, Nj, Nk);
    %â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”
    % SMOOTHED BEFORE - F_tsnr_img{13}
    clear m stddev combined;
    combined = reshape(S_combined_img_s{run}, Nvox, Nt);
    m = mean(combined(I_mask,:), 2);
    stddev = std(combined(I_mask,:), 0, 2);
    F_tsnr{run}{13} = zeros(Nvox,1);
    F_tsnr{run}{13}(I_mask) = m./stddev;
    F_tsnr{run}{13}(isnan(F_tsnr{run}{13}))=0;
    F_tsnr_img{run}{13} = reshape(F_tsnr{run}{13}, Ni, Nj, Nk);
 
    % SMOOTHED AFTER - F_tsnr_img{14}
    clear m stddev combined;
    combined = S_combined_dyn_smoothed{run};
    combined(isnan(combined))=0;
    m = mean(combined(I_mask,:), 2);
    stddev = std(combined(I_mask,:), 0, 2);
    F_tsnr{run}{14} = zeros(Nvox,1);
    F_tsnr{run}{14}(I_mask) = m./stddev;
    F_tsnr{run}{14}(isnan(F_tsnr{run}{14}))=0;
    F_tsnr_img{run}{14} = reshape(F_tsnr{run}{14}, Ni, Nj, Nk);
 
    %â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”â€”
 
    %% plot tSNR images 
    run = 3;
    figure;
    slice = 10;
    ax10 = subplot(2,2,1); imagesc(ax10, squeeze(F_tsnr_img{run}{1}(:,:,slice))); colormap(ax10, 'hot'); colorbar(ax10); ax10.CLim = [0 400]; title('tSNR - middle echo');
    ax11 = subplot(2,2,2); imagesc(ax11, squeeze(F_tsnr_img{run}{2}(:,:,slice))); colormap(ax11, 'hot'); colorbar(ax11); ax11.CLim = [0 400]; title('tSNR - combined echo (pre-calculated T2*)');
    ax12 = subplot(2,2,3); imagesc(ax12, squeeze(F_tsnr_img{run}{3}(:,:,slice))); colormap(ax12, 'hot'); colorbar(ax12); ax12.CLim = [0 400]; title('tSNR - combined echo (real time T2*)');
    ax13 = subplot(2,2,4); imagesc(ax13, squeeze(F_tsnr_img{run}{4}(:,:,slice))); colormap(ax13, 'hot'); colorbar(ax13); ax13.CLim = [0 400]; title('tSNR - combined echo (pre-calculated SNR)');
 
    figure;
    slice = 10;
    ax21 = subplot(2,2,1); imagesc(ax21, squeeze(F_tsnr_img{run}{5}(:,:,slice))); colormap(ax21, 'hot'); colorbar(ax21); ax21.CLim = [0 400]; title('a','fontsize',14);
%     title({'unsmoothed + without HMP - tSNR for Combined echo'});
    ax22 = subplot(2,2,2); imagesc(ax22, squeeze(F_tsnr_img{run}{6}(:,:,slice))); colormap(ax22, 'hot'); colorbar(ax22); ax22.CLim = [0 400]; title('b','fontsize',14);
%     title({'smoothed + without HMP - tSNR for Combined echo'});
    ax23 = subplot(2,2,3); imagesc(ax23, squeeze(F_tsnr_img{run}{7}(:,:,slice))); colormap(ax23, 'hot'); colorbar(ax23); ax23.CLim = [0 400]; title('c','fontsize',14);
%     title({'unsmoothed + with HMP - tSNR for Combined echo'});
    ax24 = subplot(2,2,4); imagesc(ax24, squeeze(F_tsnr_img{run}{8}(:,:,slice))); colormap(ax24, 'hot'); colorbar(ax24); ax24.CLim = [0 400]; title('d','fontsize',14);
%     title({'smoothed + with HMP - tSNR for Combined echo'});
 
    figure;
    slice = 10;
    ax31 = subplot(2,2,1); imagesc(ax31, squeeze(F_tsnr_img{run}{9}(:,:,slice))); colormap(ax31, 'hot'); colorbar(ax31); ax31.CLim = [0 400]; title({'unsmoothed + without HMP - tSNR for T2*'});
    ax32 = subplot(2,2,2); imagesc(ax32, squeeze(F_tsnr_img{run}{10}(:,:,slice))); colormap(ax32, 'hot'); colorbar(ax32); ax32.CLim = [0 400]; title({'smoothed + without HMP - tSNR for T2*'});
    ax33 = subplot(2,2,3); imagesc(ax33, squeeze(F_tsnr_img{run}{11}(:,:,slice))); colormap(ax33, 'hot'); colorbar(ax33); ax33.CLim = [0 400]; title({'unsmoothed + with HMP - tSNR for T2*'});
    ax34 = subplot(2,2,4); imagesc(ax34, squeeze(F_tsnr_img{run}{12}(:,:,slice))); colormap(ax34, 'hot'); colorbar(ax34); ax34.CLim = [0 400]; title({'smoothed + with HMP - tSNR for T2*'});
 
    figure;
    slice = 10;
    ax41 = subplot(1,2,1); imagesc(ax41, squeeze(F_tsnr_img{run}{13}(:,:,slice))); colormap(ax41, 'hot'); colorbar(ax41); ax41.CLim = [0 400]; title({'Smooth each single echo and then combination'});
    ax42 = subplot(1,2,2); imagesc(ax42, squeeze(F_tsnr_img{run}{14}(:,:,slice))); colormap(ax42, 'hot'); colorbar(ax42); ax42.CLim = [0 400]; title({'Only smooth the conbined echo'});
 
    % A. Middle echo vs combined echo (4 images in total)
% end
 
%% Combined method comparison
F_tsnr_img_average = cell(4,1);
sum = cell(4,1);
comparision{k} = cell(3,1);
for i = 1:4
    sum{i} = zeros(Ni,Nj,Nk);
    for run = 1:5
        sum{i} = sum{i} + F_tsnr_img{run}{i};
    end
end
 
for j = 1:4
    F_tsnr_img_average{j} = sum{j}./5;
end
 
for k = 1:3
    comparision{k} = zeros(Ni,Nj,Nk);
end
 
comparision{1} =  F_tsnr_img_average{2} - F_tsnr_img_average{1};
comparision{2} =  F_tsnr_img_average{3} - F_tsnr_img_average{1};
comparision{3} =  F_tsnr_img_average{4} - F_tsnr_img_average{1};
 
figure;
slice = 10;
% ax51 = subplot(2,2,1); imagesc(ax51, squeeze(F_tsnr_img_average{1}(:,:,slice))); colormap(ax51, 'hot'); colorbar(ax51); ax51.CLim = [0 400]; title({'unsmoothed + without HMP - tSNR for Combined echo'});
% 'Difference between Combined echo (using '  'pre-calculated T2*) and Middel echo' 
ax52 = subplot(1,3,1); imagesc(ax52, squeeze(comparision{1}(:,:,slice))); colormap(ax52, 'hot'); colorbar(ax52); ax52.CLim = [0 50]; title('a','fontsize',14);
% Difference between Combined echo (using '  'real time T2*) and Middel echo
ax53 = subplot(1,3,2); imagesc(ax53, squeeze(comparision{2}(:,:,slice))); colormap(ax53, 'hot'); colorbar(ax53); ax53.CLim = [0 50]; title('b','fontsize',14);
% Difference between Combined echo (using '  'pre-calculated SNR) and Middel echo
ax54 = subplot(1,3,3); imagesc(ax54, squeeze(comparision{3}(:,:,slice))); colormap(ax54, 'hot'); colorbar(ax54); ax54.CLim = [0 50]; title('c','fontsize',14);
 
 
    %%
    montage1 = createMontage(F_tsnr_img{3}{1}, 4, 1, 'a) tSNR (realigned echo 2)', 'hot', 0, [0 400]);
    montage2 = createMontage(F_tsnr_img{3}{2}, 4, 1, 'b) tSNR (realigned combined - pre-calculated T2*)', 'hot', 0, [0 400]);
    montage3 = createMontage(F_tsnr_img{3}{3}, 4, 1, 'c) tSNR (realigned combined - real time T2*)', 'hot', 0, [0 400]);
    montage4 = createMontage(F_tsnr_img{3}{4}, 4, 1, 'd) tSNR (realigned combined - pre-calculated SNR)', 'hot', 0, [0 400]);
    % tSNR (realigned combined - smooth before)
    montage5 = createMontage(F_tsnr_img{3}{13}, 4, 1, 'a', 'hot', 0, [0 400]);
    % tSNR (realigned combined - smmoth after)
    montage6 = createMontage(F_tsnr_img{3}{14}, 4, 1, 'b', 'hot', 0, [0 400]);
 %%
    % unsmoothed + without HMP - tSNR for Combined echo
    montage7 = createMontage(F_tsnr_img{3}{5}, 4, 1, 'a', 'hot', 0, [0 300]);
    % smoothed + without HMP - tSNR for Combined echo
    montage8 = createMontage(F_tsnr_img{3}{6}, 4, 1, 'b', 'hot', 0, [0 300]);
    % unsmoothed + with HMP - tSNR for Combined echo
    montage9 = createMontage(F_tsnr_img{3}{7}, 4, 1, 'c', 'hot', 0, [0 300]);
    % smoothed + with HMP - tSNR for Combined echo
    montage10 = createMontage(F_tsnr_img{3}{8}, 4, 1, 'd', 'hot', 0, [0 300]);

    
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
%     %% ________________________________________________________________________
% 
%     % IGNORE THE REST BELOW THIS LINE
% 
% 
%     %% ________________________________________________________________________
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
%     %% ________________________________________________________________________
%     % %%
%     % % image = S_preestimation_denoised{2};
%     % % S_combined_2D = reshape(S_combined_img, Ni*Nj*Nk, Nt);
%     % % image = S_combined_2D;
%     % T2star_2D = reshape(T2star_img, Ni*Nj*Nk, Nt);
%     % image = T2star_2D;
%     % 
%     % % Statistical measures
%     % F2D_mean = mean(image, 2);
%     % F2D_stddev = std(image, 0, 2);
%     % F2D_zstat = (image - F2D_mean)./F2D_stddev;
%     % F2D_zstat(isnan(F2D_zstat))=0;
%     % F2D_zstat_mean = mean(abs(F2D_zstat),1);
%     % Zstat_mean = mean(F2D_zstat_mean);
%     % F2D_psc = 100*(image./repmat(F2D_mean, 1, Nt)) - 100;
%     % F2D_psc(isnan(F2D_psc))=0;
%     % F2D_diff = [zeros(1, Ni*Nj*Nk); diff(image')]';
%     % F2D_DVARS = var(F2D_diff);
%     % % tSNR
%     % tSNR_2D = F2D_mean./F2D_stddev;
%     % tSNR_2D(isnan(tSNR_2D))=0;
%     % tSNR_brain = mean(tSNR_2D(I_mask));
%     % tSNR_GM = mean(tSNR_2D(I_GM));
%     % tSNR_WM = mean(tSNR_2D(I_WM));
%     % tSNR_CSF = mean(tSNR_2D(I_CSF));
%     % % Metrics
%     % disp(['Mean Zscore: ' num2str(Zstat_mean)])
%     % disp(['tSNR (brain): ' num2str(tSNR_brain)])
%     % disp(['tSNR (GM): ' num2str(tSNR_GM)])
%     % disp(['tSNR (WM): ' num2str(tSNR_WM)])
%     % disp(['tSNR (CSF): ' num2str(tSNR_CSF)])
%     % 
%     % % 3D and 4D images
%     % mask_3D = reshape(mask_reshaped, Ni, Nj, Nk);
%     % tSNR_3D = reshape(tSNR_2D, Ni, Nj, Nk);
%     % F3D_mean = reshape(F2D_mean, Ni, Nj, Nk);
%     % F3D_stddev = reshape(F2D_stddev, Ni, Nj, Nk);
%     % 
%     % montage2 = createMontage(F3D_mean, 5, 1, 'Mean EPI (whole image)', 'gray');
%     % montage3 = createMontage(F3D_stddev, 5, 1, 'Standard deviation (whole image)', 'parula');
%     % montage1 = createMontage(tSNR_3D, 5, 1, 'tSNR (whole image)', 'hot');
% 
%     %% PSC
%     % echo2_resliced = reshape(F_dyn_resliced_masked{2}, Ni, Nj, Nk, Nt);
%     % echo2_denoised = reshape(F_dyn_denoised{2}, Ni, Nj, Nk, Nt);
%     echo2_resliced = F_dyn_smoothed_masked{2};
%     echo2_denoised = F_dyn_denoised{2};
%     PSC_echo2r = getPercentageSignalChange(echo2_resliced, 6, 7, 16, [17;49;81;113;145;177], 208);
%     PSC_echo2d = getPercentageSignalChange(echo2_denoised, 6, 7, 16, [17;49;81;113;145;177], 208);
% 
%     S1 = reshape(S_combined_img, Nvox, Nt);
%     S2 = reshape(S_combined_img, Nvox, Nt);
%     PSC_Sr = getPercentageSignalChange(S1, 6, 7, 16, [17;49;81;113;145;177], 208);
%     PSC_Sd = getPercentageSignalChange(S2, 6, 7, 16, [17;49;81;113;145;177], 208);
% 
%     PSC_blocks_maskedr = zeros(Ni*Nj*Nk, Nt);
%     PSC_blocks_maskedr(I_mask, :) = PSC_echo2r.PSC_blocks(I_mask, :);
%     PSC_blocks_4Dr = reshape(PSC_blocks_maskedr, Ni, Nj, Nk, Nt);
% 
%     PSC_blocks_maskedd = zeros(Ni*Nj*Nk, Nt);
%     PSC_blocks_maskedd(I_mask, :) = PSC_echo2d.PSC_blocks(I_mask, :);
%     PSC_blocks_4Dd = reshape(PSC_blocks_maskedd, Ni, Nj, Nk, Nt);
% 
%     PSC_blocks_maskedSr = zeros(Ni*Nj*Nk, Nt);
%     PSC_blocks_maskedSr(I_mask, :) = PSC_Sr.PSC_blocks(I_mask, :);
%     PSC_blocks_4DSr = reshape(PSC_blocks_maskedSr, Ni, Nj, Nk, Nt);
% 
%     PSC_blocks_maskedSd = zeros(Ni*Nj*Nk, Nt);
%     PSC_blocks_maskedSd(I_mask, :) = PSC_Sd.PSC_blocks(I_mask, :);
%     PSC_blocks_4DSd = reshape(PSC_blocks_maskedSd, Ni, Nj, Nk, Nt);
% 
%     temporal_mask = ((PSC_echo2d.task_blocks==5)|(PSC_echo2d.task_blocks==6));
% 
%     PSC_ave_resliced = mean(PSC_blocks_maskedr(:,temporal_mask), 2);
%     PSC_ave_resliced3D = reshape(PSC_ave_resliced, Ni, Nj, Nk);
%     PSC_ave_denoised = mean(PSC_blocks_maskedd(:,temporal_mask), 2);
%     PSC_ave_denoised3D = reshape(PSC_ave_denoised, Ni, Nj, Nk);
% 
%     PSC_Save_resliced = mean(PSC_blocks_maskedSr(:,temporal_mask), 2);
%     PSC_Save_resliced3D = reshape(PSC_Save_resliced, Ni, Nj, Nk);
%     PSC_Save_denoised = mean(PSC_blocks_maskedSd(:,temporal_mask), 2);
%     PSC_Save_denoised3D = reshape(PSC_Save_denoised, Ni, Nj, Nk);
% 
%     figure;
%     slice = 10;
%     ax21 = subplot(2,2,1); imagesc(ax21, squeeze(PSC_ave_resliced3D(:,:,slice))); colormap(ax21, 'hot'); colorbar(ax21); ax21.CLim = [0 10]; 
%     ax22 = subplot(2,2,2); imagesc(ax22, squeeze(PSC_ave_denoised3D(:,:,slice))); colormap(ax22, 'hot'); colorbar(ax22); ax22.CLim = [0 10]; 
%     ax23 = subplot(2,2,3); imagesc(ax23, squeeze(PSC_Save_resliced3D(:,:,slice))); colormap(ax23, 'hot'); colorbar(ax23); ax23.CLim = [0 10]; 
%     ax24 = subplot(2,2,4); imagesc(ax24, squeeze(PSC_Save_denoised3D(:,:,slice))); colormap(ax24, 'hot'); colorbar(ax24); ax24.CLim = [0 10]; 
% 
% 
% 
%     %%
%     montage3 = createMontage(PSC_ave_denoised3D, 4, 1, 'PSC (denoised echo 2)', 'hot', 0, [-5 5]);
%     montage4 = createMontage(PSC_Save_denoised3D, 4, 1, 'PSC (denoised combined)', 'hot', 0, [-5 5]);
% 
% 
% 
% 
% % PSC_all_masked = zeros(Ni*Nj*Nk, Nt);
% % PSC_all_masked(I_mask, :) = F2D_psc(I_mask, :);
% % PSC_all_4D = reshape(PSC_all_masked, Ni, Nj, Nk, Nt);
% 
% % new_nii = make_nii(PSC_all_4D, [3.5 3.5 4.5], [27 37 17]);
% % save_nii(new_nii, 'T2S_PSC_all.nii')
% % new_nii = make_nii(PSC_blocks_4D, [3.5 3.5 4.5], [27 37 17]);
% % save_nii(new_nii, 'T2S_PSC_blocks.nii')
% % 
% % draw_brain_views_gui([20 20 20], PSC_blocks_4D)
% % draw_brain_views_gui([20 20 20], PSC_all_4D)
% 
% 
% 
% %% cluster analysis
% stats_dir = [data_dir filesep subj_dir filesep 'stats'];
% tmap_fn = [stats_dir filesep 'spmT_0001.nii'];
% xSPM = load([stats_dir filesep 'xSPM.mat']);
% threshold = xSPM.xSPM.u;
% SPM = load([stats_dir filesep 'SPM.mat']);
% % convolved_task_timecourse = SPM.SPM.xX.X(:,1);
% 
% [clusters, num] = findClustersMNI(tmap_fn, threshold);
% 
% % [Tmax, Tind] = max(clusters{:,2})
% Tind = 0;
% max_val = 0;
% val = 0;
% for c = 1:num
%     if clusters{c,2} > val
%         val = clusters{c,2};
%         Tind = c;
%     end
% end
% 
% I_cluster = clusters{Tind,1}(:,1);
% tmap = spm_read_vols(spm_vol(tmap_fn));
% [Ni,Nj,Nk] = size(tmap);
% cluster_map = zeros(Ni,Nj,Nk);
% cluster_map = cluster_map(:);
% cluster_map(I_cluster) = 1;
% cluster_map_img = reshape(cluster_map, Ni,Nj,Nk);
% % cluster_nifti = [data_dir_new filesep subj_dir filesep 'main_cluster.nii'];
% % % disp('Creating cluster nifti  ...')
% % if ~exist(cluster_nifti, 'file')
% %     new_nii = make_nii(cluster_map_img, [3.5 3.5 4.5], [27 37 17]);
% %     save_nii(new_nii, cluster_nifti)
% % end
% 
% %% display clusters
% fff = displayMaskContour(tmap, cluster_map_img, 0, 4);
% 
% 
% %%
% % print(montage1.f, 'tsnr_whole', '-dpng')
% % montage4 = createMontage(tSNR_3D_masked, 5, 1, 'tSNR (brain)', 'hot');
% 
% 
% 
% 
% % % The plot (with FD, DVARS, and mean Zscore per volume)
% % GM_img = F2D_psc(I_GM, :);
% % WM_img = F2D_psc(I_WM, :);
% % CSF_img = F2D_psc(I_CSF, :);
% % all_img = [GM_img; WM_img; CSF_img];
% % line1_pos = numel(I_GM);
% % line2_pos = numel(I_GM) + numel(I_WM);
% % tf = figure;
% % fontsizeL = 14;
% % fontsizeM = 11;
% % ax1 = subplot(7,1,4:7);
% % imagesc(ax1, all_img); colormap(gray); caxis(intensity_scale);
% % title(ax1, 'thePlotSpm','fontsize',fontsizeL)
% % ylabel(ax1, 'Voxels','fontsize',fontsizeM)
% % xlabel(ax1, 'fMRI volumes','fontsize',fontsizeM)
% % hold on; line([1 Nt],[line1_pos line1_pos],  'Color', 'b', 'LineWidth', 2 )
% % line([1 Nt],[line2_pos line2_pos],  'Color', 'r', 'LineWidth', 2 )
% % hold off;
% % ax2 = subplot(7,1,1);
% % plot(ax2, FD_measures.FD, 'LineWidth', 2); grid;
% % set(ax2,'Xticklabel',[]);
% % title(ax2, 'FD','fontsize',fontsizeL)
% % ylabel(ax2, 'mm','fontsize',fontsizeM)
% % ax3 = subplot(7,1,2);
% % plot(ax3, DVARS, 'LineWidth', 2); grid;
% % set(ax3,'Xticklabel',[]);
% % title(ax3, 'DVARS','fontsize',fontsizeL)
% % ylabel(ax3, 'a.u.','fontsize',fontsizeM)
% % ax4 = subplot(7,1,3);
% % plot(ax4, F2D_zstat_mean, 'LineWidth', 2); grid;
% % set(ax4,'Xticklabel',[]);
% % title(ax4, 'Z-score','fontsize',fontsizeL)
% % ylabel(ax4, 'a.u.','fontsize',fontsizeM)
% % print(tf, 'timeseries_summary', '-dpng')
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %% PIPELINE 3: COMBINED ECHO PROCESSING (STANDARD GLM DENOISING OF COMBINED ECHO IMAGES)
% disp('PIPELINE 3: COMBINED ECHO PROCESSING ...')
% 
% disp('Combining smoothed multi-echo images ...')
% S_combined_img = combine_echoes(TE, T2star_img, I_mask, 1, T1_img, T2_img, T3_img);
% % S_combined_img = combine_echoes([echo time vector], T2star_img, I_overlap, [method (1 = timeseries, 2 = average, 3 = average blocks)], T1_img, T2_img, T3_img);
% 
% %Then create new nifti
% disp('Creating cleaned combined nifti  ...')
% if ~exist(analysis_nii_fn{3}, 'file')
%     new_nii = make_nii(S_combined_img, [3.5 3.5 4.5], [27 37 17]);
%     save_nii(new_nii, analysis_nii_fn{3})
% end
% 
% %Then run stats
% disp('Running SPM stats on cleaned combined nifti ...')
% if ~exist([mni_stats_dir{3} filesep 'spmT_0001.nii'], 'file')
%     runSpmStats(jobs_dir, mni_stats_dir{3}, analysis_nii_fn{3}, Nt, Nregr, regressors_fn{2})
%     pause(5);
%     save([mni_stats_dir{3} filesep 'xSPM.mat'], 'xSPM')
% end
% 
% %% correlation stuff
% disp('Correlation calcs')
% tmap_fn = ['/Users/jheunis/Documents/MATLAB/ME/Volunteer_data/V12all/' subj_dir '/Task/Pipeline 2 - TE2/MNIstats/spmT_0001.nii'];
% tmap = spm_read_vols(spm_vol(tmap_fn));
% load(['/Users/jheunis/Documents/MATLAB/ME/Volunteer_data/V12all/' subj_dir '/Task/Pipeline 2 - TE2/MNIstats/xSPM.mat'])
% tmap_thresholded = (tmap > xSPM.u);
% tmap_thresholded0 = (tmap > 0);
% 
% [vals, inds] = sort(tmap(:), 'descend');
% vals_thresholded = vals > xSPM.u;
% Nu = numel(find(vals_thresholded));
% vals_thresholded0 = vals > 0;
% N0 = numel(find(vals_thresholded0));
% 
% [r0, p0] = corr(T2star(inds(1:N0), :)', S0(inds(1:N0), :)');
% 
% r_diag0 = zeros(1,N0);
% for j = 1:numel(r_diag0)
%     r_diag0(j) = r0(j,j);
% end
% new_vals = flipud(vals(1:N0));
% new_r = fliplr(r_diag0);
% figure; plot(new_vals, new_r, '.r'); grid; hold on;
% line([xSPM.u xSPM.u], [-1 1]);
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% title(['Subject-' num2str(subj) '-' num2str(run) ': Correlation between T2star and S0 (per voxel) for Tmap voxels where T>0 (TE2)']);
% xlabel(['T-values from Tmap (SPM threshold at ' num2str(xSPM.u) ')']);
% ylabel('Pearson coefficient R per voxel');
% 
% 
% 
% disp('Correlation calcs 2')
% tmap_fn = ['/Users/jheunis/Documents/MATLAB/ME/Volunteer_data/V12all/' subj_dir '/Task/Pipeline 3 - Combined_TE/MNIstats/spmT_0001.nii'];
% tmap = spm_read_vols(spm_vol(tmap_fn));
% load(['/Users/jheunis/Documents/MATLAB/ME/Volunteer_data/V12all/' subj_dir '/Task/Pipeline 3 - Combined_TE/MNIstats/xSPM.mat'])
% tmap_thresholded = (tmap > xSPM.u);
% tmap_thresholded0 = (tmap > 0);
% 
% [vals, inds] = sort(tmap(:), 'descend');
% vals_thresholded = vals > xSPM.u;
% Nu = numel(find(vals_thresholded));
% vals_thresholded0 = vals > 0;
% N0 = numel(find(vals_thresholded0));
% 
% [r0, p0] = corr(T2star(inds(1:N0), :)', S0(inds(1:N0), :)');
% 
% r_diag0 = zeros(1,N0);
% for j = 1:numel(r_diag0)
%     r_diag0(j) = r0(j,j);
% end
% new_vals = flipud(vals(1:N0));
% new_r = fliplr(r_diag0);
% figure; plot(new_vals, new_r, '.r'); grid; hold on;
% line([xSPM.u xSPM.u], [-1 1]);
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% title(['Subject-' num2str(subj) '-' num2str(run) ': Correlation between T2star and S0 (per voxel) for Tmap voxels where T>0 (Combined_TE)']);
% xlabel(['T-values from Tmap (SPM threshold at ' num2str(xSPM.u) ')']);
% ylabel('Pearson coefficient R per voxel');
% 

