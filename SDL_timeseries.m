function SDL_timeseries(SDL)

% Extract mean time series per ROI from preprocessed rs-fMRI images
% 14 ROIs (networks) are from the findlab (https://findlab.stanford.edu/functional_ROIs.html)
% The input images should be *Preproc1.nii.gz that has motion parameters from aCompCor regressed out
% need to set path for SPM12 first
% by Delin Sun (Duke University Medical Center, 03/30/2021)
%
% Input
% --- SDL, a structure containing the important paths & parameters
% Output
% --- save timeseries into a .mat file per SITE


%% ROIs

ROIs = dir(fullfile(SDL.path_ROI,'*.nii')); % list all 90 ROIs from FindLab, files were SPM-coregistration-reslice to the resolution of functional images

% loading all ROIs
for i = 1:size(ROIs,1)
    tic; froi(i).name = ROIs(i).name;
    froi(i).fullname = fullfile(ROIs(i).folder,   ROIs(i).name);
    froi(i).img3D = spm_read_vols(spm_vol(froi(i).fullname));% load mask image (3D)
    fprintf('Loading ROIs: ROI=%s\t',froi(i).name); toc;
end


%% Functions
for j = 1:size(SDL.SiteName,2) % per StudySite
    % path
    SDL.path_in = fullfile('/mnt/BIAC/duhsnas-pri.dhe.duke.edu/dusom_morey/Data/Lab',...
        'beta6_outputs',SDL.SiteName{j},'derivatives','halfpipe');
    
    % list all subjects
    fins = dir(SDL.path_in);
    
    M = []; % structure to contain the timeseries per ROI per subject
    for i = 3:size(fins,1) % extract signals per subject, i=1,2 are . and ..
        tic; M(i-2).name = fins(i).name; % subject ID
        fprintf('\nExtracting: SbjID=%s\n',fins(i).name);
        
        fstr = dir(fullfile(fins(i).folder,fins(i).name,'func','*preproc1_bold.nii.gz')); % searching the file with specified extension
        fin = fullfile(fullfile(fstr.folder,fstr.name)); % fullname of the inout image
        % e.g. 'Z:\Data\Lab\beta6_outputs\AMC\derivatives\halfpipe\sub-1138\func\sub-1138_task-rest_setting-preproc1_bold
        
        M(i-2).data = SDL_sub_timeseries(fin,froi); % extract mean timeseries per ROI per subject
        toc;
        
    end
    
    % save outputs
    fout = fullfile(SDL.path_out,[SDL.SiteName{j},'.mat']); % data from the same site into the same file
    tic; fprintf('Saving: %s\n',fout);
    save(fout,'M','-v7.3'); toc;
    
end
% End
end


% extract mean timeseries per ROI from a particular subject
function D = SDL_sub_timeseries(fin,froi)
% Input
% --- fin,  full filename of the 4D functional image
% --- froi, a structure containing all ROIs
% Output
% --- D, Matrix of the averaged timeseries per ROI

% tic;
fprintf('Loading Preprocessed functional image:\n %s\n',fin);
f4D = spm_read_vols(spm_vol(fin)); % load preprocessed resting-state functional image (4D)
[~,~,~,Nt] = size(f4D); % size of the time dimension

Nroi = size(froi,2); % number of rois
D = zeros(Nt,Nroi); % matrix to contain outputs
for i = 1:Nroi % per roi
    fm4D = repmat(froi(i).img3D,1,1,1,Nt); % repeat 3D mask in each frame to make a 4D mask
    f = f4D .* fm4D; % mask the functional image
    M = mean(f,[1,2,3],'omitnan'); % average acros x, y, and z dimensions, removing NaN values
    D(:,i) = reshape(M,Nt,1); % save data
end

end