% Main script for HYDRA analysis of rsFC data
% need to set path for SPM12 toolbox first
% by Delin Sun (Duke University Medical Center, 03/30/2021)

%% Paths & parameters
SDL.path_project = fileparts(pwd); % path of the whole project
SDL.path_script  = fullfile(SDL.path_project,'Scripts'); % path of Scripts
SDL.path_ROI     = fullfile(SDL.path_script,'ROIs','FindLab90'); % path of rois
SDL.SiteName     = {}; % 'AMC','Beijing','Masaryk','Stanford','Milwaukee'
SDL.path_out     = fullfile(SDL.path_project,'Output','timeseries'); mkdir(SDL.path_out); % path of outputs
SDL.path_jobs    = fullfile(SDL.path_project,'jobs'); mkdir(SDL.path_jobs); % path for log files from BIAC cluster

fn = '/mnt/BIAC/duhsnas-pri.dhe.duke.edu/dusom_morey/Data/Lab'; % Lab path through Duke BIAC cluster
if isdir(fn) % if we are using BIAC cluster to do remote calculations
    SDL.path_in = fullfile(fn,'beta6_outputs',SDL.SiteName{1},'derivatives','halfpipe'); % path of inputs
else % if we are testing the code using personla PC from home
    SDL.path_in = 'Z:\Data\Lab\beta6_outputs\AMC\derivatives\halfpipe'; % path of inputs
end

addpath(fullfile(SDL.path_script,'spm12')); % add path of spm12 (please revise it according to your spm12 path)

%% Functions
SDL_timeseries(SDL); % extract & save mean time series per ROI per subject
SDL_corr_MDD(SDL); % ROI-ROI correlaton coefficients, & R-to-Z transformation
SDL_corr_PTSD(SDL); % ROI-ROI correlaton coefficients, & R-to-Z transformation
SDL_HYDRA(SDL); % hydra analyses
