function SDL_HYDRA(SDL)
% HYDRA analyses
% hydra input
% -i, input data file, col#1st = subjects' ID, col#last = Group labels (-1=controls,1=patients)
% -o, output directory
% -z, input covariates file, col#1st = subjects' ID
% -k, the maximum number of clusters,1~10
% -f, the number of cross-validations, 2~10


% % load harmonized data
% T  = readtable(fullfile(SDL.path_out,'data_harmonized.csv')); % load harmonized data
% T0 = readtable(fullfile(SDL.path_out,'All.csv')); % load covs & raw data
% hasMatch = ~cellfun('isempty', regexp(T0.Properties.VariableNames, 'V', 'once')); % index of V1 to V91
% T0{:, T0.Properties.VariableNames(hasMatch)} = table2array(T); % update data using data_harmonized
% T0{T0.Group==0,'Group'} = -1; % controls = -1, patients = 1, for HYDRA
% % save full data (harmonized)
% writetable(T0,fullfile(SDL.path_out,'All_harmonized.csv'));
% % save data_for_HYDRA
% writetable([T0(:,'fID'),T,T0(:,'Group')],fullfile(SDL.path_out,'data_for_HYDRA.csv'));
% % save cov_for_HYDRA
% T0.Age2 = T0.Age.^2; % add a column of age2
% writetable(T0(:,{'fID','Age','Age2','Sex'}),fullfile(SDL.path_out,'cov_for_HYDRA.csv'));

%% data preparation
% load harmonized data
T = readtable(fullfile(SDL.path_out,'data_harmonized.csv')); % load data harmonized by ComBat-GAM
a = []; for i=1:size(T,2), a{i} = ['V',num2str(i)]; end % make the list of feature names, i.e. V1 to V4005
T.Properties.VariableNames = a; % change neiroimaging features to begin with "V"
% load demographic/clinical data
T1 = readtable(fullfile(SDL.path_out,'All_harmonized.csv')); % using previously generated datasheet
idx = find(cellfun(@(s) ~isempty(strfind('V1', s)), T1.Properties.VariableNames)); % locate the column of 'V1'
T1(:,idx:end) = []; % remove previously stored features, e.g. V1 to V91 (or V4005)
% combine 2 tables
T0 = [T1,T];
% write data for HYDRA
T = [T0(:,'fID'), T0(:,a), T0(:,'Group')]; % fID, Feature1, ..., FeatureN, Group
writetable(T,fullfile(SDL.path_out,'data_for_HYDRA.csv'));
% write cov for HYDRA
T0.Age2 = T0.Age .* T0.Age;
T1 = T0(:,{'fID','Age','Age2','Sex'}); % fID, Age, Age2, Sex
writetable(T1,fullfile(SDL.path_out,'cov_for_HYDRA.csv'));



%% HYDRA function
cd(fullfile(SDL.path_script,'HYDRA-hydra_libsvm')); % get into the folder of HYDRA-hydra_libsvm

tic; fprintf('Running: HYDRA\n');
% hydra('-i','test.csv','-o','.','-k',3,'-f',3);

% run hydra
hydra('-i',fullfile(SDL.path_out,'data_for_HYDRA.csv'),...
    '-o',SDL.path_out,...
    '-z',fullfile(SDL.path_out,'cov_for_HYDRA.csv'),...
    '-k',10,'-f',3); % max 10 subgroups, 3 cross-validations

fprintf('Completed\t');toc;

cd(SDL.path_script); % go back to Scripts folder

% % plot the ARI curve
% load('HYDRA_result.mat');
[val,I] = max(ARI); % find the peak
figure; plot(ARI,'b--o','LineWidth',2,'MarkerSize',5); xlabel('K'); ylabel('ARI'); hold on; % plot ARI
plot(I,val,'ro','markerfacecolor','red','MarkerSize',5); % plot the peak
set(gcf,'color','w'); % white background of the whole figure
box off; % remove borderlines




%% End
end