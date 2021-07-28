function SDL_corr_MDD(SDL)

% (1) ROI-ROI correlation coefficients per subject
% (2) correlation coefficients R-to-Z transofrmation
% (3) save R-to-Z transformed correlation coefficients into a .mat file
% Input
% --- SDL, a structure containing the important paths & parameters
% Output
% --- save correlation coefficients across subjects and sites into a .mat file

%% Parameters
fn = dir(fullfile(SDL.path_out,'*.mat'));

D = []; % structure to contain sitename, sbjID, & corr. coef.
k = 0; % index of subjects in D
for i = 1:size(fn,1) % per SITE
    clear M;
    load(fullfile(fn(i).folder, fn(i).name),'M'); % load the M structure that contains time series
    
    for j = 1:size(M,2) % per subject
        tic; fprintf('Data: SiteName=%s, SbjID=%s\t',fn(i).name(1:end-4),M(j).name);
        T = M(j).data; % row = time point, col = ROI
        R = corr(T); % correlation coefficients per ROI-ROI pair
        Z = atanh(R); % R-to-Z Fisher transformation
        
        % vectorize the upper triangle of matrix Z (with diagonal)
        k = k + 1; % index of subject in D
        V(k,:)        = Z(triu(true(size(Z)),1)); % vector of corr. coef.
        D(k).SiteName = fn(i).name(1:end-4); % SiteName
        D(k).SbjID    = M(j).name; % SbjID
        toc;
    end
end

%% Organize data & cov file
T = [struct2table(D),array2table(V)]; % concatenate StudyName, SubID & data of neuroimaging measures
T.fID = T.SbjID; % fID column is used to merge two tables


% MDD data only
T1 = readtable(fullfile(SDL.path_project,'Origins','MDD','ALL_covariates.csv')); % load all covariates
% make fID, the keys to merge data & covs
for i = 1:size(T1,1) % per subject in the covariates file
    tmp = extractAfter(T1.SubjID{i},'_'); % extract the name after '_'
    if ~isempty(tmp) % if the raw name looks like 'ABC_123'
        T1.fID{i} = ['sub-',tmp];
    else % if the raw name looks like '123'
        T1.fID{i} = ['sub-',T1.SubjID{i}];
    end
end

% merge neuroimaging data & covs
T0 = innerjoin(T1,T,'Keys','fID'); % merge covars & data by fID
% change variablenames
T0.Properties.VariableNames{'Dx'} = 'Group';
T0.Properties.VariableNames{'Site'} = 'SITE';
T0.Properties.VariableNames{'SiteName'} = 'StudyName';

% final cleaned data & cov
hasMatch = ~cellfun('isempty', regexp(T0.Properties.VariableNames, 'V', 'once')); % index of V1 to V91 (or V4005)
T = T0(:, T0.Properties.VariableNames(hasMatch)); % only data, i.e. V1 to V91 (or V4005)
T1 = T0(:,{'fID','SITE','Group','Age','Sex','Recur','AD','Rem','AO','Sev','PTSD_Current','StudyName'});
T0 = [T1,T]; % the final cleaned data & covs
% replace NaN with median value by (SITE,Group)
ListSite = unique(T0.SITE); % list of site names
ListGroup = unique(T0.Group); % list of Group name
T0(T0.Sex==3,:) = []; % Sex can only be 1 (male) and female (2)
for i = 1:size(ListSite) % per site
    for j = 1:size(ListGroup) % per group
        for k = 1:size(V,2) % per connection
            vec = T0{(T0.SITE==ListSite(i))&(T0.Group==ListGroup(j)),['V',num2str(k)]};
            if any(isnan(vec)) % if there is at least one NaN
                vec(isnan(vec),1) = nanmedian(vec); % replace NaN with group median value
                T0{(T0.SITE==ListSite(i))&(T0.Group==ListGroup(j)),['V',num2str(k)]} = vec;
            end
        end
    end
end

sum(sum(isnan(T0{:,13:4017}))); % 0 if there is no NaN in V1 to V4005

% show the number of subjects per site
B = unique(T0.SITE); % unique site names
disp('Number of subjects per SITE:');
[B,histc(T0.SITE,B)]


% save data for Haronization (ComBat-GAM)
hasMatch = ~cellfun('isempty', regexp(T0.Properties.VariableNames, 'V', 'once')); % index of V1 to V91
T = T0(:, T0.Properties.VariableNames(hasMatch)); % only data, i.e. V1 to V91
writetable(T,fullfile(SDL.path_out,'data.csv'));
%save cov for Harmonization (ComBat-GAM)
T = T0(:,{'SITE','Group','Age','Sex'});
writetable(T,fullfile(SDL.path_out,'cov.csv'));
% save full data
writetable(T0,fullfile(SDL.path_out,'All.csv'));



%% End
end