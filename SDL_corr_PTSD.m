function SDL_corr_PTSD(SDL)

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
% make fID to merge neuroimaging & demographic/clinical tables
for i = 1:size(T,1) % per subject
    T.fID{i} = [T.SiteName{i},'_',T.SbjID{i}(5:end)];
end
% load covariates, this file is the final file for R statistical analyses
T1 = readtable(fullfile(SDL.path_out,'All_harmonized.csv')); % using previously gerated 
% remove existed V1 to V91
for j = 1:91 % per neurimaging feature
    T1(:,['V',num2str(j)]) = [];
end
% merge neuroimaging data & covs
T0 = innerjoin(T1,T,'Keys','fID'); % merge covars & data by fID

% %% PTSD data only
% T1 = readtable(fullfile(SDL.path_out,'SbjDemoClinScanInfo1.csv')); % load all covariates, this file is R's output to normalize Group, Age & Sex values
% T1 = T1(strcmp(T1.Group,'PTSD') | strcmp(T1.Group,'Control'),:); % remove those without clear group labels
% % make fID, the keys to merge data & covs
% for i = 1:size(T,1) % per subject in the covariates file
%     tmp = extractAfter(T.SbjID{i},'-'); % extract the name after '_'
%     T.fID{i} = [T.SiteName{i},'_',tmp];
% end
% 
% % merge neuroimaging data & covs
% T0 = innerjoin(T1,T,'Keys','fID'); % merge covars & data by fID
% % change variablenames
% T0.Properties.VariableNames{'ScannerSite'} = 'SITE';
% T0.Properties.VariableNames{'x___StudySite'} = 'StudyName';
% T0(find(strcmp(T0.SITE,'CausCon')),:) = []; % remove the scanner site of no interest, i.e. 'Causon' from Stanford

% final cleaned data & cov
a = []; for i=1:size(V,2), a{i} = ['V',num2str(i)]; end % make the list of feature names, i.e. V1 to V91
T = T0(:,a); % only data, i.e. V1 to V91
T1 = T0(:,{'fID','SITE','Group','Age','Sex'});
T0 = [T1,T]; % the final cleaned data & covs
% replace NaN with median value by (SITE,Group)
ListSite = unique(T0.SITE); % list of site names
ListGroup = unique(T0.Group); % list of Group name
for i = 1:size(ListSite) % per site
    for j = 1:size(ListGroup) % per group
        for k = 1:size(V,2) % per connection
            try
                idx = intersect(find(strcmp(T0.SITE,ListSite{i})), find(strcmp(T0.Group,ListGroup{j}))); % index of the subjects of interested SITE and Group
            catch
                idx = find((T0.SITE==ListSite(i))&(T0.Group==ListGroup(j)));
            end
            vec = T0{idx,['V',num2str(k)]};
            if any(isnan(vec)) % if there is at least one NaN
                vec(isnan(vec),1) = nanmedian(vec); % replace NaN with group median value
                T0{idx,['V',num2str(k)]} = vec;
            end
        end
    end
end


% save data for Harmonization (ComBat-GAM)
T = T0(:, a); % only data, i.e. V1 to V91
writetable(T,fullfile(SDL.path_out,'data.csv'));


% %save cov for Harmonization (ComBat-GAM)
% T = T0(:,{'SITE','Group','Age','Sex'});
% T{strcmp(T.SITE,'BOOSTER study'),            'SITE1'} = 6; % change SITE name, 6~10, to concatenate with MDD data (1~5 sites)
% T{strcmp(T.SITE,'BRAINS'),                   'SITE1'} = 7;
% T{strcmp(T.SITE,'Holocaust survivors'),      'SITE1'} = 8;
% T{strcmp(T.SITE,'Larson-Milwaukee'),         'SITE1'} = 9;
% T{strcmp(T.SITE,'Wenchuan Earthquake Study'),'SITE1'} = 10;
% 
% T{strcmp(T.Group,'PTSD'),            'Group1'} = 1; % relabel Group, PTSD=1, Control=0
% T{strcmp(T.Group,'Control'),         'Group1'} = 0; % relabel Group, PTSD=1, Control=0
% 
% T{strcmp(T.Sex,'Male'),              'Sex1'} = 1; % relabel Sex, Male=1, Female=2
% T{strcmp(T.Sex,'Female'),            'Sex1'} = 2; % relabel Sex, Male=1, Female=2
% 
% T.SITE = T.SITE1;
% T.Group = T.Group1;
% T.Sex = T.Sex1;

writetable(T0(:,{'SITE','Group','Age','Sex'}),fullfile(SDL.path_out,'cov.csv'));
% save full data
writetable(T0,fullfile(SDL.path_out,'All.csv'));

% %% ComBat-GAM (Matlab runs Python seems do not work)
% copyfile(fullfile(SDL.path_script, 'SDL_ComBatGAM.py'),...
%     fullfile(SDL.path_out, 'SDL_ComBatGAM.py')); 
% cd(SDL.path_out);
% system('python3 SDL_ComBatGAM.py');
% cd(SDL.path_script);

% %% Combine MDD & PTSD data for ComBat-GAM harmonization
% 
% data_PTSD = readtable(fullfile(SDL.path_project,'Output','FindLab90','timeseries_PTSD','data.csv'));
% cov_PTSD  = readtable(fullfile(SDL.path_project,'Output','FindLab90','timeseries_PTSD','cov.csv'));
% data_MDD  = readtable(fullfile(SDL.path_project,'Output','FindLab90','timeseries_MDD','data.csv'));
% cov_MDD   = readtable(fullfile(SDL.path_project,'Output','FindLab90','timeseries_MDD','cov.csv'));
% 
% cov_MDD{cov_MDD.Group==0,'Group'} = -1; % MDD,  group 0->1
% cov_PTSD{cov_PTSD.Group==1,'Group'} = 2;% PTSD, group 1->2
% 
% data = [data_MDD; data_PTSD];
% cov  = [cov_MDD;  cov_PTSD];
% 
% writetable(data,fullfile(SDL.path_out,'data.csv'));
% writetable(cov, fullfile(SDL.path_out,'cov.csv'));

%% End
end