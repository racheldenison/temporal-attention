% rd_plotTemporalAttentionAdjustFit.m

% standard_model = StandardMixtureModel_SD;
% @(data,g,sd)((1-g).*vonmisespdf(data.errors(:),0,deg2k(sd))+(g).*1/360)

%% group i/o
subjectIDs = {'bl','rd','id','ec','ld'};
run = 29;
nSubjects = numel(subjectIDs);

saveFigs = 1;

groupFigTitle = [sprintf('%s ',subjectIDs{:}) sprintf('(N=%d), run %d', nSubjects, run)];

% load the fit results from all subjects
% load data/adjust_workspace_20141225.mat
load(sprintf('data/adjust_workspace_run%02d_20150106.mat', run))

%% get data and plot data and fits
for iSubject = 1:nSubjects
    %% indiv i/o
    subjectID = subjectIDs{iSubject};
    subject = sprintf('%s_a1_tc100_soa1000-1250', subjectID);

    expName = 'E3_adjust';
    % dataDir = 'data';
    % figDir = 'figures';
    dataDir = pathToExpt('data');
    figDir = pathToExpt('figures');
    dataDir = sprintf('%s/%s/%s', dataDir, expName, subject(1:2));
    figDir = sprintf('%s/%s/%s', figDir, expName, subject(1:2));
    
    %% load data
    dataFile = dir(sprintf('%s/%s_run%02d*', dataDir, subject, run));
    load(sprintf('%s/%s', dataDir, dataFile(1).name))
    
    % setup
    df = 4;
    xEdges = -90:df:90;
    xgrid = xEdges(1:end-1) + df/2; % bin centers
    errorIdx = strcmp(expt.trials_headers, 'responseError');
    
    fit = indivResults(iSubject).fit;
    
    % get and plot data and model pdfs
    targetNames = {'T1','T2'};
    validityNames = {'valid','invalid','neutral'};
    for iEL = 1:2
        for iV = 1:3
            % get errors for this condition
            errors = results.totals.all{iV,iEL}(:,errorIdx);
            n = histc(errors, xEdges);
            n(end-1) = n(end-1) + n(end); % last element of n contains the count of values exaclty equal to xEdges(end), so just combine it with the previous bin
            n(end) = [];
            
            % get fit parameters for this condition
            p = fit(iV,iEL).posteriorMean;
            mu = p(1);
            g = p(2);
            sd = p(3);
            
            % store fit parameters
            paramsData.mu(iV,iEL,iSubject) = mu;
            paramsData.g(iV,iEL,iSubject) = g;
            paramsData.sd(iV,iEL,iSubject) = sd;
            
            % generate data and model pdfs (and find residuals) using a common
            % x-axis
            pdfData = (n/sum(n*df))';
            pdfModel = (1-g).*vonmisespdf(xgrid,mu,deg2k(sd))+(g).*1/180;
            resid = pdfData - pdfModel;
            
            % store residuals
            resids(iV,iEL,iSubject,:) = resid;
            
            % also generate smooth model pdf for plotting
            x = -90:90;
            y = (1-g).*vonmisespdf(x,mu,deg2k(sd))+(g).*1/180;
            
            ylims = [-0.02 0.06];
            validityOrder = [1 3 2];
            figure(iSubject)
            subplot(3,2,(validityOrder(iV)-1)*2+iEL)
            hold on
            plot(xgrid,pdfData)
            plot(x,y,'r','LineWidth',1.5)
            %         plot(xgrid,pdfModel,'.r')
            plot(xgrid, resid, 'g')
            ylim(ylims)
            title(sprintf('%s %s', targetNames{iEL}, validityNames{iV}))
        end
    end
    rd_supertitle(sprintf('%s, run %d', subjectID, run));
    if saveFigs
        print(gcf, '-depsc2', ...
            sprintf('%s/%s_run%02d_TemporalAttentionAdjust_fit', figDir, subject, run))
    end
end

%% plot average residuals
validityOrder = [1 3 2];
ylims = [-0.02 0.02];
figNames{1} = 'residsByCond';
f(1) = figure;
for iV = 1:3
    for iEL = 1:2
        subplot(3,2,(validityOrder(iV)-1)*2+iEL)
        hold on
        plot([-100 100], [0 0], '-k');
        plot(xgrid, squeeze(resids(iV,iEL,:,:)), 'g')
        plot(xgrid, mean(squeeze(resids(iV,iEL,:,:))), 'k', 'LineWidth', 2)
        ylim(ylims)
        title(sprintf('%s %s', targetNames{iEL}, validityNames{iV}))
        if iV==3 && iEL==1
            ylabel('residuals (data-model)');
        end
    end
end
rd_supertitle(groupFigTitle);
rd_raiseAxis(gca);

% all conditions on same plot
residsMean = squeeze(mean(resids,3));
figNames{2} = 'residsAllConds';
f(2) = figure;
hold on
plot(xgrid,squeeze(residsMean(:,1,:))')
plot(xgrid,squeeze(residsMean(:,2,:))')
plot(xgrid,squeeze(mean(mean(residsMean,1),2))','k','LineWidth',2)
plot([-100 100], [0 0], '-k');
legend(validityNames)
ylabel('p(error) residual mean')
rd_supertitle(groupFigTitle);
rd_raiseAxis(gca);


%% param summary
fieldNames = fields(paramsData);
for iField = 1:numel(fieldNames)
    fieldName = fieldNames{iField};
    paramsMean.(fieldName) = mean(paramsData.(fieldName),3);
    paramsSte.(fieldName) = std(paramsData.(fieldName),0,3)./sqrt(nSubjects);
end

%% plot fit parameters
validityOrder = [1 3 2];
fieldNames = fields(paramsMean);

% indiv subjects
ylims = [];
ylims.mu = [-8 8];
ylims.g = [0 0.3];
ylims.sd = [0 30];
for iField = 1:numel(fieldNames)
    fieldName = fieldNames{iField};
    figNames{2+iField} = [fieldName 'Indiv'];
    f(2+iField) = figure;
    for iEL = 1:2
        subplot(1,2,iEL)
        bar(squeeze(paramsData.(fieldName)(validityOrder,iEL,:))')
        set(gca,'XTickLabel',subjectIDs)
        colormap(flag(3))
        ylim(ylims.(fieldName))
        if iEL==1
            ylabel(fieldName)
            legend(validityNames(validityOrder))
        end
        title(targetNames{iEL})
    end
    rd_supertitle(groupFigTitle);
    rd_raiseAxis(gca);
end

% group
ylims.mu = [-4 4];
ylims.g = [0 0.16];
ylims.sd = [0 25];
for iField = 1:numel(fieldNames)
    fieldName = fieldNames{iField};
    figNames{5+iField} = [fieldName 'Group'];
    f(5+iField) = figure;
    for iEL = 1:2
        subplot(1,2,iEL)
        hold on
        b1 = bar(1:3, paramsMean.(fieldName)(validityOrder,iEL),'FaceColor',[.5 .5 .5]);
        p1 = errorbar(1:3, paramsMean.(fieldName)(validityOrder,iEL)', ...
            paramsSte.(fieldName)(validityOrder,iEL)','k','LineStyle','none');
        ylim(ylims.(fieldName))
        ylabel(fieldName)
        set(gca,'XTick',1:3)
        set(gca,'XTickLabel', validityNames(validityOrder))
        title(targetNames{iEL})
    end
    rd_supertitle(groupFigTitle);
    rd_raiseAxis(gca);
end

%% save figures
if saveFigs
    groupFigPrefix = sprintf('gE3_N%d_run%02d', nSubjects, run);
    rd_saveAllFigs(f, figNames, groupFigPrefix, [], '-depsc2');
end
    



