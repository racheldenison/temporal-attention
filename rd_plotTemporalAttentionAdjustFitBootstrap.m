% rd_plotTemporalAttentionAdjustFitBootstrap.m

% standard_model = StandardMixtureModel_SD;
% @(data,g,sd)((1-g).*vonmisespdf(data.errors(:),0,deg2k(sd))+(g).*1/360)

%% group i/o
subjectIDs = {'bl','rd','id','ec','ld','en','sj','ml','ca','jl','ew','jx'};
% subjectIDs = {'ew'};
run = 9;
nSubjects = numel(subjectIDs);

plotDistributions = 0;
saveFigs = 0;

groupFigTitle = [sprintf('%s ',subjectIDs{:}) sprintf('(N=%d), run %d', nSubjects, run)];

modelName = 'mixtureNoBias'; % 'mixtureWithBias','mixtureNoBias','swapNoBias', 'swapWithBias'
bootstraps = 1:100;
nBoots = numel(bootstraps);

%% get data
for iSubject = 1:nSubjects    
    %% indiv i/o
    subjectID = subjectIDs{iSubject};
    subject = sprintf('%s_a1_tc100_soa1000-1250', subjectID);
    fprintf('\n%s\n', subjectID)
    
    expName = 'E3_adjust';
    dataDir = pathToExpt('data');
    figDir = pathToExpt('figures');
    dataDir = sprintf('%s/%s/%s/bootstrap/%s', dataDir, expName, subject(1:2), modelName);
    figDir = sprintf('%s/%s/%s/bootstrap/%s', figDir, expName, subject(1:2), modelName);
    
    %% load data
    for iBoot = 1:nBoots
        fprintf('.')
        
        bootRun = bootstraps(iBoot);
        dataFile = dir(sprintf('%s/%s_run%02d_%s_boot%04d.mat', dataDir, subject, run, modelName, bootRun));
        load(sprintf('%s/%s', dataDir, dataFile.name))
        
        % get data
        for iEL = 1:2
            for iV = 1:3
                % get fit parameters for this condition
                p = fit(iV,iEL).maxPosterior;
                switch modelName
                    case 'mixtureWithBias'
                        mu = p(1);
                        g = p(2);
                        sd = p(3);
                    case 'mixtureNoBias'
                        mu = 0;
                        g = p(1);
                        sd = p(2);
                    case 'swapNoBias'
                        mu = 0;
                        g = p(1);
                        B = p(2);
                        sd = p(3);
                    case 'swapWithBias'
                        mu = p(1);
                        g = p(2);
                        B = p(3);
                        sd = p(4);
                    otherwise
                        error('modelName not recognized')
                end
                
                % store fit parameters
                paramsData.absMu(iV,iEL,iSubject,iBoot) = abs(mu);
                paramsData.mu(iV,iEL,iSubject,iBoot) = mu;
                paramsData.g(iV,iEL,iSubject,iBoot) = g;
                paramsData.sd(iV,iEL,iSubject,iBoot) = sd;
                if exist('B','var')
                    paramsData.B(iV,iEL,iSubject,iBoot) = B;
                end
                
            end
        end
    end
end
fprintf('\n\n')

%% param differences (invalid - valid)
fieldNames = fields(paramsData);
for iField = 1:numel(fieldNames)
    fieldName = fieldNames{iField};
    paramsDiff.(fieldName) = squeeze(paramsData.(fieldName)(2,:,:,:) - ...
        paramsData.(fieldName)(1,:,:,:));
end

%% param summary
fieldNames = fields(paramsData);
for iField = 1:numel(fieldNames)
    fieldName = fieldNames{iField};
    paramsMean.(fieldName) = mean(paramsData.(fieldName),4);
    paramsStd.(fieldName) = std(paramsData.(fieldName),0,4);
    paramsSte.(fieldName) = std(paramsData.(fieldName),0,4)./sqrt(nBoots);
    
    paramsMedian.(fieldName) = median(paramsData.(fieldName),4);
    paramsConfInt.(fieldName) = prctile(paramsData.(fieldName),[2.5 97.5],4);
    
    paramsDiffMedian.(fieldName) = median(paramsDiff.(fieldName),3);
    paramsDiffConfInt.(fieldName) = prctile(paramsDiff.(fieldName),[2.5 97.5],3);
end

%% plot fit parameters
targetNames = {'T1','T2'};
validityNames = {'valid','invalid','neutral'};
validityOrder = [1 3 2];
fieldNames = fields(paramsMean);
f = [];

%% bootstrapped parameter distributions
xlims.g = [-0.1 1];
edges.g = 0:.01:1;

fieldName = 'g';
for iSubject = 1:nSubjects
    figure
    for iEL = 1:2
        for iCV = 1:3
            vals = squeeze(paramsData.(fieldName)(iCV,iEL,iSubject,:));
            subplot(3,2,2*(validityOrder(iCV)-1)+iEL)
            n = histc(vals, edges.(fieldName));
            bar(edges.(fieldName), n);
            xlim(xlims.(fieldName))
            title(validityNames{iCV})
        end
    end
    rd_supertitle(subjectIDs{iSubject});
    rd_raiseAxis(gca);
end
        
%% indiv subjects
ylims = [];
ylims.absMu = [-1 8];
ylims.mu = [-8 8];
ylims.g = [0 0.7];
ylims.sd = [0 70];
ylims.B = [0 0.06];
colors = {'b','g','r'};
for iField = 1:numel(fieldNames)
    fieldName = fieldNames{iField};
%     figNames{end+1} = [fieldName 'Indiv'];
    f(end+1) = figure;
    for iEL = 1:2
        subplot(1,2,iEL)
        hold on
        vals = squeeze(paramsMedian.(fieldName)(:,iEL,:));
        confInt = squeeze(paramsConfInt.(fieldName)(:,iEL,:,:));
        
        for iCV = 1:3
            nudge = -0.4 + 0.2*validityOrder(iCV);
            errorbar((1:nSubjects) + nudge, vals(iCV,:), confInt(iCV,:,1), confInt(iCV,:,2), ...
                '.', 'Color', colors{iCV});
        end
        set(gca,'XTick',1:nSubjects)
        set(gca,'XTickLabel',subjectIDs)
        colormap(flag(3))
        xlim([0 nSubjects+1])
        ylim(ylims.(fieldName))
        if iEL==1
            ylabel(fieldName)
            legend(validityNames)
        end
        title(targetNames{iEL})
    end
    rd_supertitle(groupFigTitle);
    rd_raiseAxis(gca);
end

%% indiv subjects, param tradeoffs
xlims = [-32 32];
ylims = [-.4 .4];
figure
for iEL = 1:2
    subplot(1,2,iEL)
    hold on
    errorbar(paramsDiffMedian.sd(iEL,:), paramsDiffMedian.g(iEL,:), ...
        paramsDiffConfInt.g(iEL,:,1), paramsDiffConfInt.g(iEL,:,2), '.r');
    herrorbar(paramsDiffMedian.sd(iEL,:), paramsDiffMedian.g(iEL,:), ...
        paramsDiffConfInt.sd(iEL,:,1), paramsDiffConfInt.sd(iEL,:,2), '.');
    plot(paramsDiffMedian.sd(iEL,:), paramsDiffMedian.g(iEL,:), '.k')
    xlim(xlims)
    ylim(ylims)
    vline(0,'k')
    plot(xlims,[0 0],'k')
    axis square
    xlabel('standard deviation (invalid-valid)')
    ylabel('guess rate (invalid-valid)')
    title(targetNames{iEL})
end

%% group
ylims.absMu = [-1 4];
ylims.mu = [-4 4];
ylims.g = [0 0.16];
ylims.sd = [0 25];
ylims.B = [0 0.06];
for iField = 1:numel(fieldNames)
    fieldName = fieldNames{iField};
%     figNames{end+1} = [fieldName 'Group'];
    f(end+1) = figure;
    for iEL = 1:2
        subplot(1,2,iEL)
        hold on
        b1 = bar(1:3, squeeze(mean(paramsMedian.(fieldName)(validityOrder,iEL,:),3)),'FaceColor',[.5 .5 .5]);
        p1 = errorbar(1:3, squeeze(mean(paramsMedian.(fieldName)(validityOrder,iEL,:),3)), ...
            squeeze(std(paramsMedian.(fieldName)(validityOrder,iEL,:),0,3)./sqrt(nSubjects)),'k','LineStyle','none');
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
    turnallwhite
    groupFigPrefix = sprintf('gE3_N%d_run%02d_%sMAPBootstrap', nSubjects, run, modelName);
    rd_saveAllFigs(f, figNames, groupFigPrefix, [], '-pdf'); %-depsc2, -dpng
end




