% rd_plotTemporalAttentionAdjustFit.m

% standard_model = StandardMixtureModel_SD;
% @(data,g,sd)((1-g).*vonmisespdf(data.errors(:),0,deg2k(sd))+(g).*1/360)

%% group i/o
subjectIDs = {'bl','rd','id','ec','ld','en','sj','ml','ca','jl','ew','jx'};
% subjectIDs = {'bl'};
run = 9;
nSubjects = numel(subjectIDs);

plotDistributions = 1;
saveFigs = 0;

groupFigTitle = [sprintf('%s ',subjectIDs{:}) sprintf('(N=%d), run %d', nSubjects, run)];

modelName = 'mixtureKurtosis'; % 'fixedNoBias','mixtureWithBias','mixtureNoBias','swapNoBias', 'swapWithBias', 'variablePrecision', 'variablePrecisionGammaSD', 'variablePrecisionNoGuess', 'mixtureKurtosis'

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
    dataFile = dir(sprintf('%s/%s_run%02d_%s.mat', dataDir, subject, run, modelName));
    load(sprintf('%s/%s', dataDir, dataFile.name))
    
    % setup
    df = 4;
    xEdges = -90:df:90;
    xgrid = xEdges(1:end-1) + df/2; % bin centers
    paramNames = model.paramNames;
    nParams = numel(paramNames);
    
    % get and plot data and model pdfs
    targetNames = {'T1','T2'};
    validityNames = {'valid','invalid','neutral'};
    for iEL = 1:2
        for iV = 1:3
            % get errors for this condition
            errors = err{iV,iEL};
            n = histc(errors, xEdges);
            n(end-1) = n(end-1) + n(end); % last element of n contains the count of values exaclty equal to xEdges(end), so just combine it with the previous bin
            n(end) = [];
            
            % get fit parameters for this condition
            p = fit(iV,iEL).maxPosterior; % posteriorMean
            
            for iP = 1:nParams
                pName = paramNames{iP};
                paramsData.(pName)(iV,iEL,iSubject) = p(iP);
                paramsLowerCred.(pName)(iV,iEL,iSubject) = fit(iV,iEL).lowerCredible(iP);
                paramsUpperCred.(pName)(iV,iEL,iSubject) = fit(iV,iEL).upperCredible(iP);
            end
            
            mu = p(strcmp(paramNames,'mu'));
            g = p(strcmp(paramNames,'g'));
            sd = p(strcmp(paramNames,'sd'));
            if isempty(mu)
                mu = 0;
            else
                paramsData.absMu(iV,iEL,iSubject) = abs(mu);
            end
            if isempty(g)
                g = 0;
            end

%             switch modelName
%                 case 'mixtureWithBias'
%                     mu = p(1);
%                     g = p(2);
%                     sd = p(3);
%                 case 'mixtureNoBias'
%                     mu = 0;
%                     g = p(1);
%                     sd = p(2);
%                 case 'swapNoBias'
%                     mu = 0;
%                     g = p(1);
%                     B = p(2);
%                     sd = p(3);
%                 case 'swapWithBias'
%                     mu = p(1);
%                     g = p(2);
%                     B = p(3);
%                     sd = p(4);
%                 case 'variablePrecision'
%                     g = p(1);
%                     mnSTD = p(2);
%                     sdSTD = p(3);
%                 case 'variablePrecisionGammaSD'
%                     g = p(1);
%                     modeSTD = p(2);
%                     sdSTD = p(3);
%                 otherwise
%                     error('modelName not recognized')
%             end
%             
%             % store fit parameters
%             switch modelName
%                 case 'variablePrecision'
%                     paramsData.g(iV,iEL,iSubject) = g;
%                     paramsData.mnSTD(iV,iEL,iSubject) = mnSTD;
%                     paramsData.sdSTD(iV,iEL,iSubject) = sdSTD;
%                 case 'variablePrecisionGammaSD'
%                     paramsData.g(iV,iEL,iSubject) = g;
%                     paramsData.modeSTD(iV,iEL,iSubject) = modeSTD;
%                     paramsData.sdSTD(iV,iEL,iSubject) = sdSTD;
%                 otherwise
%                     paramsData.absMu(iV,iEL,iSubject) = abs(mu);
%                     paramsData.mu(iV,iEL,iSubject) = mu;
%                     paramsData.g(iV,iEL,iSubject) = g;
%                     paramsData.sd(iV,iEL,iSubject) = sd;
%                     if exist('B','var')
%                         paramsData.B(iV,iEL,iSubject) = B;
%                     end
%             end
            
            % generate data and model pdfs (and find residuals) using a common
            % x-axis
            if strfind(modelName,'variablePrecision') | strfind(modelName,'Kurtosis')
                paramsAsCell = num2cell(p);
                model = GetModelPdfForPlot(model);
                p1 = model.pdfForPlot(xgrid, [], paramsAsCell{:});
                pdfModel = (p1/sum(p1*df))';
            else
                pdfModel = (1-g).*vonmisespdf(xgrid,mu,deg2k(sd))+(g).*1/180;
            end
            pdfData = (n/sum(n*df))';
            resid = pdfData - pdfModel;
            
            % store residuals
            resids(iV,iEL,iSubject,:) = resid;
            if exist('mu','var')
                residsShift(iV,iEL,iSubject,:) = circshift(resid,[0 round(mu/df)]);
            else
                residsShift = resids;
            end
            
            % also generate smooth model pdf for plotting
            x = -90:90;
            if strfind(modelName,'variablePrecision') | strfind(modelName,'Kurtosis')
                y1 = model.pdfForPlot(x, [], paramsAsCell{:});
                y = y1/sum(y1);
            else
                y = (1-g).*vonmisespdf(x,mu,deg2k(sd))+(g).*1/180;
            end
            
            if plotDistributions
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
    end
    if plotDistributions
        rd_supertitle(sprintf('%s, run %d', subjectID, run));
        if saveFigs
            print(gcf, '-depsc2', ...
                sprintf('%s/%s_run%02d_TemporalAttentionAdjust_fit_%s', figDir, subject, run, modelName))
        end
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
colors = {[1 0 0],[0 0 1],[0 0 0]};
for iCol = 1:numel(colors)
    barColors{iCol} = desaturate(colors{iCol}, 0.5);
end
barColors{3} = [.5 .5 .5];
fieldNames = fields(paramsMean);

% indiv subjects
ylims = [];
ylims.absMu = [-1 8];
ylims.mu = [-8 8];
ylims.g = [0 0.3];
ylims.sd = [0 30];
ylims.B = [0 0.06];
ylims.modeSTD = [0 30];
ylims.mnSTD = [0 30];
ylims.sdSTD = [0 20];
ylims.stdSTD = [0 20];
for iField = 1:numel(fieldNames)
    fieldName = fieldNames{iField};
    figNames{end+1} = [fieldName 'Indiv'];
    f(end+1) = figure;
    for iEL = 1:2
        subplot(1,2,iEL)
        bar(squeeze(paramsData.(fieldName)(validityOrder,iEL,:))')
        set(gca,'XTickLabel',subjectIDs)
        colormap(flag(3))
        xlim([0 nSubjects+1])
%         ylim(ylims.(fieldName))
        if iEL==1
            ylabel(fieldName)
            legend(validityNames(validityOrder))
        end
        title(targetNames{iEL})
    end
    rd_supertitle(groupFigTitle);
    rd_raiseAxis(gca);
end

% with indiv error bars
ylimsE = [];
ylimsE.g = [-0.1 0.7]; % E = with individual error bars
ylimsE.sd = [0 70];
ylimsE.modeSTD = [0 80];
ylimsE.mnSTD = [0 80];
ylimsE.sdSTD = [-10 60];
ylimsE.stdSTD = [-10 60];
offsets = [-.2 0 .2];
barWidth = 0.15;
for iField = 1:numel(fieldNames)
    fieldName = fieldNames{iField};
    figNames{end+1} = [fieldName 'IndivError'];
    f(end+1) = figure('Position',[100 100 700 400]);
    for iEL = 1:2
        subplot(1,2,iEL)
        hold on
%         bar(squeeze(paramsData.(fieldName)(validityOrder,iEL,:))')
        for iV = 1:3
            b1(iV) = bar((1:nSubjects)+offsets(iV), squeeze(paramsData.(fieldName)(validityOrder(iV),iEL,:)), ...
                barWidth, 'FaceColor', barColors{validityOrder(iV)}, 'EdgeColor', 'none');
            errorbar((1:nSubjects)+offsets(iV), squeeze(paramsData.(fieldName)(validityOrder(iV),iEL,:)), ...
                squeeze(paramsLowerCred.(fieldName)(validityOrder(iV),iEL,:)), ...
                squeeze(paramsUpperCred.(fieldName)(validityOrder(iV),iEL,:)), ...
                'LineStyle','none','color',colors{validityOrder(iV)});
        end
        set(gca,'XTick',1:nSubjects)
        set(gca,'XTickLabel',subjectIDs)
        colormap(flag(3))
        xlim([0 nSubjects+1])
%         ylim(ylimsE.(fieldName))
        if iEL==1
            ylabel(fieldName)
            legend(b1, validityNames(validityOrder))
            legend('boxoff')
        end
        title(targetNames{iEL})
    end
    rd_supertitle(groupFigTitle);
    rd_raiseAxis(gca);
end

% group
ylims.absMu = [-1 4];
ylims.mu = [-4 4];
ylims.g = [0 0.16];
ylims.sd = [0 25];
ylims.B = [0 0.06];
ylims.modeSTD = [0 18];
ylims.mnSTD = [0 22];
ylims.sdSTD = [0 6];
ylims.stdSTD = [0 14];
for iField = 1:numel(fieldNames)
    fieldName = fieldNames{iField};
    figNames{end+1} = [fieldName 'Group'];
    f(end+1) = figure;
    for iEL = 1:2
        subplot(1,2,iEL)
        hold on
        b1 = bar(1:3, paramsMean.(fieldName)(validityOrder,iEL),'FaceColor',[.5 .5 .5]);
        p1 = errorbar(1:3, paramsMean.(fieldName)(validityOrder,iEL)', ...
            paramsSte.(fieldName)(validityOrder,iEL)','k','LineStyle','none');
%         ylim(ylims.(fieldName))
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
    groupFigPrefix = sprintf('gE3_N%d_run%02d_%sMAP', nSubjects, run, modelName);
    rd_saveAllFigs(f, figNames, groupFigPrefix, [], '-pdf'); %-depsc2, -dpng
end
    



