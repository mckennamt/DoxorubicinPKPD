function [allParms, allParms_ci, allInitGuess, allAICs, numParms] = ...
    fitEachTreatmentModel_regularized(...
    resampd_TrainingData, fixedParms, postTreat_TP, allTspan_training,...
    regConstant, doxPss, resampd_TrainingWeights)

% Initialize fits
kp = fixedParms(1);
theta = fixedParms(2);

switch fixedParms(3)
    case 1
        numParms = 2;
        allParms_names = {'kd','theta'};
        lBound = [-2*kp .7*theta];
        uBound = [10*kp 2*theta];
        normParms = [2*kp 1.1*theta];
    case 2
        numParms = 3;
        allParms_names = {'kd','delay','delay2'};
        lBound = [-Inf 0 0];
        uBound = [Inf 1e6 1e6];
    case 3
        numParms = 3;
        allParms_names = {'kd','delay','delay2'};
        lBound = [-Inf 0 0];
        uBound = [Inf 1e6 1e6];
    case 4
        numParms = 4;
        allParms_names = {'kd','r','r/p','theta'};
        lBound = [-1 0 1 .7*theta];
        uBound = [1 1 Inf 1.1*theta];
        
    case 5
        numParms = 3;
        allParms_names = {'kd', 'r', 'theta'};
        lBound = [-2*kp .001 .7*theta];
        uBound = [5*kp 0.05 2*theta];
        normParms = [2*kp .05 1.2*theta];
    case 6
        numParms = 3;
        allParms_names = {'kd','r', 'theta'};
        lBound = [-1 0 .7*theta];
        uBound = [5 0.5 1.1*theta];
        
    case 7
        numParms = 4;
        allParms_names = {'kd','r','tau', 'theta'};
        lBound = [-1 0 0 .7*theta];
        uBound = [5 0.5 1000 1.1*theta];
        
    case 8
        numParms = 3;
        allParms_names = {'kd','r', 'theta'};
        lBound = [-1 0 .7*theta];
        uBound = [1 0.5 1.1*theta];
end

numTreatConds = length(resampd_TrainingData);

%fixedParms = [kp theta functionSelectIter];
%postTreat_TP = cellData.DrugAdded_Tp + 1;

% place-holder variables for old response model
Yd = [];
time_vector = [];
maxDrugConc = [];

% since simultaneously updating parameters, need to define lower and upper
% bound matrices
lBound = repmat(lBound, numTreatConds, 1);
uBound = repmat(uBound, numTreatConds, 1);

if isempty(doxPss)
    regularizer{1} = 'small';
    % when biasing toward small values
    hld_normalizer(:,1) = lBound(:);
    hld_normalizer(:,2) = uBound(:)-lBound(:);
    regularizer{2} = hld_normalizer;
    sz_Jacob = 1;
else
    regularizer{1} = 'smooth';
    % biasing toward smooth
    regularizer{2} = doxPss;
    regularizer{3} = normParms;
    sz_Jacob = size(doxPss, 2);
end

% calculate size of jacobian to speed optimization
% Jacobian has # of parameter columns
% columns: # treatment conditions * # samples per condition + 2 (2
% predictor variables)
numSamples_perConc = numel(resampd_TrainingData{1}(postTreat_TP:end,:));
jacobPattern = zeros(numSamples_perConc*numTreatConds+2*numel(lBound),...
    numel(lBound));
for pIter = 1:numTreatConds
    sr = (pIter-1)*numSamples_perConc + 1;
    er =  pIter*numSamples_perConc;
    jacobPattern(sr:er, pIter:numTreatConds:end) = 1;
end
% first regularizer term (sorted)

% diagonal (only first and last)
cols = repmat([1 numTreatConds]', numParms, 1);
adder = sort(repmat(numTreatConds*(0:(numParms-1)), 1, 2))';
ioi = sub2ind(size(jacobPattern), er+cols+adder,...
    cols+adder)';

% plus_1
cols = repmat(2:numTreatConds, 1, numParms)';
rows = repmat(1:(numTreatConds-1), 1, numParms)';
adder = sort(repmat(numTreatConds.*(0:(numParms-1)), 1, numTreatConds-1))';
ioi = [ioi sub2ind(size(jacobPattern), er+rows+adder,...
    cols+adder)'];

% minus_1
cols = repmat(1:(numTreatConds-1), 1, numParms)';
rows = repmat(2:(numTreatConds), 1, numParms)';
ioi = [ioi sub2ind(size(jacobPattern), er+rows+adder,...
    cols+adder)'];

ioi = [ioi sub2ind(size(jacobPattern), er+(1:numel(lBound)), 1:size(jacobPattern,2))];
%ioi = [ioi sub2ind(size(jacobPattern), er+(1:numel(lBound)-1), 2:size(jacobPattern,2))];
%ioi = [ioi sub2ind(size(jacobPattern), er+(2:numel(lBound)), 1:size(jacobPattern,2)-1)];

if sz_Jacob == 2
    % second regularizer term
    [~, srtedInds] = sort(doxPss(:,2));
    er = er+numel(lBound);
    
    % diagonal (only first and last)
    cols = repmat([1 numTreatConds]', numParms, 1);
    adder = sort(repmat(numTreatConds*(0:(numParms-1)), 1, 2))';
    ioi = [ioi sub2ind(size(jacobPattern), er+cols+adder,...
        cols+adder)'];
    ioi = [ioi sub2ind(size(jacobPattern), er+(1:numel(lBound)), 1:size(jacobPattern,2))];
    
    % plus_1
    cols = repmat(srtedInds(2:end), numParms, 1);
    rows = repmat([1:(numTreatConds-1)]', numParms, 1);
    adder = sort(repmat(numTreatConds.*(0:(numParms-1)), 1, numTreatConds-1))';
    ioi = [ioi sub2ind(size(jacobPattern), er+rows+adder,...
        cols+adder)'];
    
    % minus_1
    cols = repmat(srtedInds(1:end-1), numParms, 1);
    rows = repmat([2:(numTreatConds)]', numParms, 1);
    ioi = [ioi sub2ind(size(jacobPattern), er+rows+adder,...
        cols+adder)'];
end


jacobPattern(ioi) = 1;
jacobPattern = sparse(jacobPattern);

% set up multi-start problem
lsqOpts = optimoptions('lsqnonlin','MaxFunEvals',10000,'Display','none',...
    'TolFun', 1e-50,'MaxIter', 5000,...
    'TolX', 1e-6, 'JacobPattern',jacobPattern);%, 'Jacobian','on');
%lsqOpts = optimoptions(lsqOpts, 'Algorithm', 'levenberg-marquardt');

calcJacobian = false;

prob = createOptimProblem('lsqnonlin','x0',zeros(size(lBound(:))),...
    'objective',@(x) runDamageModel_regularized(x,resampd_TrainingData,...
    allTspan_training, fixedParms, postTreat_TP, regConstant, regularizer, calcJacobian,...
    resampd_TrainingWeights),...
    'lb',lBound(:),'ub',uBound(:),'options',lsqOpts);

ms = MultiStart('UseParallel',true,'Display','none');
%ms = MultiStart('PlotFcns',@gsplotbestf);

% create parameter estimates
parmEsts = cell(1,numParms);
for iter = 1:numParms
    parmEsts{iter} = linspace(lBound(1,iter),uBound(1,iter),5);
end

hldGrid = cell(1,numel(parmEsts));

[hldGrid{:}] = ndgrid(parmEsts{:});
hldGrid = cellfun(@(x) x(:), hldGrid,'UniformOutput',false);

allEstimates = cat(2,hldGrid{:});

% randomly sample 50 start points
%rng('default'); % set to make sure results are repeatable
numStartPoints = 50;
startPointList = zeros(numStartPoints,numTreatConds*numParms);
for spi = 1:numStartPoints
    stPts = datasample(allEstimates,numTreatConds,1,'Replace',true);
    startPointList(spi,:) = stPts(:)';
end

% vectorize to match parameters
custpts = CustomStartPointSet(startPointList);
[allParms_vector,~,~,~,ms_Solutions] = run(ms,prob,custpts);

allParms = reshape(allParms_vector, numTreatConds, numParms);
bestGuess = ms_Solutions(1).X0{1};
allInitGuess = reshape(bestGuess, numTreatConds, numParms);

% calculate confidence interval for each parameter
lsqOpts = optimoptions(lsqOpts, 'MaxIter',2, 'MaxFunEvals',100);
[optParm,~,res,~,~,~,J] = ...
    lsqnonlin(@(x) runDamageModel_regularized(x,resampd_TrainingData,...
    allTspan_training, fixedParms, postTreat_TP, regConstant, regularizer, calcJacobian,...
    resampd_TrainingWeights),...
    allParms_vector, lBound(:), uBound(:), lsqOpts);
optParm_ci_vector = nlparci(optParm,res,'jacobian',J);

op_ci_rs = zeros(numParms,2,numTreatConds);
for reshapeCI = 1:numTreatConds
    op_ci_rs(:,1,reshapeCI) = optParm_ci_vector(reshapeCI:numTreatConds:end,1);
    op_ci_rs(:,2,reshapeCI) = optParm_ci_vector(reshapeCI:numTreatConds:end,2);
end
allParms_ci = op_ci_rs;

% calculate AIC for each parameter estimate
allAICs = zeros(numTreatConds,1);
k = numParms + 1;
for treatmentIter = 1:numTreatConds
    
    cYall = resampd_TrainingData{treatmentIter};
    cYall = cYall(postTreat_TP:end,:);
    
    Tspan = allTspan_training{treatmentIter};
    Tspan_modified = Tspan(postTreat_TP:end);
    
    if isempty(resampd_TrainingWeights)
        cWeights = length(Tspan_modified) .* ones(1, size(cYall,2));
    else
        cWeights = resampd_TrainingWeights{treatmentIter};
        cWeights = cWeights - (postTreat_TP-1);
    end
    
    optParm = allParms(treatmentIter,:);
    [resNorm,n] = runDamageModel_AICCalc(optParm,cYall,Tspan_modified,...
        fixedParms,maxDrugConc,Yd,time_vector,cWeights);
    
    allAICs(treatmentIter) = n*log(resNorm/n)+(2*k)+(((2*k)*(k+1))/(n-k-1));
    
end


end