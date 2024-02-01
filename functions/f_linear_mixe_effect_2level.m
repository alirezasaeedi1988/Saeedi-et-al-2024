
function [stimulusPValue,anovaResults]= f_linear_mixe_effect_2level(A1,A2,MID,NANflag)
% %  This code compare neuronal reponse properties to 2 different stimuli and uses the mouse ID as random effect
% %   to account for the nested data structure
% %  Response variable: Neuronal responses e.g. firing rate or latency ....
% %  Fixed effect: Type of stimulus.
% %  Random effect: Mouse ID (to account for inter-mouse variability).
% %  inputs
% %        A1 = 1XN vector of response data to stimulus 1 
% %        A2 = 1XN vector of response data to stimulus 2
% %        MID = 1XN vector of actual Mouse IDs
% %         NANflag = "pair" to omit rows (neurons) if only one reponse is nan, "single", this option
% only delete the nan valuse 

% Pairwise deletion: Remove rows where either A1 or A2 has NaN
% remove nan and its pair
if strcmp(NANflag,"pair")
   
    validIndices = ~isnan(A1) & ~isnan(A2);
    A1 = A1(validIndices);
    A2 = A2(validIndices);
    MID = MID(validIndices,:);
end
MID1 = MID; MID2 = MID1; 
% Organize data into a table
data1 = table(repmat((1:length(A1))',1,1), A1, MID1(:,1),MID1(:,2), ones(length(A1), 1), ...
    'VariableNames', {'Neuron', 'Response', 'MouseID','Day', 'Stimulus'});

data2 = table(repmat((1:length(A2))',1,1), A2, MID2(:,1),MID2(:,2), 2*ones(length(A2), 1), ...
    'VariableNames', {'Neuron', 'Response', 'MouseID','Day', 'Stimulus'});

data = [data1; data2];
% Remove rows with NaN values in the 'Response' column 
% this will do the "single" option in the nanflag anyway!
data = data(~isnan(data.Response), :);

% Fit the linear mixed-effects model
mdl = fitlme(data, 'Response ~ 1 + Stimulus + (1|MouseID)+ (1|MouseID:Day)');
% mdl = fitglme(data, 'Response ~ 1 + Stimulus + (1|MouseID)', 'Distribution', 'Poisson');

% % Display the results
% disp(mdl);

% Assuming mdl is  fitted LinearMixedModel
anovaResults = anova(mdl);

% stimulusEstimate  = mdl.Coefficients.Estimate(2); 
stimulusPValue     = anovaResults.pValue(2);


% Here are some of the key results and properties can be extracted:

% 1. **Coefficients**: Estimated coefficients of the fixed effects, including their standard errors, t-values, and p-values.

%    disp(mdl.Coefficients)

% 2. **Formula**: The formula used in the model.

%    mdl.Formula

% 3. **Random Effects**: Estimates of the random effects.

%    randEff = randomEffects(mdl);

% 4. **Fixed, Random, and Conditional R-squared values**: 

%    mdl.Rsquared

% 5. **Model Information**: Includes AIC (Akaike Information Criterion), BIC (Bayesian Information Criterion), logLik (log-likelihood), and more.

%    mdl.ModelCriterion

% 6. **Diagnostics**: You can extract residuals, leverage, and other diagnostics to check the assumptions of the model.

%    residuals = residuals(mdl);

% 7. **ANOVA Table**: A table of hypothesis tests for the fixed effects and their interactions.

%    anovaTbl = anova(mdl);


% 8. **Variance Components**: Information about the variance attributed to each random effect and the residuals.

%    disp(mdl.VarComp)

% 9. **Model Fitting Statistics**: Details of the optimization algorithm used, number of iterations, convergence status, etc.
%    mdl.FitMethod
%    mdl.LogLikelihood
%    mdl.ConvergenceInfo

% 10. **Covariance Parameters**: The estimated covariance parameters for the random effects.

%    covParms = covarianceParameters(mdl);

% 
% You can use various plotting functions like `plotResiduals`, `plotDiagnostics`, etc., to visually explore different aspects of the model fit.
