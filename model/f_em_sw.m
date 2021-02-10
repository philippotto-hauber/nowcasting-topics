function [Fhat, x2] = f_em_sw(x,k,DEMEAN)
% =========================================================================
% =========================================================================
% This function is an adjusted version of the EM-algorithm by Stock and Watson (2002)
% as implemented by McCracken and Ng (2016). Compared to their original code
% this implementation fixes the number of factors k instead of choosing it
% based on an information criterion and hardcodes the transformation so that
% all series are demeaned and standardized. In terms of output only the factors
% and complete data set with missings filled in are returned. Note that the
% output matrices are transposed so that Fhat is an k x T matrix and 
% x2 is an N x T matrix. The rest is unaltered. PH 10.2.2021
% =========================================================================
% Orginal code begins here
% =========================================================================
% DESCRIPTION
% This program estimates a set of factors for a given dataset using
% principal component analysis. The number of factors estimated is
% determined by an information criterion specified by the user. Missing
% values in the original dataset are handled using an iterative
% expectation-maximization (EM) algorithm.
%
% -------------------------------------------------------------------------
% INPUTS
%           x       = dataset (one series per column)
%           k       = number of factors  %
% OUTPUTS
%           Fhat    = set of factors
%           x2      = x with missing values replaced from the EM algorithm
%
% -------------------------------------------------------------------------
% SUBFUNCTIONS
%
% baing() - selects number of factors
%
% pc2() - runs principal component analysis
%
% minindc() - finds the index of the minimum value for each column of a
%       given matrix
%
% transform_data() - performs data transformation
%
% -------------------------------------------------------------------------
% BREAKDOWN OF THE FUNCTION
%
% Part 1: Check that inputs are specified correctly.
%
% Part 2: Setup.
%
% Part 3: Initialize the EM algorithm -- fill in missing values with
%         unconditional mean and estimate factors using the updated
%         dataset.
%
% Part 4: Perform the EM algorithm -- update missing values using factors,
%         construct a new set of factors from the updated dataset, and
%         repeat until the factor estimates do not change.
% 
% -------------------------------------------------------------------------
% NOTES
% Authors: Michael W. McCracken and Serena Ng
% Date: 9/5/2017
% Version: MATLAB 2014a
% Required Toolboxes: None
%
% Details for the three possible information criteria can be found in the
% paper "Determining the Number of Factors in Approximate Factor Models" by
% Bai and Ng (2002).
%
% The EM algorithm is essentially the one given in the paper "Macroeconomic
% Forecasting Using Diffusion Indexes" by Stock and Watson (2002). The
% algorithm is initialized by filling in missing values with the
% unconditional mean of the series, demeaning and standardizing the updated
% dataset, estimating factors from this demeaned and standardized dataset,
% and then using these factors to predict the dataset. The algorithm then
% proceeds as follows: update missing values using values predicted by the
% latest set of factors, demean and standardize the updated dataset,
% estimate a new set of factors using the demeaned and standardized updated
% dataset, and repeat the process until the factor estimates do not change.
%
% =========================================================================
% PART 1: CHECKS

% Check that x is not missing values for an entire row
if sum(sum(isnan(x),2)==size(x,2))>0
    error('Input x contains entire row of missing values.');
end

% Check that x is not missing values for an entire column
if sum(sum(isnan(x),1)==size(x,1))>0
    error('Input x contains entire column of missing values.');
end

% =========================================================================
% PART 2: SETUP

% Maximum number of iterations for the EM algorithm
maxit=50;

% Number of observations per series in x (i.e. number of rows)
T=size(x,1);

% Number of series in x (i.e. number of columns)
N=size(x,2);

% Set error to arbitrarily high number
err=999;

% Set iteration counter to 0
it=0;

% Locate missing values in x
x1=isnan(x);

% =========================================================================
% PART 3: INITIALIZE EM ALGORITHM
% Fill in missing values for each series with the unconditional mean of
% that series. Demean and standardize the updated dataset. Estimate factors
% using the demeaned and standardized dataset, and use these factors to
% predict the original dataset.

% Get unconditional mean of the non-missing values of each series
mut=repmat(nanmean(x),T,1);

% Replace missing values with unconditional mean
x2=x;
x2(isnan(x))=mut(isnan(x));

% Demean and standardize data using subfunction transform_data()
%   x3  = transformed dataset
%   mut = matrix containing the values subtracted from x2 during the
%         transformation
%   sdt = matrix containing the values that x2 was divided by during the
%         transformation
[x3,mut,sdt]=transform_data(x2,2);

% Set number of factors equal to k 
icstar=k;

% Run principal components on updated dataset using subfunction pc2()
%   chat   = values of x3 predicted by the factors
%   Fhat   = factors scaled by (1/sqrt(N)) where N is the number of series
%   lamhat = factor loadings scaled by number of series
%   ve2    = eigenvalues of x3'*x3 
[chat,Fhat,lamhat,ve2]  = pc2(x3,icstar);

% Save predicted series values
chat0=chat;

% =========================================================================
% PART 4: PERFORM EM ALGORITHM
% Update missing values using values predicted by the latest set of
% factors. Demean and standardize the updated dataset. Estimate a new set
% of factors using the updated dataset. Repeat the process until the factor
% estimates do not change.

% Run while error is large and have yet to exceed maximum number of
% iterations
while err> 0.000001 && it <maxit
    
    % ---------------------------------------------------------------------
    % INCREASE ITERATION COUNTER
    
    % Increase iteration counter by 1
    it=it+1;
    
    % Display iteration counter, error, and number of factors
    fprintf('Iteration %d: obj %10f IC %d \n',it,err,icstar);

    % ---------------------------------------------------------------------
    % UPDATE MISSING VALUES
    
    % Replace missing observations with latest values predicted by the
    % factors (after undoing any transformation)
    for t=1:T;
        for j=1:N;
            if x1(t,j)==1
                x2(t,j)=chat(t,j)*sdt(t,j)+mut(t,j);    
            else
                x2(t,j)=x(t,j);
            end
        end
    end
    
    % ---------------------------------------------------------------------
    % ESTIMATE FACTORS
    
    % Demean/standardize new dataset and recalculate mut and sdt using
    % subfunction transform_data()
    %   x3  = transformed dataset
    %   mut = matrix containing the values subtracted from x2 during the
    %         transformation
    %   sdt = matrix containing the values that x2 was divided by during 
    %         the transformation
    [x3,mut,sdt]=transform_data(x2,2);
        
    % Set number of factors equal to k 
    icstar=k;

    % Run principal components on the new dataset using subfunction pc2()
    %   chat   = values of x3 predicted by the factors
    %   Fhat   = factors scaled by (1/sqrt(N)) where N is the number of 
    %            series
    %   lamhat = factor loadings scaled by number of series
    %   ve2    = eigenvalues of x3'*x3 
    [chat,Fhat,lamhat,ve2]  = pc2(x3,icstar);

    % ---------------------------------------------------------------------
    % CALCULATE NEW ERROR VALUE
    
    % Caclulate difference between the predicted values of the new dataset
    % and the predicted values of the previous dataset
    diff=chat-chat0;
    
    % The error value is equal to the sum of the squared differences
    % between chat and chat0 divided by the sum of the squared values of
    % chat0
    v1=diff(:);
    v2=chat0(:);
    err=(v1'*v1)/(v2'*v2);

    % Set chat0 equal to the current chat
    chat0=chat;
end

% Produce warning if maximum number of iterations is reached
if it==maxit
    warning('Maximum number of iterations reached in EM algorithm');
end

% -------------------------------------------------------------------------
% TRANSPOSE RETURN ARGUMENTS Fhat AND x2 
Fhat = Fhat';
x2 = x2'; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [chat,fhat,lambda,ss]=pc2(X,nfac)
% =========================================================================
% DESCRIPTION
% This function runs principal component analysis.
%
% -------------------------------------------------------------------------
% INPUTS
%           X      = dataset (one series per column)
%           nfac   = number of factors to be selected
%
% OUTPUTS
%           chat   = values of X predicted by the factors
%           fhat   = factors scaled by (1/sqrt(N)) where N is the number of
%                    series
%           lambda = factor loadings scaled by number of series
%           ss     = eigenvalues of X'*X 
%
% =========================================================================
% FUNCTION

% Number of series in X (i.e. number of columns)
N=size(X,2);

% Singular value decomposition: X'*X = U*S*V'
[U,S,~]=svd(X'*X);

% Factor loadings scaled by sqrt(N)
lambda=U(:,1:nfac)*sqrt(N);

% Factors scaled by 1/sqrt(N) (note that lambda is scaled by sqrt(N))
fhat=X*lambda/N;

% Estimate initial dataset X using the factors (note that U'=inv(U))
chat=fhat*lambda';

% Identify eigenvalues of X'*X
ss=diag(S);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function pos=minindc(x)
% =========================================================================
% DESCRIPTION
% This function finds the index of the minimum value for each column of a
% given matrix. The function assumes that the minimum value of each column
% occurs only once within that column. The function returns an error if
% this is not the case.
%
% -------------------------------------------------------------------------
% INPUT
%           x   = matrix 
%
% OUTPUT
%           pos = column vector with pos(i) containing the row number
%                 corresponding to the minimum value of x(:,i)
%
% =========================================================================
% FUNCTION

% Number of rows and columns of x
nrows=size(x,1);
ncols=size(x,2);

% Preallocate memory for output array
pos=zeros(ncols,1);

% Create column vector 1:nrows
seq=(1:nrows)';

% Find the index of the minimum value of each column in x
for i=1:ncols
    
    % Minimum value of column i
    min_i=min(x(:,i));
    
    % Column vector containing the row number corresponding to the minimum
    % value of x(:,i) in that row and zeros elsewhere
    colmin_i= seq.*((x(:,i)-min_i)==0);
    
    % Produce an error if the minimum value occurs more than once
    if sum(colmin_i>0)>1
        error('Minimum value occurs more than once.');
    end
    
    % Obtain the index of the minimum value by taking the sum of column
    % vector 'colmin_i'
    pos(i)=sum(colmin_i);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [x22,mut,sdt]=transform_data(x2,DEMEAN)
% =========================================================================
% DESCRIPTION
% This function transforms a given set of series based upon the input
% variable DEMEAN. The following transformations are possible:
%
%   1) No transformation.
%   
%   2) Each series is demeaned only (i.e. each series is rescaled to have a
%   mean of 0).
%   
%   3) Each series is demeaned and standardized (i.e. each series is
%   rescaled to have a mean of 0 and a standard deviation of 1).
%   
%   4) Each series is recursively demeaned and then standardized. For a
%   given series x(t), where t=1,...,T, the recursively demeaned series
%   x'(t) is calculated as x'(t) = x(t) - mean(x(1:t)). After the
%   recursively demeaned series x'(t) is calculated, it is standardized by
%   dividing x'(t) by the standard deviation of the original series x. Note
%   that this transformation does not rescale the original series to have a
%   specified mean or standard deviation.
%
% -------------------------------------------------------------------------
% INPUTS
%           x2      = set of series to be transformed (one series per
%                     column); no missing values;
%           DEMEAN  = an integer indicating the type of transformation
%                     performed on each series in x2; it can take on the
%                     following values:
%                           0 (no transformation)
%                           1 (demean only)
%                           2 (demean and standardize)
%                           3 (recursively demean and then standardize) 
%
% OUTPUTS
%           x22     = transformed dataset
%           mut     = matrix containing the values subtracted from x2
%                     during the transformation
%           sdt     = matrix containing the values that x2 was divided by
%                     during the transformation
%
% =========================================================================
% FUNCTION

% Number of observations in each series (i.e. number of rows in x2)
T=size(x2,1);

% Number of series (i.e. number of columns in x2)
N=size(x2,2);

% Perform transformation based on type determined by 'DEMEAN'
switch DEMEAN
    
    % ---------------------------------------------------------------------
    % No transformation
    case 0
        mut=repmat(zeros(1,N),T,1);
        sdt=repmat(ones(1,N),T,1);
        x22=x2;
        
    % ---------------------------------------------------------------------
    % Each series is demeaned only
    case 1
        mut=repmat(mean(x2),T,1);
        sdt=repmat(ones(1,N),T,1);
        x22=x2-mut;
        
    % ---------------------------------------------------------------------
    % Each series is demeaned and standardized 
    case 2
        mut=repmat(mean(x2),T,1);
        sdt=repmat(std(x2),T,1);
        x22=(x2-mut)./sdt;
        
    % ---------------------------------------------------------------------
    % Each series is recursively demeaned and then standardized
    case 3
        mut=NaN(size(x2));
        for t=1:T
            mut(t,:)=mean(x2(1:t,:),1);
        end
        sdt=repmat(std(x2),T,1);
        x22=(x2-mut)./sdt; 
end





