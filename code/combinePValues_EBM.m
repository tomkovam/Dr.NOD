function pCombined = combinePValues_EBM(p1,p2)
% Expects column vectors of the same length
% p1 = pM; p2 = pE;

% Resources:
% https://github.com/broadinstitute/getzlab-PCAWG-pvalue_combination/tree/master/pvalue_combination
% https://github.com/IlyaLab/CombiningDependentPvaluesUsingEBM/blob/master/Matlab/EmpiricalBrownsMethod.m

nRows = length(p1);
pCombined = NaN*ones(nRows, 1);
pMatrix = [p1,p2];

isOK = sum(isnan(pMatrix), 2)==0; % We keep only rows without NaNs
pMatrix = pMatrix(isOK,:);
nRows = size(pMatrix, 1);
%%
covar_matrix = CalculateCovariances(pMatrix');
m = size(covar_matrix,1);
df_fisher = 2.0*m;
Expected = 2.0*m;
cov_sum = sum(sum(covar_matrix))-sum(diag(covar_matrix));
Var = 4.0*m+cov_sum;
c = Var/(2*Expected);
df_brown = 2*(Expected^2)/Var;
if df_brown > df_fisher
    df_brown = df_fisher;
    c = 1;
end

pCombined_clean = NaN*ones(nRows, 1);
for iRow = 1:nRows
    p_values = pMatrix(iRow,:);
    x = 2*sum(-log(p_values));
    p_brown = chi2cdf(1.0*x/c,df_brown,'upper');
    %p_fisher = chi2cdf(1.0*x,df_fisher,'upper');
    pCombined_clean(iRow) = p_brown;

end
pCombined(isOK) = pCombined_clean;

%% From https://github.com/IlyaLab/CombiningDependentPvaluesUsingEBM/blob/master/Matlab/EmpiricalBrownsMethod.m
%Calculate Covariances
    function covar_matrix = CalculateCovariances(data_matrix)
        m = size(data_matrix, 1);
        transformed_data_matrix = NaN(size(data_matrix));
        for f = 1:m
            transformed_data_matrix(f,:) = TransformData(data_matrix(f,:));
        end
        covar_matrix = cov(transformed_data_matrix',0,'omitrows');
    end

% Transform data
    function transformed_data_vector = TransformData(data_vector)
        isOK2 = ~isnan(data_vector);
        data_vector_clean = data_vector(isOK2);
        dvm = mean(data_vector_clean, 'omitnan');
        dvsd = std(data_vector_clean,1, 'omitnan');
        if (isnan(dvsd))
            transformed_data_vector = NaN*data_vector_clean;
        else
            if (dvsd==0)
                s = data_vector_clean;
            else
                s = (data_vector_clean-dvm)/dvsd;
            end
            [~,~,i] = unique(s);
            z = ecdf(s); z(1) = []; %ecdf
            transformed_data_vector = NaN*data_vector;
            transformed_data_vector(isOK2) = -2*log(z(i));
        end
    end
end