function similarityPerSignature = computeSimilarityWithSignatures(tableSignaturesCOSMIC, normalisedMutations)
% tableSignaturesCOSMIC: columns represent signatures (each should sum to 1)
% normalisedMutations: a column of 96 values (should sum to 1)
% we compute the distance of normalisedMutations to all signatures

% if (sum(normalisedMutations)~=1)
%     fprintf('We normalise normalisedMutations to sum to 1.\n');
    normalisedMutations = normalisedMutations/sum(normalisedMutations);
% end

lstSignatures = tableSignaturesCOSMIC.Properties.VariableNames;
nSignatures = length(lstSignatures);

similarityPerSignature = NaN*ones(nSignatures, 1);

for iSignature = 1:nSignatures
    similarityPerSignature(iSignature) = cosineSimilarity(tableSignaturesCOSMIC.(lstSignatures{iSignature}), normalisedMutations);
end

% cosSim = sum(a.*b)/sqrt(sum(a.^2)*sum(b.^2));            % 0.9436
