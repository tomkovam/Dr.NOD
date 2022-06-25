function similarityPerSignature = computeSimilarityWithSignatures(tableSignaturesCOSMIC, normalisedMutations)

% Computes  distance of normalisedMutations to all signatures from tableSignaturesCOSMIC
% - tableSignaturesCOSMIC: columns represent signatures (each should sum to 1)
% - normalisedMutations: a column of 96 values (should sum to 1)

normalisedMutations = normalisedMutations/sum(normalisedMutations);

lstSignatures = tableSignaturesCOSMIC.Properties.VariableNames;
nSignatures = length(lstSignatures);

similarityPerSignature = NaN*ones(nSignatures, 1);

for iSignature = 1:nSignatures
    similarityPerSignature(iSignature) = cosineSimilarity(tableSignaturesCOSMIC.(lstSignatures{iSignature}), normalisedMutations);
end

