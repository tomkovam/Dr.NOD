function cosSim = cosineSimilarity(a,b)

cosSim = sum(a.*b)/sqrt(sum(a.^2)*sum(b.^2));            % 0.9436
