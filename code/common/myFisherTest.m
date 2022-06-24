function [p, enrichment, nObserved, nExpected] = myFisherTest(isGroup1, isGroup2, tail, verbose)

if (~exist('verbose', 'var'))
    verbose = false;
end

fisherTable = [sum(isGroup1 == 1 & isGroup2 == 1), sum(isGroup1 == 0 & isGroup2 == 1);
    sum(isGroup1 == 1 & isGroup2 == 0), sum(isGroup1 == 0 & isGroup2 == 0)];
[~,p,~] = fishertest(fisherTable, 'Tail', tail);

if (verbose)
    fisherTable
    fprintf('p-value %s\n', getPValueAsText(p));
end

nExpected = sum(isGroup1) * mean(isGroup2); % nSelectedGenes * fractionDriverGenesOutOfAllGenes
nObserved = sum(isGroup1 & isGroup2);
enrichment = nObserved/nExpected;


% https://systems.crump.ucla.edu/hypergeometric/index.php
% fisherTable = [18, 1016; 41, 14471]
% tmp = [18,1016+18; 18+41, 18+41+1016+14471]
% 18        1034
% 59       15546
% It gives identical results to the hypergeometric test
% Could we run a matlab function?
% hygecdf(x,M,K,N) % not sure what the parameters mean...