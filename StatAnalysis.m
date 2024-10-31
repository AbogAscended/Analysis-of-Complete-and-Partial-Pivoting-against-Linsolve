function [output, time] = StatAnalysis ()
LNL = 11;
LNU = 1500;
SNLd = .01;
SNUd = 1;
SNLw = 1;
SNUw = 10;
MSL = 3;
MSU = 1000;
MA = 1000000;
%parallel computing workers 12
parpool(12);
%output all of the errors calculated for each matrix, each variabel in
%which differs only by method(complete/partial) were run on the same
%matrix per index i.
tic;

[whole_partial_small, whole_complete_small, decimal_partial_small, decimal_complete_small] = smallnumMats(MA, SNLd, SNUd, SNLw, SNUw, MSL, MSU);
output{1} = {whole_partial_small, whole_complete_small, decimal_partial_small, decimal_complete_small};
disp("Small Done");
[whole_partial_large, whole_complete_large, decimal_partial_large, decimal_complete_large] = largenumMats(MA, LNL, LNU, MSL, MSU);
output{2} = {whole_partial_large, whole_complete_large, decimal_partial_large, decimal_complete_large};
disp("Large Done");
[whole_partial_mixed, whole_complete_mixed, decimal_partial_mixed, decimal_complete_mixed] = mixednumMats(MA, LNL, SNLd, SNUd, SNLw, SNUw, LNU, MSL, MSU);
output{3} = {whole_partial_mixed, whole_complete_mixed, decimal_partial_mixed, decimal_complete_mixed};
disp("All Done");
time = toc;
delete(gcp('nocreate'));
end

%function to specifically calculate the small num matricies
function [whole_errPartial, whole_errComplete, decimal_errPartial, decimal_errComplete] = smallnumMats(MA, SNL, SNU, SNLw, SNUw,  MSL, MSU)
whole_errPartial = [];
whole_errComplete = [];
decimal_errPartial = [];
decimal_errComplete = [];
parfor i = 1:MA
    rng(i + rand());

    matsize = randi([MSL MSU],1);

    A = randi([SNLw SNUw], matsize);
    B = randi([SNLw SNUw], matsize, 1);

    whole_errPartial(i) = partialPivot(A,B);
    whole_errComplete(i) = completePivot(A,B);

    A = SNL + (SNU - SNL) * rand(matsize);
    B = SNL + (SNU - SNL) * rand(matsize,1);

    decimal_errPartial(i) = partialPivot(A,B);
    decimal_errComplete(i) = completePivot(A,B);
end
end
%function to calculate specifically the large num matricies
function [whole_errPartial, whole_errComplete, decimal_errPartial, decimal_errComplete] = largenumMats(MA, LNL, LNU, MSL, MSU)
whole_errPartial = [];
whole_errComplete = [];
decimal_errPartial = [];
decimal_errComplete = [];
parfor i = 1:MA
    rng(i + rand());
    matsize = randi([MSL,MSU]);

    A = randi([LNL LNU], matsize, matsize);
    B = randi([LNL LNU], matsize, 1);

    whole_errPartial(i) = partialPivot(A,B);
    whole_errComplete(i) = completePivot(A,B);

    A = LNL + (LNU - LNL) * rand(matsize);
    B = LNL + (LNU - LNL) * rand(matsize,1);

    decimal_errPartial(i) = partialPivot(A,B);
    decimal_errComplete(i) = completePivot(A,B);
end
end
%this one was the hardest to figure out how to program as it required
%writing randomaugmat to make it reasonably sized.
function [whole_errPartial, whole_errComplete, decimal_errPartial, decimal_errComplete] = mixednumMats(MA, LNL, SNLd, SNUd, SNLw, SNUw, LNU, MSL, MSU)
whole_errPartial = [];
whole_errComplete = [];
decimal_errPartial = [];
decimal_errComplete = [];
parfor i = 1:MA
    rng(i + rand());

    n = randi([MSL MSU]);

    [A,B] = randomAugMat(n, 1, SNUw, SNLw, LNU, LNL);

    whole_errComplete(i) = completePivot(A,B);
    whole_errPartial(i) = partialPivot(A,B);

    [A,B] = randomAugMat(n, 0, SNUd, SNLd, LNU, LNL);

    decimal_errComplete(i) = completePivot(A, B);
    decimal_errPartial(i) = partialPivot(A, B);
end
end

function [A, B] = randomAugMat(n, whole, Supper, Slower, Lupper, Llower)
%generates random 0,1 uniform matrix for main and augment vector
A = rand(n);
B = rand(n,1);
%make new logical matrix and vector based on 5050 of being greater than
%.5 in order to randomly choose big or small numbers.
Abool = A > .5;
Bbool = B > .5;
if ~whole
    %make iterator list of indexes for large and small numbers
    lIndex = find(Abool);
    sIndex = find(~Abool);
    %actually generate numbers using rand and the formula for generating
    %a continuous distribution between wanted integers
    lNums = (rand(size(lIndex)) * (Lupper - Llower)) + Llower;
    sNums = (rand(size(sIndex)) * (Supper - Slower)) + Slower;
    %actually feel the numbers into the matrix we orginally created
    A(lIndex) = lNums;
    A(sIndex) = sNums;
    %same proccess for b vector
    lIndex = find(Bbool);
    sIndex = find(~Bbool);
    lNums = (rand(size(lIndex)) * (Lupper - Llower)) + Llower;
    sNums = (rand(size(sIndex)) * (Supper - Slower)) + Slower;
    B(lIndex) = lNums;
    B(sIndex) = sNums;
else
    %same thing but for whole numbers, slightly easier to do.
    lIndex = find(Abool);
    sIndex = find(~Abool);
    lNums = randi([Llower Lupper], size(lIndex));
    sNums = randi([Slower Supper], size(sIndex));
    A(lIndex) = lNums;
    A(sIndex) = sNums;
    lIndex = find(Bbool);
    sIndex = find(~Bbool);
    lNums = randi([Llower Lupper], size(lIndex, 1), 1);
    sNums = randi([Slower Supper], size(sIndex, 1), 1);
    B(lIndex) = lNums;
    B(sIndex) = sNums;
end
end



