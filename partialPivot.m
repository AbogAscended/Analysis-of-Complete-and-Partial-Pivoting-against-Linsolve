function [residualerror] = partialPivot(inA, inB)
inM = [inA inB];
[n,~] = size(inM);
r = 1:n;
for p = 1:n
    [~, mIndex] = max(abs(inM(r(p:end), p)));
    mIndex = mIndex + p - 1;
    r([p, mIndex]) = r([mIndex, p]);
    inM = GEGivesUsPartial(inM, r, p);
end
solV = backSubPartial(inM, r);
residualerror = norm(solV - linsolve(inA,inB))/norm(linsolve(inA,inB));
end

function[inM] = GEGivesUsPartial(inM, r, pivot)
[n,~] = size(inM);
for i = pivot+1:n
    inM(r(i), : ) = inM(r(i), : ) + -(inM(r(i),pivot)/inM(r(pivot),pivot)).*(inM(r(pivot),:));
    inM(r(i), pivot) = 0;
end
end

function [tempM] = backSubPartial(inM, r)
n = length(r);
tempM = zeros(n, 1);
for i = n:-1:1
    tempM(i) = (inM(r(i), end) - inM(r(i),1:n) * tempM) / inM(r(i), i);
end
end
