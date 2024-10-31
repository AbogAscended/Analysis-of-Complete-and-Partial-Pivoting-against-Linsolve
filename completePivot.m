function [residualerror] = completePivot(inA, inB)
inM = [inA inB];
[n,~] = size(inM);
r = 1:n;
c = 1:n+1;
for p = 1:n
    [newpR, newpC] =  findMaxComplete(inM,r,c,p);
    r([p,newpR]) = r([newpR,p]);
    c([p,newpC]) = c([newpC,p]);
    inM = GEGivesUsComplete(inM, r, c, p);
end
solV = backsub(inM, r, c);
residualerror = norm(solV - linsolve(inA,inB))/norm(linsolve(inA,inB));
end

function [logicalN, logicalM] = findMaxComplete(inM,r,c,p)
[n,~] = size(inM);
max = inM(r(p), r(p));
logicalN = p;
logicalM = p;
for i = p:n
    for j = p:n
        if abs(inM(i,j)) > abs(max)
            max = inM(r(i),c(j));
            logicalN = i;
            logicalM = j;
        end
    end
end
end

function[inM] = GEGivesUsComplete(inM, r, c, pivot)
[n,~] = size(inM);
for i = pivot+1:n
    inM(r(i), c) = inM(r(i), c) + -(inM(r(i),c(pivot))/inM(r(pivot),c(pivot))).*(inM(r(pivot),c));
    inM(r(i), c(pivot)) = 0;
end
end

function [outM] = backsub(inM, r, c)
n = length(r);
outM = zeros(n, 1);
tempM = zeros(n, 1);
for i = n:-1:1
    tempM(i) = (inM(r(i), end) - inM(r(i), c(1:n)) * tempM) / inM(r(i), c(i));
end
outM(c(1:n)) = tempM;
end