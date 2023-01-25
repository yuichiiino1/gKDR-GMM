function d2 = distfun(XI,XJ)
% consider NaN  2017/6/21
d2 = [];
for i=1:size(XJ,1)
    mask = ~isnan(XI) & ~isnan(XJ(i,:));
    R = corrcoef(XI(mask),XJ(i,mask));
    d2 = [d2; 1-R(1,2)];
end
end