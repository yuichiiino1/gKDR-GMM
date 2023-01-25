function y = generatesalt(n, startframe, period)

%saltdata = zeros(n,1);
x = (1:n)-startframe;
y = (abs(sin(x/period*pi))).^0.25.*sign(sin(x/period*pi));
y(1:floor(startframe)) = 0; %%%%%%%%  2021/4/19  %%%%%%%%

%{
saltdata = -ones(n+freerun_length,1)/sqrt(2);  % mean=0, variance=1
salttable = readtable('stimulation_timing.csv','PreserveVariableNames',true);
%framesec = table2array(salttable(samplei,2));
startframe = table2array(salttable(samplei,3));
period = table2array(salttable(samplei,4));
from = startframe;
while round(from) <= n+freerun_length
    saltdata(round(from):round(min(n+freerun_length, from+period))) = 1/sqrt(2);
    from = from + period*2;
end
%}

end