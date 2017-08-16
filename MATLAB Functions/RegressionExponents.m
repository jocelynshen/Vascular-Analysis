function [a_REG,b_REG] = RegressionExponents(tips,rad,len,bad_ind)
%This function calculates the scaling exponents a and b using the
%regression method.

% Remove unwanted values from data;
tips(bad_ind)=[];
rad(bad_ind)=[];
len(bad_ind)=[];

% Take log of tip, radius, length data
tips_log=log(tips);
rad_log=log(rad);
len_log=log(len);

if length(rad_log)<=1
    a_REG=NaN;
else
    p_a=polyfit(tips_log,rad_log,1); % Linear fit of log-log plot
    a_REG=p_a(1);
end

if length(len_log)<=1
    b_REG=NaN;
else
    p_b=polyfit(tips_log,len_log,1); % Linear fit of log-log plot
    b_REG=p_b(1);
end
end