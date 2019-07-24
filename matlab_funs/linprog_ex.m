function [outputArg1,outputArg2] = linprog_ex(inputArg1,inputArg2)
%LINPROG_EX Summary of this function goes here
%   Detailed explanation goes here
opt=optimoptions('linprog');
opt.Algorithm='dual-simplex';
opt.Display='off';


end

