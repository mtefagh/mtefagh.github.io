clear all; close all; clc;
cvx_solver gurobi;
cvx_quiet('true');

opts = delimitedTextImportOptions("NumVariables", 24);
opts.VariableTypes = repmat("double", 1, 24);
opts = setvaropts(opts, "FillValue", 0);
grade = readtable("Convex Optimization 98-2 - grade.csv", opts);

grade = grade(2:45,:);
homework =  table2array(grade(:, 18))./160;
project =  (table2array(grade(:, 19))+table2array(grade(:, 20)))./20;
final =  table2array(grade(:, 21))./70;
extra = table2array(grade(:, 22));
class = table2array(grade(:, 23));
grade = table2array(grade(:, [1, 2, 24]));

figure;
histogram(grade(:, 3), (0:40)./2);
savefig('before.fig');

for i = 1:length(final)
    grade(i, 2) = min(round(score(homework(i), project(i), final(i), extra(i), class(i)), 3), 20);
end

figure;
histogram(grade(:, 2), (0:40)./2);
savefig('after.fig');

writematrix(grade(:, 1:2), 'grades.csv', 'Delimiter', ',');

function Gstar = score(h, p, f, e, c)
    cvx_begin
    variables xh xp xf a
    maximize(xh*h + xp*p + xf*f + (2-a)*max(c, h) + e)
    subject to
        xh + xp + xf == 20
        xh >= 4*inv_pos(a)
        xp >= 2*inv_pos(a)
        xf >= 14*inv_pos(a)
        a <= 0.3*geo_mean([xh, xp, xf])
    cvx_end
    Gstar = cvx_optval;
end