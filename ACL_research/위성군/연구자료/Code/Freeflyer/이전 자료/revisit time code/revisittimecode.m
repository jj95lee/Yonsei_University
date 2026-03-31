clc
clear all

a = importdata("point 4x35.txt");

b = mode(a.data, 'all')

c = sum(a.data == b)

tbl = tabulate(a.data);

d = [];
for i = 1:length(tbl(:, 1))
    if tbl(i, 2) == 1
        d = [d; tbl(i, 1)];
    end
end
