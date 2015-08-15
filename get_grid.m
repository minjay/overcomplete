all_files = dir('/Users/minjay/Documents/MATLAB/Needlets/overcomplete/SS24-Jun-2015/ss*');
all_names = {all_files.name};

n = length(all_names);
ss = cell(n, 1);

for i = 1:n
    i
    ss{i} = textread(all_names{i});
end
