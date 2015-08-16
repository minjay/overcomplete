all_files = dir('/Users/minjay/Documents/MATLAB/Needlets/overcomplete/SS24-Jun-2015/ss*');
all_names = {all_files.name};

n = length(all_names);
ss = cell(n, 1);
degree_t = zeros(n, 1);

for i = 1:n
    i
    ss{i} = textread(all_names{i});
    degree_t(i) = str2num(all_names{i}(3:5));
end

save('/Users/minjay/Documents/MATLAB/Needlets/overcomplete/ss.mat', 'ss', 'degree_t')
