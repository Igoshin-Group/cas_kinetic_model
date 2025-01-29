function merged = merge_datasets(datasets)
merged = datasets(1);
for ii = 2:length(datasets)
    cur_dataset = datasets(ii);
    for jj = 1:length(cur_dataset.data)
        merged.data(end+1) = cur_dataset.data(jj);
    end
end
end