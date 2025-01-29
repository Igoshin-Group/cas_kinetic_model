function paths = get_all_files(folder)
search_path = fullfile(folder, '*.traces');
filenames = dir(search_path);
paths = {};
for ii=1:length(filenames)
    paths{ii} = fullfile(folder, filenames(ii).name);
end
end