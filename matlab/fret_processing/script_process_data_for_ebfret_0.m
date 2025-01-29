% Script to pre-process smFRET data for the ebFRET tool
% Date: 2021-07-07
% Updated: 2021-07-07

clc;
clear;
close all;

% wild-type data
wt_0_path = '../../data/fret/raw/wt/wt-0bp/';
wt_3_path = '../../data/fret/raw/wt/wt-3bp/';
wt_4_path = '../../data/fret/raw/wt/wt-4bp/';

% hf1-cas9 data
hf1_0_path = '../../data/fret/raw/hf1/hf1-0bp/';
hf1_1_path = '../../data/fret/raw/hf1/hf1-1bp/';
hf1_2_path = '../../data/fret/raw/hf1/hf1-2bp/';
hf1_3_path = '../../data/fret/raw/hf1/hf1-3bp/';
hf1_4_path = '../../data/fret/raw/hf1/hf1-4bp/';

% eSpy-cas9 data
eSpy_0_path = '../../data/fret/raw/enh/enh-0bp/';
eSpy_1_path = '../../data/fret/raw/enh/enh-1bp/';
eSpy_2_path = '../../data/fret/raw/enh/enh-2bp/';
eSpy_3_path = '../../data/fret/raw/enh/enh-3bp/';
eSpy_4_path = '../../data/fret/raw/enh/enh-4bp/';


% output paths
output_folder = '../../data/fret/preprocessed/wt/';
status = mkdir(output_folder);

wt_0_out = fullfile(output_folder, 'wt-0bp_preprocessed.dat');
wt_3_out = fullfile(output_folder, 'wt-3bp_preprocessed.dat');
wt_4_out = fullfile(output_folder, 'wt-4bp_preprocessed.dat');

output_folder = '../../data/fret/preprocessed/hf1/';
status = mkdir(output_folder);

hf1_0_out = fullfile(output_folder, 'hf1-0bp_preprocessed.dat');
hf1_1_out = fullfile(output_folder, 'hf1-1bp_preprocessed.dat');
hf1_2_out = fullfile(output_folder, 'hf1-2bp_preprocessed.dat');
hf1_3_out = fullfile(output_folder, 'hf1-3bp_preprocessed.dat');
hf1_4_out = fullfile(output_folder, 'hf1-4bp_preprocessed.dat');

output_folder = '../../data/fret/preprocessed/enh/';
status = mkdir(output_folder);

eSpy_0_out = fullfile(output_folder, 'enh-0bp_preprocessed.dat');
eSpy_1_out = fullfile(output_folder, 'enh-1bp_preprocessed.dat');
eSpy_2_out = fullfile(output_folder, 'enh-2bp_preprocessed.dat');
eSpy_3_out = fullfile(output_folder, 'enh-3bp_preprocessed.dat');
eSpy_4_out = fullfile(output_folder, 'enh-4bp_preprocessed.dat');

% create files
process_folder(wt_0_out, wt_0_path);
process_folder(wt_3_out, wt_3_path);
process_folder(wt_4_out, wt_4_path);

process_folder(hf1_0_out, hf1_0_path);
process_folder(hf1_1_out, hf1_1_path);
process_folder(hf1_2_out, hf1_2_path);
process_folder(hf1_3_out, hf1_3_path);
process_folder(hf1_4_out, hf1_4_path);

process_folder(eSpy_0_out, eSpy_0_path);
process_folder(eSpy_1_out, eSpy_1_path);
process_folder(eSpy_2_out, eSpy_2_path);
process_folder(eSpy_3_out, eSpy_3_path);
process_folder(eSpy_4_out, eSpy_4_path);

function process_folder(output, input)
% get paths for .traces files
paths = get_all_files(input);

% load all data
for ii = 1:length(paths)
    path = paths{ii};
    datasets(ii) = load_binary_data(path);
end

% merge datasets
dataset = merge_datasets(datasets);

% pre-process data
dataset = smooth_dataset(dataset);
dataset = remove_photobleach(dataset);
dataset = filter_dataset(dataset);

% export data
export_ascii_data(dataset, output);
end