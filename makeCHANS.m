file_names = dir('./Data/*.mat');
files_shape = size(file_names);
num_files = files_shape(1);

for i = 1:num_files
    current_file = matfile(file_names(i).name);
    CHANNAMES = who(current_file);
    channames_size = size(CHANNAMES);
    n_chans = channames_size(1);
    CHANACTIVE = ones(n_chans,1);
    CHANSFILENAME = replace(file_names(i).name,'LFP.mat','CHANS.mat');
    save(string(CHANSFILENAME),"CHANNAMES","CHANACTIVE");
end
