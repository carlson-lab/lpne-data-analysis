function saveJson(saveFile, labels)

data = jsonencode(labels);
saveJsonFile=strrep(saveFile,'.mat','.json');
% Create a new file and print json-encoded data to it
filename = sprintf(saveJsonFile);
fID = fopen(filename,'w+');
fprintf(fID,'%s',data);

end