function [ celldata,meta, filepath,name,ext ] = readin_atf_acoustic(folder,numfile)
%READIN Reads in multiple Excel files 
%   Function to read in multiple Excel files from an inputted filepath 
%   into cell called "celldata" and return filepath, name of file and 
%   extension separately. 
%   Louisa Brotherson, 10/11/17
%    
%   Variables:-
%
%   celldata: cell array, stores atf data for each station
%   ext: string, file extensions of each file  
%   filepath: char, file path of each file
%   fullFileName: char, full path and filenames of all data to read in 
%   i: indice for "for" loop
%   name: char, name of file without extension or filepath                
%   theFiles: structure containing all Excel files in directory to read

    
% find all files that end in ".atf" 


if numfile == 1
    
    filePattern = fullfile(folder, '*_03.atf');  
    theFiles = dir(filePattern);
  
elseif numfile == 2
    
    filePattern = fullfile(folder, {'*_01.atf';'*_03.atf'});
    
for i =1:length(filePattern)

    theFiles{i} = dir(filePattern{i});

end

theFiles = cell2mat(theFiles);

theFiles = reshape(theFiles,[],1);
% [~,theFiles_ind]= unique([theFiles.bytes],'stable');
% theFiles = theFiles(theFiles_ind);    

elseif numfile == 3
    
    filePattern = fullfile(folder, {'*_01.atf';'*_03.atf';'*_05.atf';});
    
for i = 1:length(filePattern)

    theFiles{i} = dir(filePattern{i});

end

theFiles = cell2mat(theFiles);

theFiles = reshape(theFiles,[],1);
% [~,theFiles_ind]= unique([theFiles.bytes],'stable');
% theFiles = theFiles(theFiles_ind);    

elseif numfile == 4
    
    filePattern = fullfile(folder, {'*_01.atf';'*_03.atf';'*_05.atf';'*_07.atf'});
    
for i =1:length(filePattern)

    theFiles{i} = dir(filePattern{i});

end

theFiles = cell2mat(theFiles);

theFiles = reshape(theFiles,[],1);

end



for j = 1:numel(theFiles)
  %fn{j} = theFiles(j,1).name;
  % BE: on mac it seems I need the folder adding: 
  fn{j} = char(strcat(theFiles(j).folder, "/",theFiles(j).name));

  TF = contains(fn,'~$');
  if TF == 1
    delete(fn) 
  end
end

%for loop to create cell with geomagnetic annual mean data for each
%station

%BE - ^^ as far as I can see, this just reads the data, no mean calculated?
    
    for i = 1:length(theFiles)
    
        [filepath,name{i},ext] = fileparts(fn{i});
        celldata{i,1} = importdata(strcat(fn{i}));   
        
        %BE
        fid=fopen(strcat(filepath,"/",name{i},ext));
        nCols=9;
        format = repmat('%s', [1 nCols]);
        hdr{i}=textscan(fid, format, 1, 'HeaderLines', 1,'Delimiter', ';');
        for j=1:nCols
            meta{i,j}=split(char(hdr{i}{j}),'=');
        end

    end
    


end

