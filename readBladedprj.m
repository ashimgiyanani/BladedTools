% Script to read the Bladed xml files

clear; clc; close all
addpath('c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\Matlab\xml_io_tools_2010_11_05\');
addpath('c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\Matlab\f2matlab\')
filename = 'w:\Bladed\2Benergy_MC141\MC141_rev44_WG_061118_DoFactive_NoVar6p3.$PJ';

% read the Bladed prj file according to Matlab  xmlwrite function
DOMnode = xmlread(filename);
[s] = xml2struct(filename);

% read in the Bladed prj file
[tree, treeName] = xml_read(filename); % reads an xml file
disp([treeName{1} ' ='])  % shows the treename
gen_object_display(tree)  % displays the structure of the xml file
cdata = tree.BladedData.CDATA_SECTION;   % cdata section in an xml file, contains unsorted data
idx_returnkey = find(cdata==char(10));  % find the indices for end of line character (return key)
i_start = [1 idx_returnkey(1:end)+1];   % each line is starting at (position of return key + 1)
i_end = [idx_returnkey-1 size(cdata,2)]; % end position of the line i.e. (position of return key - 1)
idx_mend = [];
for ir = 1:numel(idx_returnkey)
  varLines{ir} = cdata(i_start(ir):i_end(ir)); % reading in lines
  if ~isempty(regexp(string(varLines{ir}),'\s+', 'match')) % if there is whitespace in varLines
    strSplits = regexp(string(varLines{ir}),'\s+', 'split');  % split the varLines
        if (~isempty(char(strSplits(1)))) && (isempty(str2num(strSplits(1))))
            varName(ir) = strSplits(1);  % first variable is varName
            varValue(ir) = strjoin(strSplits(2:end)); % second variable is varValue
        elseif isempty(char(strSplits(1))) || isempty(char(strSplits(3)))
            tempValue = [varValue(ir-1),strtrim(strjoin(strSplits))] ;
            varValue(ir) = strjoin(tempValue, ';');
        elseif (~isempty(str2num(strSplits(1))))
            varValue(ir) = strjoin(strSplits(1:end));
        end
  elseif ~isempty(regexp(string(varLines{ir}),'-\d', 'match')) % condition where '-1' is found
      tempValue = [varValue(ir-1), varLines{ir}];
      varValue(ir) = strjoin(tempValue,';')
  elseif ~isempty(regexp(string(varLines{ir}),'MEND', 'match')) % if no whitespace exists, find if there is 'MEND'
    idx_mend = [idx_mend; ir]; % find indiuces for module end
  elseif ~isempty(regexp(string(varLines{ir}),'\s+', 'split')) && (numel(regexp(string(varLines{ir}),'\s+', 'split')) == 1) && ~isempty(varLines{ir})
      strSplits = regexp(string(varLines{ir}),'\s+', 'split');
      varName(ir) = strSplits(1);
      varValue(ir) = '';
  end % if else condn
end % idx_returnkey condn
idx_mstart = find(varName == 'MSTART'); % find indices for module start
moduleNames = varValue(idx_mstart); %
idx_adat = find(moduleNames == 'ADAT');
moduleNames(idx_adat) = [];
idx_mstart(idx_adat) = [];
idx_mend(idx_adat) = [];

for im = 1:numel(moduleNames)
 for jm = 1:(idx_mend(im)-idx_mstart(im) -1)
  if ~exist(char(varName(idx_mstart(im)+jm)),'var'); % if variable not in workspace already
      if ~isempty(regexp(varName(idx_mstart(im)+jm), '^\d'))
         varName(idx_mstart(im)+jm) = strcat("No_",varName(idx_mstart(im)+jm));
      elseif ~isempty(regexp(varName(idx_mstart(im)+jm), '[\'']'))
         varName(idx_mstart(im)+jm) = regexprep(varName(idx_mstart(im)+jm), '[\'']', 'Ap');
      end
      
      eval([char(varName(idx_mstart(im)+jm)) ' = varValue(idx_mstart(im)+jm);']);
      if (ismissing(varName(idx_mstart(im)+jm+1)))
          varName(idx_mstart(im)+jm+1) = varName(idx_mstart(im)+jm);
      end
  elseif exist(char(varName(idx_mstart(im)+jm)),'var') && (~ismissing(varValue(idx_mstart(im)+jm))) % if variable exists and it's value is available on the subs lines
        sizeVar = numel(eval(char(varName(idx_mstart(im)+jm))));
        eval([char(varName(idx_mstart(im)+jm)) '{' num2str(sizeVar+1) '}'  ' = char(varValue(idx_mstart(im)+jm));']);
        if (ismissing(varName(idx_mstart(im)+jm+1)))
            varName(idx_mstart(im)+jm+1) = varName(idx_mstart(im)+jm);
        end
  end
 end % jm
end % im

%% Write the read xml file
com.mathworks.xml.XMLUtils.serializeXML(DOMnode,'c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\Matlab\MC141_rev44_WG_061118_DoFactive_NoVar6p3_xmlwrite.prj','ISO-8859-1'); 
filename1 = 'WindWise_2p3MW_141m_xml_write.prj';
DOMnode1 = xml_write(filename1, tree, treeName);
filename2 = 'WindWise_2p3MW_141m_xmlwrite.prj'
chr = xmlwrite(filename2, DOMnode);

%% write the cdata into a seperate file for reading in using f2matlab function
%formatSpec = repmat('%c',size(cdata));
%fid = fopen('example1.f90', 'w');
%fprintf(fid, formatSpec, cdata);
%fclose(fid);

%% read the fortran file
%filename = 'c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\Bladed\example1.f90';
%[filestr,numErrors,extraFunctions,localVar,varp,typeDefs]=f2matlab(filename);
%
%% get the individual lines from the saved fortran file
%ForM = 1;
%funstr=filestr2funstr_2(filestr,ForM);
%
%% check if first string exists or a number
%str1=breakOffFirstFormat(string(funstr(1)))
%
%
%%% Extras
%fid = fopen('Trial.ini', 'w');
%fprintf(fid, '%s', cdata);
%fclose(fid);
%
%% ini = IniConfig();
%% ini.ReadFile('Trial.ini')
%% ini.ToString()
%
%% ini = IniConfig();
%% ini.ReadFile('Trial.ini')
%% sections = ini.GetSections()
%% [keys, count_keys] = ini.GetKeys(sections{1})
%% values = ini.GetValues(sections{1}, keys)
%% new_values(:) = {rand()};
%% ini.SetValues(sections{1}, keys, new_values, '%.3f')
%% ini.WriteFile('example1.ini')
%
%% Creation of an example file called "Configuration.ini"
%% Example_write();
%
%% Initialization
%File = 'Trial.ini';
%I = INI('File',File);
%
%% INI file reading
%I.read();
%
%% Sections from INI file
%Sections = I.get('Sections');
%
%% Sections names
%fprintf(1,'Sections of "%s" file :\n',File);
%
%% Sections data
%for s = 1:numel(Sections)
    %fprintf(1,'- Section "%s"\n',Sections{s});
    %disp(I.get(Sections{s}));
%end

