function data = write_output_report(filename)
%   data = write_output_report(filename)
%
%   Reads computation time data from NC.out file.
%
%   Input:
%       filename ... path and name of the output file from NumCalc 
%                    (usually some NC.out file)
%
%   Output:
%       data ....... computation time in seconds, with columns:
%                    'Frequency index', 'Frequency', 'Building', 'Solving',
%                    'Postprocessing', 'Total', 'relative error', 'iterations',
%                    'maximum number of iterations reached' (1 = yes, 0 = no)

% This file is part of the Mesh2HRTF software package developed by the
% Mesh2HRTF Developer Team (https://mesh2hrtf.org) and licensed under the 
% EUPL, Version 1.2, or, as soon as approved by the European Commission, 
% subsequent versions of the EUPL. Details on the license can be found 
% in the file "license.txt" provided with Mesh2HRTF package
% or at https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12
%
% You may not use this work except in compliance with the license.
% Unless required by applicable law or agreed to in writing, software 
% distributed under the license is distributed on an "AS IS" basis,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

% #Author: Fabian Brinkmann (TU-Berlin): 2022, original implementation
% #Author: Katharina Pollack (ARI, ÖAW): 2022, various improvements
% #Author: Piotr Majdak (ARI, ÖAW): 2023, help text, license boiler plate


fid=fopen(filename);
count=0;
data=zeros(1,9);
while ~feof(fid)
    line=fgets(fid);
    if strfind(line,'Frequency')
        count=count+1;
        idx1=strfind(line,'Step');
        idx2=strfind(line,'Frequency');
        data(count,1)=sscanf(line(idx1+4:idx2-1),'%d');
        idx2=strfind(line,'=');
        data(count,2)=sscanf(line(idx2+1:end),'%f');
    end
    if strfind(line,'Assembling the equation system  ')
        idx=strfind(line,':');
        data(count,3)=sscanf(line(idx+1:end),'%d');
    end
    if strfind(line,'Solving the equation system  ')
        idx=strfind(line,':');
        data(count,4)=sscanf(line(idx+1:end),'%d');
    end
    if strfind(line,'Post processing  ')
        idx=strfind(line,':');
        data(count,5)=sscanf(line(idx+1:end),'%d');
    end
    if strfind(line,'Total  ')
        idx=strfind(line,':');
        data(count,6)=sscanf(line(idx+1:end),'%d');
    end
    if strfind(line, 'relative error')
        [~, endIdx] = regexp(line, '(relative error = )');
        data(count,7)=sscanf(line(endIdx+1:end-1), '%f');
    end
    if contains(line, 'iterations') && ~contains(line, 'Warning')
        [startIdx, endIdx] = regexp(line, '(number of iterations = )\S+[,]');
        data(count,8) = sscanf(line(startIdx+23:endIdx-1), '%d');
    end
%     if ~isempty(regexp(line, '^\d+\s\d+'))
%         [idx1, idx2] = regexp(line, '\s');
%         data(count,7) = sscanf(line(idx2+1:end-1), '%f');
%         data(count,8) = sscanf(line(1:idx1-1), '%d');
%         data(count,9) = sscanf('0', '%d');
%     end
    if strfind(line, 'Warning: Maximum number of iterations')
        data(count,9) = sscanf('1', '%d');
    end
    if strfind(line,'Address computation ')
        fclose(fid);
        break
    end
end

end %of function