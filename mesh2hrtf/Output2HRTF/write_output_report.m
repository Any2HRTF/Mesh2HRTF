function data = write_output_report(filename)
%   [] = Output2HRTF_ReadBEMPerformance(filename)
%
%   Reads computation time data from NC.out file.
%
%   Input:
%       filename ... read computation time from NC.out (output file from NumCalc)
%
%   Output:
%       data ....... computation time in seconds, with columns:
%                    'Frequency index', 'Frequency', 'Building', 'Solving',
%                    'Postprocessing', 'Total', 'relative error', 'iterations',
%                    'maximum number of iterations reached' (1 = yes, 0 = no)

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