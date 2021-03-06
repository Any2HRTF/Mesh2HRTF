%                                Mesh2HRTF
%                Copyright (C) 2015 by Harald Ziegelwanger,
%        Acoustics Research Institute, Austrian Academy of Sciences
%                        mesh2hrtf.sourceforge.net
% 
% Mesh2HRTF is licensed under the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% Mesh2HRTF is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
% You should have received a copy of the GNU LesserGeneral Public License along with Mesh2HRTF. If not, see <http://www.gnu.org/licenses/lgpl.html>.
% 
% If you use Mesh2HRTF:
% - Provide credits:
%   "Mesh2HRTF, H. Ziegelwanger, ARI, OEAW (mesh2hrtf.sourceforge.net)"
% - In your publication, cite both articles:
%   [1] Ziegelwanger, H., Kreuzer, W., and Majdak, P. (2015). "Mesh2HRTF: Open-source software package for the numerical calculation of head-related transfer functions," in Proceedings of the 22nd ICSV, Florence, IT.
%   [2] Ziegelwanger, H., Majdak, P., and Kreuzer, W. (2015). "Numerical calculation of listener-specific head-related transfer functions and sound localization: Microphone model and mesh discretization," The Journal of the Acoustical Society of America, 138, 208-222.

function data=Output2HRTF_ReadBEMPerformance(filename)
%OUTPUT2HRTF_READBEMPERFORMANCE
%   []=Output2HRTF_ReadBEMPerformance(filename) reads
%   computation time data.
%
%   Input:
%       filename
%
%   Output:
%       data

fid=fopen(filename);
count=0;
while ~feof(fid)
    line=fgets(fid);
    if strfind(line,'Frequency')
        count=count+1;
        idx1=strfind(line,'Step');
        idx2=strfind(line,'Frequency');
        data(count,1)=sscanf(line(idx1+4:idx2-1),'%d');
        idx2=strfind(line,'=');
        data(count,2)=sscanf(line(idx2+1:end),'%d');
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
    if strfind(line,'Address computation ')
        break
    end
end

end %of function