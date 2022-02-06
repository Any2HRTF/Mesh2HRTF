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
%
% Author: Harald Ziegelwanger (Acoustics Research Institute, Austrian Academy of Sciences)
% Co-Author(s): Katharina Pollack (Acoustics Research Institute, Austrian Academy of Sciences)

function [data,frequency]=Output2HRTF_Load(foldername,filename, numFrequencies)
%   [data,frequency] = OUTPUT2HRTF_LOAD(foldername,filename)
%   
%   Loads results of the BEM-HRTF calculation from the NumCalc folder.
%
%   Input:
%       foldername ... path to NumCalc output
%       filename ..... either 'pBoundary' or 'pEvalGrid'
%
%   Output:
%       data ......... Matrix of complex values
%           dim1: frequencies
%           dim2: datapoints (complex pressure values)

%% --------------------check number of header lines------------------------
for ii=1:1000
    temp1=importdata([foldername, 'be.1', filesep, filename], ' ', ii);
    if ~isempty(temp1)
        if size(temp1.data,2)==3
            if temp1.data(2,1)-temp1.data(1,1)==1
                numHeaderlines_BE=ii;
                break
            end
        end
    end
end

%% -----------------------------load data----------------------------------
temp1 = importdata([foldername, 'be.1', filesep, filename], ' ', numHeaderlines_BE);

data=zeros(numFrequencies,size(temp1.data,1));
frequency=zeros(numFrequencies,1);

for ii=1:numFrequencies
    tmpData = importdata([foldername, 'be.', num2str(ii), filesep, filename], ' ', numHeaderlines_BE);
    if ~isempty(tmpData)
        tmpFrequency = importdata([foldername, '..', filesep, 'fe.out', filesep, 'fe.', num2str(ii), filesep, 'load'], ' ', 2);

        data(ii,:)=transpose(tmpData.data(:,2)+1i*tmpData.data(:,3));
        frequency(ii,1)=tmpFrequency.data(1,1);
    else
        data(ii,:)      = NaN;
        frequency(ii,1) = NaN;
    end
end

end %of function
