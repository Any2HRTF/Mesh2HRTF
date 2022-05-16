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
