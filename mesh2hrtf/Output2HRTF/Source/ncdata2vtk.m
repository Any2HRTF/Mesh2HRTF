function Ncdata2vtk(name,datatype,rootdir,datadir1, ...
                                  meshdir, evaldir)
% exports data from NumCalc results to vtk files for display with paraview
% 
%
% input: 
%   name: Name of the vtk file(s) to be produced
%   datatype: type of data to be presented
%       0    pressure boundary
%       1    part. velocity boundary
%       2    pressure evalgrid
%       3    part. velocity evalgrid
%       4    pressure everywhere
%       5    part. velocity everywhere
%       6    display all (default)
%   rootdir: name of the root directory containing the be.out file,
%        default the current directory
%   datadir: name of the directory containing the data: default
%            be.out
%   meshdir: directory that contains the mesh data. Default 
%            rootdir/../../ObjectMeshes/Reference/
%   evaldir: directory that contains the mesh data of the
%            evalgrid. Default: rootdir/../../EvaluationGrids/3_ARI/
%    
% the data will be stored in the vtk file as Amplitude and Phase values
%
% written by kreiza, if you have any questions, suggestions 
%   send me an email, depending on my workload I will try to reply
%   in time: wolfgang.kreuzer@oeaw.ac.at
%
% if you start this script in the CPU_*_Core_* directory everything
% should work without the need for input parameter, (see below for
% the default input parameters)
%
%
%   This script is distributed in the hope that it will be useful,              
%   but WITHOUT ANY WARRANTY; without even the implied warranty of            
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

if(nargin < 4) 
    datadir1 = 'be.out';
    if(nargin < 3)
        rootdir = '.';
        if(nargin < 2)
            datatype = 6;   
            if(nargin < 1)
                name = 'output';
            end
        end
    end
end

if(nargin < 6) 
    evaldir = sprintf('%s/../../EvaluationGrids/3_ARI', ...
                      rootdir);
    if(nargin < 5)
        meshdir = sprintf('%s/../../ObjectMeshes/Reference', ...
                          rootdir);
    end
end

% check consistency between datatype and evalname
if( isempty(evaldir) && (datatype > 1) )
    datatype = 1;
end

if( isempty(meshdir) && (datatype < 2 || datatype > 3) )
    datatype = 3;
end

%% read the grid data

if(~isempty(meshdir))
    filename = sprintf('%s/Nodes.txt', ...
                       meshdir);
    fh = fopen(filename,'r');
    if(fh == -1)
        error('Sorry cannot open the objects node file');
    end
    
    nodes = dlmread(filename);  % matlab cannot read fh, octave can
    nnodes = nodes(1,1);    % first entry should be the number of
                            % nodes
    nodes(1,:) = [];
    fclose(fh);
    %  read the elements
    filename = sprintf('%s/Elements.txt', ...
                       meshdir);
    % I wish it would be as easy as above, but NumCalc allows us
    % the use of triangles and quadrilaterals, so we may have different
    % numbers of columns for different elements and filling those up
    % with 0s may lead to wrong results, which means we have to
    % parse the file, see read_elements
    fh = fopen(filename,'r');
    if(fh == -1)
        error('Sorry cannot open the objects element file');
    end
    [Elem3,Elem4] = read_elements(fh);
    fclose(fh);
else
    Elem3 = []; Elem4 = []; nodes = []; nnodes = 0;
end

% Do the same with the Evalgrid
if(~isempty(evaldir))
    filename = sprintf('%s/Nodes.txt', evaldir);
    fh = fopen(filename,'r');
    if(fh == -1)
        error('Sorry cannot open the objects node file');
    end
    
    enodes = dlmread(filename); % matlab can only read filenames
                                % not handles
    ennodes = enodes(1,1);    % first entry should be the number of
                             % nodes
    enodes(1,:) = [];
    fclose(fh);
    

    filename = sprintf('%s/Elements.txt', evaldir);
    fh = fopen(filename,'r');
    if(fh == -1)
        error('Sorry cannot open the objects element file');
    end
    
    [EElem3,EElem4] = read_elements(fh);
    fclose(fh);
else
    EElem3 = []; EElem4 = []; ennodes = 0; enodes = [];
end

% NumCalc would allow the nodes not to be ordered and to start with
% an index different to 0 (which is used for the evaluation grid)
% theoretically NumCalc would also allow to have jumps in the
% indices (as it has been used for the nodes between eval and
% object boundary), however, index jumps of eval nodes or boundary
% nodes will not be considered here as a possibility, which means
% that the results for meshes with index jumps will be diplayed incorrectly
%
% In case of jumps (which may be given, if you build up your mesh
% of different parts) the script needs to be extended, or
% even better NumCalc should be rewritten, to produce a better
% output file, if one would want to extend the script one possible
% way would be to compare the ordered indices with [1:nnodes] and
% if jumps are detected change Elem* accordingly, a second option
% would be a find and replace routine that replaces the index of
% the node with its actual position in the ordered list, but a bit
% of care has to be taken here, because you would not want the
% overwrite already changed nodes
% please remember, if you use all the mesh2hrtf routines
% starting from a blender mesh, this should all be *no problem*,
% everything is as it should be, however if you use your own
% created input files, you may want to be careful
%
% Things to do for the future: From the NumCalc side it would not
% be necessary to have the nodes order, because they can be
% identified via its node index, however this has not been
% considered here, it would make things a bit difficult for the vtk file


% nodes in vtk files do not have explicit node numbers and are
% assumed to be sorted and starting with 0 

if(nnodes > 0) 
    [nodenum, order] = sort(nodes(:,1));
    nodes = nodes(order,2:4);
    if( size(Elem3,1) > 0)
        Elem3 = Elem3(:,1:4);
        Elem3(:,2:4) = Elem3(:,2:4) - nodenum(1);     %% vtk starts with 0 
    end
    
    if( size(Elem4,1) > 0)
        Elem4 = Elem4(:,1:5);
        Elem4(:,2:5) = Elem4(:,2:5) - nodenum(1);     %% vtk starts with 0 
    end
end

%do the same for the evalgrid, now the nodes definitely do not
%start with 0

if(ennodes > 0)
    [nodenum, order] = sort(enodes(:,1));
    enodes = enodes(order,2:4);

    if(size(EElem3,1) > 0)
        EElem3 = EElem3(:,1:4);
        EElem3(:,2:4) = EElem3(:,2:4) - nodenum(1);     %% vtk starts with 0 
    end
    
    if(size(EElem4,1) > 0)
        EElem4 = EElem4(:,1:5);
        EElem4(:,2:5) = EElem4(:,2:5) - nodenum(1);     %% vtk starts with 0 
    end
end

% get the number of the frequency steps
% if you know the relation between frequency and frequency step you
% can create a time filter in paraview for displaying the right
% frequency, just a remark

% nf = length(ls(sprintf('%s/be.out',rootdir))); works only in
% octave
nf = length( dir( sprintf('%s/%s',rootdir,datadir1) ) ) - 2;


printpb = 0;
printpe = 0;
printvb = 0;
printve = 0;
pB = [];
pE = [];
vB = [];
vE = [];
%
% we look at the data for all frequency steps so one can make an
% animation out of the data, if you want things just for one
% specific frequency you'll need to change nf below
%
for i = 1:nf,
    if( ismember(datatype,[0,4,6]) )  
        %% read the pressure data at the boundary
        filename = sprintf('%s/%s/be.%d/pBoundary',rootdir, ...
                           datadir1, i);
        % now that is a bit spooky because of the string in file
        % but it seems to work, however just on octave
        %pB = dlmread(filename);
        %pB(1:3,:) = [];
        printpb = 1;
        % should work in both matlab and octave
        % throw away the first 3 lines

        D = importdata(filename, ' ', 3);
        pB = D.data;
        pB(:,2) = pB(:,2) + pB(:,3)*1.0i;
        pB(:,3) = [];

        % please NOTE: It is assumed that the Elements are ordered
        % and have no index jumps
    end
    if( ismember(datatype,[2,4,6]) )
        %% read the pressure data at the eval Grid
        filename = sprintf('%s/%s/be.%d/pEvalGrid',rootdir, ...
                            datadir1, i);
        %        pE = dlmread(filename2);
        %pE(1:3,:) = [];
        %pE(:,2) = pE(:,2) + pE(:,3)*1.0i;
        %pE(:,3) = [];
        D = importdata(filename, ' ', 3);
        pE = D.data;
        pE(:,2) = pE(:,2) + pE(:,3)*1.0i;
        pE(:,3) = [];
        printpe = 1;
    end
    if( ismember(datatype, [1,5,6]) )
        %% read the velocity at the boundary
        filename = sprintf('%s/%s/be.%d/vBoundary',rootdir, ...
                            datadir1, i);
        %vB = dlmread(filename3);
        D = importdata(filename, ' ', 3);
        vB = D.data;
        %vB(1:3,:) = [];
        vB(:,2) = vB(:,2) + vB(:,3)*1.0i;
        printvb = 1;

    end
    if( ismember(datatype, [2,5,6]) )
        %% read the velocity at the eval grid
        filename = sprintf('%s/%s/be.%d/vEvalGrid',rootdir, ...
                            datadir1, i);
        %        vE = dlmread(filename4);
        
        %vE(1:3,:) = [];
        D = importdata(filename, ' ', 3);
        vE = D.data;
        vE(:,2) = vE(:,2) + 1.0i*vE(:,3);
        vE(:,3) = [];
        printve = 1;
    end
    
    % the value at the boundary is given at the collocations nodes which is the
    % midpoint of each element, however in order to get a smoother
    % graph we decided to assign a value to each nodes of the grid
    % if you do not like the idea, just reformulate the vtk file
    % with cell data instead of point data, this should work too
    if(printpb)
        pBound = zeros(nnodes,1);
    end
    if(printvb) 
        vBound = zeros(nnodes,1);
    end
    for j = 1:nnodes, 
        if( ~isempty(Elem3) )
            % find the elements that contain the node
            [n31,n32] = find( j-1 == Elem3(:,2:4) );
        else
            n31 = 0;
        end

        if( ~isempty(Elem4) )
            [n41,n42] = find( j-1 == Elem4(:,2:5) );
        else
            n41 = 0;
        end
        n = length(n41) + length(n31);
        if( n31 )
            if(printpb)
                % take the mean value of all elements ==
                % collocation nodes containing the vertex
                % 
                % the index juggling may seem a bit unnecessary
                % right now, but remember, it is possible to have
                % jumps in the node numbers and element numbers so
                % this can be helpful if the script will be
                % extended for that case
                for j1 = 1:length(n31),
                    pBound(j) = pBound(j) + pB( find( pB(:,1) == Elem3(n31(j1),1) ...
                                                  ), 2);
                end
            end
            if(printvb)
                for j1 = 1:length(n31),
                    vBound(j) = vBound(j) + vB( find( vB(:,1) == Elem3(n31(j1),1) ...
                                                      ), 2);
                end
            end 
        end
        if( n41 )
            if(printpb)
                for j1 = 1: length(n41),
                    pBound(j) = pBound(j) + pB( find( pB(:,1) == Elem4(n41(j1),1) ...
                                                      ), 2);
                end
            end
            if(printvb)
                for j1 = 1:length(n41),
                    vBound(j) = vBound(j) + vB( find( vB(:,1) == Elem4(n41(j1),1) ...
                                                      ), 2);
                end
            end 
            
        end
        if(printpb)
            pBound(j) = pBound(j) / n;
        end
        if(printvb)
            vBound(j) = vBound(j) / n;
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %          write the vtk-files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    filename = sprintf('%s_bound_%d.vtk', name, i);
    if(printpb + printvb == 2)
        writefile(filename, nodes, Elem3, Elem4, pBound, vBound);
    else
        if(printpb)
            writefile(filename, nodes, Elem3, Elem4, pBound);
        end
        if(printvb)
            writefile(filename, nodes, Elem3, Elem4, vBound);
        end
    end
    
    if(~isempty(pE))
        pE(:,1) = [];
    end
    if(~isempty(vE))
        vE(:,1) = [];
    end
    
    filename = sprintf('%s_eval_%d.vtk', name, i);
    if(printpe + printve == 2)
        writefile(filename, enodes, EElem3, EElem4, pE, vE);
    else
        if(printpe)
            writefile(filename, enodes, EElem3, EElem4, pE);
        end
        if(printve)
            writefile(filename, enodes, EElem3, EElem4, vE);
        end
    end
end  % loop over freqsteps



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                 Function writefile
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writefile(filename, N, E3, E4, data1, data2)    
%
% write the mesh data and the values at the mesh nodes to a vtk file
% 
% input: filename  : name of the output file
%        N         : list of nodes (matrix nx3)
%        E3        : list of all triangular elements (matrix mx3)
%        E4        : list of all quadrilateral elements (matrix mx4)
%        data1     : vector containing the comlex pressure or the
%                    complex velocity
%        data2     : optional data vector if both pressure and
%                    velocity shall be displayed
if(nargin < 6)
    data2 = [];
end
% E3 and E4 still contain the indices of the elements
% get rid of them

if(~isempty(E3))
    E3(:,1) = [];
end
if(~isempty(E4))
    E4(:,1) = [];
end

nnodes = size(N,1);

fh = fopen(filename,'w');
if (fh == -1) 
    error('Sorry cound not open the output file.')
end
%
%  write the header
%
fprintf(fh, '# vtk DataFile Version 3.0\n');
fprintf(fh, 'Generated by data2vtk\n');
fprintf(fh, 'ASCII\n');
%    fprintf(fh, '\n');
fprintf(fh, 'DATASET POLYDATA\n');
fprintf(fh, 'POINTS %d double\n',nnodes);
fprintf(fh, '%e %e %e\n', N');
%
%  write the triangle geo data
%

fprintf(fh, 'POLYGONS %d %d\n',size(E3,1)+size(E4,1), ...
        size(E3,1)*4+size(E4,1)*5);

if(size(E3,1) > 0)
    fprintf(fh, '3 %d %d %d\n', E3');
end
%
% write 4-sided geometry
%
if(size(E4,1) > 0)
    fprintf(fh, '4 %d %d %d %d\n', E4');
end

fprintf(fh, 'POINT_DATA %d\n', size(N,1) );
fprintf(fh, 'SCALARS Amplitude1(dB) double\n');
fprintf(fh, 'LOOKUP_TABLE default\n');
data = 20*log10(abs(data1)/2e-5);
fprintf(fh, '%e\n', data);

fprintf(fh, 'SCALARS Phase1(rad) double\n');
fprintf(fh, 'LOOKUP_TABLE default\n');
data = angle(data1);
fprintf(fh, '%e\n', data);

if(~isempty(data2))
    fprintf(fh, 'SCALARS Amplitude2(dB) double\n');
    fprintf(fh, 'LOOKUP_TABLE default\n');
    i = find(abs(data2) < 1e-10);
    data2(i) = 1e-10;
    data = 20*log10(abs(data2)/2e-5);
    fprintf(fh, '%e\n', data);

    fprintf(fh, 'SCALARS Phase2(rad) double\n');
    fprintf(fh, 'LOOKUP_TABLE default\n');
    data = angle(data2);
    fprintf(fh, '%e\n', data);
end
fclose(fh);
% function end





function [Elem3,Elem4] = read_elements(fh)
% reads element number and element nodes from filehandle
% assumes mesh2hrtf format
nelems = str2num(fgets(fh));
Elem3 = [];
Elem4 = [];
for j = 1:nelems
    str = fgets(fh);
    el = sscanf(str,'%f');
    n = length(el);
    if(n == 7)
        % we have a triangle
        el = el(1:4);
        Elem3 = [Elem3;el'];  
    else
        % quadrilateral
        el = el(1:5);
        Elem4 = [Elem4;el'];  
    end
end

