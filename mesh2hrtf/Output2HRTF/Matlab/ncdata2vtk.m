function ncdata2vtk(name,datatype,rootdir,datadir1, ...
                                  meshdir, evaldir)
% exports data from NumCalc results to vtk files for display with paraview
%
% we have to distinguish between two types of data
% a) data on the scatterer. For that the values for the BEM nodes are
% set by interpolating the data of all elements containing this
% node. Alternatively one could use Cell_DATA entry to display the
% constant data on each element, lets call this version elementwise =
% 1, default is elementwise = 0, but for large meshes this may take
% lot of time

% b) Evalnodes. This time the value is given on the specific node, and
%    we can use POINT_DATA
%
% Warning: This script assumed that the information about nodes is
%          given in the file Nodes.txt,
%         information about elements is given in Elements.txt
%
% Additional Warning: It is assumed that there are no jumps and holes
%    in the nodenumbering and the elementnumbering for bem and eval
%    mesh, thus nojumps = 1
%
% In the future, we will generalize these names and set them as input 
  
% input: 
%   name: Name of the vtk file(s) to be produced
%   datatype: type of data to be presented
%       0    pressure boundary
%       1    part. velocity boundary
%       2    pressure evalgrid
%       3    part. velocity evalgrid
%       4    pressure everywhere (default)
%       5    part. velocity everywhere
%       6    display all 
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
% written by Wolfgang Kreuzer, if you have any questions, suggestions 
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
%
% Licensed under the EUPL, Version 1.2 or – as soon they will be
%   approved by the European Commission – subsequent versions of the
%   EUPL (the "Licence"); you may not use this work except in
%   compliance   with the Licence. You may obtain a copy of the
%   Licence at:   http://joinup.ec.europa.eu/software/page/eupl
%
% Unless required by applicable law or agreed to in writing,
%   software distributed under the licence is distributed on an "AS IS"
%   basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
%   express or implied. See the Licence for the specific language
%   governing permissions and limitations under the Licence.
  
if(nargin < 4) 
  datadir1 = 'be.out';
  if(nargin < 3)
    rootdir = '.';
    if(nargin < 2)
      datatype = 4;   
      if(nargin < 1)
        name = 'output';
      end
    end
  end
end
%%
% information about the geometry, these are the default blender export locations
%%
if(nargin < 6) 
  evaldir = sprintf('%s/../../EvaluationGrids/3_ARI', ...
                    rootdir);
  if(nargin < 5)
    meshdir = sprintf('%s/../../ObjectMeshes/Reference', ...
                      rootdir);
  end
end

isoctave = is_octave();
nojumps = 1;
elementwise = 1;

if ( (datatype < 2) )
  evaldir = [];
end
%% check consistency between datatype and evalname, no need to
% try to find values for the evaluation grid if no evaluation grid is
% given
if( isempty(evaldir) && (datatype > 1) )
    datatype = 1;
end

if( isempty(meshdir) && (datatype < 2 || datatype > 3) )
    datatype = 3;
end

%% read the grid data of the scatterer, only one scattererfile is assumed
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
  if( isoctave )
    %% importdata works like a charm here for octave,
    % however for MATLAB you're out of luck and need the slow path
    Elems = importdata(filename);
    Elems(1,:) = []; % first line contains the number of elements
    if( columns(Elems) == 7 )
      Elem3 = Elems(:,1:4);
      Elem4 = [];
    else
      e3 = find (isnan (Elems(:,8) ) );
      if(~isempty(e3))
	Elem3 = Elems(e3,1:4);
	Elem(e3,:) = [];
      end
      Elem4 = Elem(:,1:5);
    end
  else
    [Elem3,Elem4] = read_elements(fh); % triagles and quadrilaterals
  end

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

  if(isoctave)
    EElems = importdata(filename);
    EElems(1,:) = [];
    if(columns(EElems) == 7)
      EElem3 = EElems(:,1:4);
      EElem4 = [];
    else
      e3 = find (isnan (EElems(:,8) ) );
      if(~isempty(e3))
	EElem3 = EElems(e3,1:4);
	EElems(e3,:) = [];
      else
	EElem3 = [];
      end
      EElem4 = EElems(:,1:5);
    end
  else
    [EElem3,EElem4] = read_elements(fh);
  end
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
% in different parts) the script needs to be extended, or
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

if(nnodes > 0)   % there are nodes for the scatterer
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

%%do the same for the evalgrid, now the nodes definitely do not
%%start with 0

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
    %% now that is a bit spooky because of the string in file
    % but it seems to work, however just on octave
    %pB = dlmread(filename);
    %pB(1:3,:) = [];
    printpb = 1;
    %% should work in both matlab and octave
    % throw away the first 3 lines
    
    D = importdata(filename, ' ', 3);
    pB = D.data;
    pB(:,2) = pB(:,2) + pB(:,3)*1.0i;
    pB(:,3) = [];
    
    %% please NOTE: It is assumed that the Elements are ordered
    % and have no index jumps
  end
  if( ismember(datatype,[2,4,6]) )
    %% read the pressure data at the eval Grid
    filename = sprintf('%s/%s/be.%d/pEvalGrid',rootdir, ...
                       datadir1, i);
    %%        pE = dlmread(filename2);
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
    %%vB = dlmread(filename3);
    D = importdata(filename, ' ', 3);
    vB = D.data;
    %%vB(1:3,:) = [];
    vB(:,2) = vB(:,2) + vB(:,3)*1.0i;
    vB(:,3) = [];
    printvb = 1;
	
  end
  if( ismember(datatype, [3,5,6]) )
    %% read the velocity at the eval grid
    filename = sprintf('%s/%s/be.%d/vEvalGrid',rootdir, ...
                       datadir1, i);
    D = importdata(filename, ' ', 3);
    vE = D.data;
  
    %% this should ba a matix n times 7
    vE(:,2) = vE(:,2) + 1.0i*vE(:,3);
    vE(:,3) = vE(:,4) + 1.0i*vE(:,5);
    vE(:,4) = vE(:,6) + 1.0i*vE(:,7);
    vE(:,7) = [];
    vE(:,6) = [];
    vE(:,5) = [];
    printve = 1;
  end
    
% the value at the boundary is given at the collocations nodes which is the
% midpoint of each element, however in order to get a smoother
% graph we decided to assign a value to each nodes of the grid
% if you do not like the idea, just reformulate the vtk file
% with cell data instead of point data, this should work too
  if(printpb && ~elementwise) 
    pBound = zeros(nnodes,1);
  end
  if(printvb && ~elementwise) 
    vBound = zeros(nnodes,1);
  end
  
  if( ~elementwise ) % get the values at the BEM nodes by interpolation
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
          %% take the mean value of all elements ==
          % collocation nodes containing the vertex
          % 
          % the index juggling may seem a bit unnecessary
          % right now, but remember, it is possible to have
          % jumps in the node numbers and element numbers so
          % this can be helpful if the script will be
         % extended for that case, howevert that takes a long time for
	 % large meshes, thus assume nojumps and do it faster
	  if(nojumps)
	    pBound(j) = mean( pB(n31,2) );
	  else
            for j1 = 1:length(n31),
              pBound(j) = pBound(j) + pB( find( pB(:,1) == Elem3(n31(j1),1) ...
                                              ), 2);
            end
	  end
	end
	%% particle velocity on the boundary
	if(printvb)
	  if( nojumps )
	    vBound(j) = mean( vb(n31,2) );
	  else
            for j1 = 1:length(n31),
              vBound(j) = vBound(j) + vB( find( vB(:,1) == Elem3(n31(j1),1) ...
                                              ), 2);
            end
          end
	end
      end
      %% quadrilaterals
      if( n41 )
        if(printpb)
	  if( nojumps )
	    pBound(j) = mean( pB(n41,2) );
	  else
            for j1 = 1: length(n41),
              pBound(j) = pBound(j) + pB( find( pB(:,1) == Elem4(n41(j1),1) ...
                                              ), 2);
            end
          end
	end
        if(printvb)
	  if( nojumps )
	    vBound(j) = mean ( vB(n41,2) );
	  else
            for j1 = 1:length(n41),
              vBound(j) = vBound(j) + vB( find( vB(:,1) == Elem4(n41(j1),1) ...
                                              ), 2);
            end
          end 
        end
      end
      if(printpb) % use the average over all elements containing
		  % this node
        pBound(j) = pBound(j) / n;
      end
      if(printvb)
        vBound(j) = vBound(j) / n;
      end
    end % loop nodes
  end % if elementwise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          write the vtk-files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  filename = sprintf('%s_bound_%d.vtk', name, i);


  %% in case of elementwise == 1
  % pB, vB, pE, and vE contain number in column 1
  % and the respective complex value in the rest of the columns
  if(printpb + printvb == 2)
    if( elementwise )
      writefile(filename, 0, nodes, Elem3, Elem4, pB, vB);
    else
      writefile(filename, 0, nodes, Elem3, Elem4, pBound, vBound);
    end
  else
    if(printpb)
      if( elementwise )
	writefile(filename, 0, nodes, Elem3, Elem4, pB);
      else
	writefile(filename, 0, nodes, Elem3, Elem4, pBound);
      end
    end
    if(printvb)
      if( elementwise )
	writefile(filename, 1, nodes, Elem3, Elem4, vB);
      else
	writefile(filename, 1, nodes, Elem3, Elem4, vBound);
      end
    end
  end

  
 %   if(~isempty(pE))
 %     pE(:,1) = [];
 %   end
 %   if(~isempty(vE))
 %     vE(:,1) = [];
 %   end
  
  filename = sprintf('%s_eval_%d.vtk', name, i);

  if(printpe + printve == 2)
    writefile(filename, 2, enodes, EElem3, EElem4, pE, vE);
  else
    if(printpe)
      writefile(filename, 2, enodes, EElem3, EElem4, pE);
    end
    if(printve)
      writefile(filename, 3, enodes, EElem3, EElem4, vE);
    end
  end
end  % loop over freqsteps



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                 Function writefile
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writefile(filename, datatype, N, E3, E4, data1, data2)    
%
% write the mesh data and the values at the mesh nodes to a vtk file
% 
% input: filename  : name of the output file
%        datatype  : what should be shown
%                    0 .... pressure (and velocity) boundary
%                    1 .... velocity boundary
%                    2 .... pressure (and velocity) eval
%                    3 .... velocity eval
%        N         : list of nodes (matrix nx3)
%        E3        : list of all triangular elements (matrix mx3)
%        E4        : list of all quadrilateral elements (matrix mx4)
%        data1     : vector containing the comlex pressure or the
%                    complex velocity
%        data2     : optional data vector if both pressure and
%                    velocity shall be displayed, can only be the
%                    velocity
%
% One of the problems is that the velocity on the evaluation grid has
% three components, we have to take care of that
%
% also there are two different kind of data1 types
% 1) data1 contains the data on the nodes of BE and Evalelem
%    in this case the first columns of data1 (and data2) already
%    contains the value
% 2) data1 contains the data on the collocation nodes, i.e. the
%    midpoints of the elements: in this case data1 (and data2)
%    contains the indexnumber of the element in the first columns
% In the first case no ordering is necessary, that has already been
% done when taking the mean value. In the second case we have to split
% between triangular and quadrilateral elements
if(nargin < 7)
    data2 = [];
end
% E3 and E4 still contain the indices of the elements
% get rid of them



if(~isempty(E3))
  i3 = E3(:,1);
  E3(:,1) = [];
else
  i3 = [];
end
if(~isempty(E4))
  i4 = E4(:,1);
  E4(:,1) = [];
else
  i4 = [];
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
%%%%%%%
%
%  geometry is written, lets think about the data we have
%
%%%%%%%%

if( size(N,1) == size(data1,1) )
  %
  % write point data, this could either be the mean value over the
  % boundary elements or data from the evaluationgrid
  %
  fprintf(fh, 'POINT_DATA %d\n', size(N,1) );
  if (datatype == 0) 
    fprintf(fh, 'SCALARS Abs(BEM_Pressure)_(dB) double\n');
    data = 20*log10(abs(data1)/2e-5);
    theta = angle(data1);
    fprintf(fh, 'LOOKUP_TABLE default\n');
    fprintf(fh, '%e\n', data);
				%
    fprintf(fh, 'SCALARS Phase(BEM_Pressure) double \n');
    fprintf(fh, 'LOOKUP_TABLE default\n');
    fprintf(fh, '%e\n', theta);
  end

  if(datatype == 1)
    fprintf(fh, 'SCALARS Abs(BEM_Velocity) double\n');
    data = abs(data1);
    theta = angle(data1);
    fprintf(fh, 'LOOKUP_TABLE default\n');
    fprintf(fh, '%e\n', data);
    fprintf(fh, 'SCALARS Phase(BEM_Velocity) double\n');
    fprintf(fh, 'LOOKUP_TABLE default\n');
    fprintf(fh, '%e\n', theta);
  end

  if(datatype == 2)
    fprintf(fh, 'SCALARS Abs(Eval_Pressure)_(dB) double\n');
    data = 20*log10(abs(data1(:,2)/2e-5));
    theta = angle(data1(:,2));
    fprintf(fh, 'LOOKUP_TABLE default\n');
    fprintf(fh, '%e\n', data);
    fprintf(fh, 'SCALARS Phase(Eval_Pressure) double\n');
    fprintf(fh, 'LOOKUP_TABLE default\n');
    fprintf(fh, '%e\n', theta);
  end

  if(datatype == 3)
    data(:,1:3) = abs(data1(:,2:4));
    theta(:,1:3) = angle(data1(:,2:4));
    fprintf(fh, 'SCALARS Abs(Eval_V_x) double\n');
    fprintf(fh, 'LOOKUP_TABLE default\n');
    fprintf(fh, '%e\n', data(:,1));
    fprintf(fh, 'SCALARS Phase(Eval_V_x) double\n');
    fprintf(fh, 'LOOKUP_TABLE default\n');
    fprintf(fh, '%e\n', theta(:,1));
				%
    fprintf(fh, 'SCALARS Abs(Eval_V_y) double\n');
    fprintf(fh, 'LOOKUP_TABLE default\n');
    fprintf(fh, '%e\n', data(:,2));
    fprintf(fh, 'SCALARS Phase(Eval_V_y) double\n');
    fprintf(fh, 'LOOKUP_TABLE default\n');
    fprintf(fh, '%e\n', theta(:,2));
				%
    fprintf(fh, 'SCALARS Abs(Eval_V_z) double\n');
    fprintf(fh, 'LOOKUP_TABLE default\n');
    fprintf(fh, '%e\n', data(:,3));
    fprintf(fh, 'SCALARS Phase(Eval_V_z) double\n');
    fprintf(fh, 'LOOKUP_TABLE default\n');
    fprintf(fh, '%e\n', theta(:,3));
  end
  %

  %

else
  %% so we have data for the elements, split between triangles + quadril
  %
  % do the triangles first
  i = ismember( data1(:,1) , i3);
  j3 = find( i == 1 );
  j4 = find( i == 0 );
  % there may be a slight problem when we look at the particle
  % velocity at the evalgrid, because there are three components

  c = data1(j3,2:end);
  c = [c;data1(j4,2:end)];

  fprintf(fh, 'Cell_DATA %d\n', length(c));
  [~,i] = size(c);
%  if( i > 1 )   % we have a velocity
%    theta = angle(c(:,1) + 1.0i*c(:,2));
%    keyboard
%  else
  theta = angle(c);
%  end
%    if ( (datatype == 2) )
%      %% velocity data
%      if( i > 1 )
%        fprintf(fh, 'SCALARS Eval_Abs(v_x) double 1\n');
%        fprintf(fh, 'LOOKUP_TABLE default\n');
%        fprintf(fh, '%e\n', abs(c(:,1)));
%   				%
%        fprintf(fh, 'SCALARS Eval_Abs(v_y) double 1\n');
%        fprintf(fh, 'LOOKUP_TABLE default\n');
%        fprintf(fh, '%e \n', abs(c(:,2)));
%   				%	      
%        fprintf(fh, 'SCALARS Eval_Abs(v_z) double 1\n');
%        fprintf(fh, 'LOOKUP_TABLE default\n');
%        fprintf(fh, '%e \n', abs(c(:,3)));
%      else
%        fprintf(fh, 'SCALARS Eval_Abs(v) double 1\n');
%        fprintf(fh, 'LOOKUP_TABLE default\n');
%        fprintf(fh, '%e\n',abs(c));
%      end
%    end
  %% this needs to be the sound pressure
  if(datatype == 0)
    fprintf(fh, 'SCALARS BEM_Abs(Pressure)_(db) double 1\n');
    fprintf(fh, 'LOOKUP_TABLE default\n');
    fprintf(fh, '%e\n', 20 * log10(abs(c)/2.0e-5));

    fprintf(fh, 'SCALARS BEM_Phase(p) double 1\n');
    fprintf(fh, 'LOOKUP_TABLE default\n');
    fprintf(fh, '%e\n', theta);
    
  end

  if(datatype == 1)
    fprintf(fh, 'SCALARS Abs(V_BEM) double 1\n');
    fprintf(fh, 'LOOKUP_TABLE default\n');
    fprintf(fh, '%e\n', abs(c));

    fprintf(fh, 'SCALARS BEM_Phase(v) double 1\n');
    fprintf(fh, 'LOOKUP_TABLE default\n');
    fprintf(fh, '%e\n', theta);
  end
				% now the quadrilaterals
%  i = ismember( data1(:,1) , i4);
%  j = find( i == 1 );
%  c = data1(j,2);
%  theta = [theta;angle(c)];
%  fprintf(fh, '%e\n', 20 * log10(abs(c)/2.0e-5));
	%  if( i > 1 )
	%    fprintf(fh, 'SCALARS Phase_x double 1\n');
	%    fprintf(fh, 'LOOKUP_TABLE default\n');
	%    fprintf(fh, '%e\n', theta(:,1));
	% 				%
	%    fprintf(fh, 'SCALARS Phase_y double 1\n');
	%    fprintf(fh, 'LOOKUP_TABLE default\n');
	%    fprintf(fh, '%e\n', theta(:,2));
	% 				%
	%    fprintf(fh, 'SCALARS Phase_z double 1\n');
	%    fprintf(fh, 'LOOKUP_TABLE default\n');
	%    fprintf(fh, '%e\n', theta(:,3));
	%  else


 end

if(~isempty(data2))
  % this can only be a velocity
  if(datatype == 0)   %% BEM 
    fprintf(fh, 'SCALARS BEM_Abs(Velo) double\n');
    fprintf(fh, 'LOOKUP_TABLE default\n');
    i = find(abs(data2(:,2)) < 1e-10);
    data2(i,2) = 1e-10;
    data = abs(data2(:,2));
    theta = angle(data2(:,2));
    fprintf(fh, '%e\n',data);
%
    fprintf(fh, 'SCALARS BEM_Phase(Velo) double\n');
    fprintf(fh, 'LOOKUP_TABLE default\n');
    fprintf(fh, '%e\n',theta);
  end
%
  if(datatype == 2) % eval
    fprintf(fh, 'SCALARS Eval_Abs(v_x) double\n');
    fprintf(fh, 'LOOKUP_TABLE default\n');
    c = abs(data2(:,2));
    theta_x = angle(data2(:,2));
    fprintf(fh, '%e\n', c);
				%
    fprintf(fh, 'SCALARS Eval_Abs(v_y) double\n');
    fprintf(fh, 'LOOKUP_TABLE default\n');
    c = abs(data2(:,3));
    theta_y = angle(data2(:,3));
    fprintf(fh, '%e\n', c);
				%     
    fprintf(fh, 'SCALARS Eval_Abs(v_z) double\n');
    fprintf(fh, 'LOOKUP_TABLE default\n');
    c = abs(data2(:,4));
    theta_z = angle(data2(:,4));
    fprintf(fh, '%e\n', c);
%%%%%%%%%%%%%%%%%%%%      
    fprintf(fh, 'SCALARS Eval_Phase_x double\n');
    fprintf(fh, 'LOOKUP_TABLE default\n');
    fprintf(fh, '%e\n', theta_x);
    %%
    fprintf(fh, 'SCALARS Eval_Phase_y double\n');
    fprintf(fh, 'LOOKUP_TABLE default\n');
    fprintf(fh, '%e\n', theta_y);
				%
    fprintf(fh, 'SCALARS Eval_Phase_z double\n');
    fprintf(fh, 'LOOKUP_TABLE default\n');
    fprintf(fh, '%e\n', theta_z);
  end    
end
fclose(fh);
% function end





function [Elem3,Elem4] = read_elements(fh)
% This is slow and just for MATLAB users
% reads element number and element nodes from filehandle
% assumes mesh2hrtf format
% there should be a faster version for MATLAB but importdata fills the
% matrix incorrectly, and dlmread fills it with 0 which is not good
% too
  
 
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

function r = is_octave ()
  persistent x;
  if (isempty (x))
    x = exist ('OCTAVE_VERSION', 'builtin');
  end
  r = x;

