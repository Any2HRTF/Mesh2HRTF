function ncdata2vtkparse(name,datatype,fstep,rootdir,datadir, filename)
% exports data from NumCalc results to vtk files for display with paraview
% 
%% input: 
%   name: string: Name of the vtk file(s) to be produced
%            the outputfiles will be in ascii format and named
%                      <name>_boundary_<frequencystep>.vtk
%                      <name>_eval_<frequencystep>.vtk 
%   datatype: integer: type of data to be presented
%       0    pressure boundary
%       1    part. velocity boundary
%       2    pressure evalgrid
%       3    part. velocity evalgrid (not implemented yet)
%       4    pressure everywhere (default)
%       5    part. velocity everywhere (not implemented yet)
%       6    display all (not implemented yet)
%   fstep:   integer vector: numbers of the frequencysteps of interest,
%            default all steps in datadir: fstep = 0
%   rootdir: string. name of the root directory containing the be.out file,
%        default the current directory
%   datadir: string. name of the directory containing the data: default
%            be.out
%   inpfile: strin. name of the input file, containing the mesh data: default
%            NC.inp
% the data will be stored in the vtk file as Amplitude and Phase values
%
% written by kreiza, if you have any questions, suggestions 
%   send me an email, depending on my workload I will try to reply
%   in time: wolfgang.kreuzer@oeaw.ac.at
%
%
%  The script should work with Octave as well as Matlab, however
%  Matlab is a pain in the ass, because
%    a) drmwrite is not encouraged by mathworks
%    b) their alternative writematrix does not know how to append to a
%    file in the 2019 version of matlab
%    c) fprintf was also relatively slow
%
%   FYI: We stil use the drmwrite option
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

  if(nargin < 6)
    filename = "NC.inp";
    if(nargin < 5) 
      datadir = 'be.out';
      if(nargin < 4)
        rootdir = '.';
	if(nargin < 3)
	  fstep = 0;
          if(nargin < 2)
            datatype = 4;   
            if(nargin < 1)
              name = 'output';
            end
          end
	end
      end
    end
  end


  
  if( ismember(datatype,[5,6]) )
    error("Sorry not implemented yet");
  end
  ifp = 0;
  ifp = fopen(filename,'r');
  if(ifp < 0)
    error(sprintf("Sorry could not open %s\n",filename));
  end

  nodes = [];

  inpline = [];
  while ~strncmp(inpline,"NODES",5);
    inpline = fgets(ifp);
  end 
  %% the next line should be either a comment or a filename
  filename = strtrim(fgets(ifp)); % to get rid of leading and
                                  % trailing whitespaces
  while ~strncmp(filename,"ELEMENTS",4)
    if(filename(1)~="#")
      %% not comment thus, a real line
      %% since the nodes have all the same columns,
      %% importdata also works with matlab
      %% the 1 gets rid of the first entry in the file, which is the
      %% number of nodes in the file

      % lets assume that there is a trailing whitespace, and get rid of that
      % I haven't found a way to check if there is one
      N = importdata(filename,' ',1);
      N = N.data;
      nodes = [nodes;N(1:end,:)];
    end
    filename = strtrim(fgets(ifp));
  end			       % there may be more then one node files

  [nnodes,~] = size(nodes);
  
  E = [];
  Elem3 = [];
  Elem4 = [];
  %% do the same thing as above
  filename = strtrim(fgets(ifp));
  while(~strncmp(filename,"BOUNDARY",5))
    if(filename(1)~="#")
      if is_octave()
	Elems = importdata(filename,' ',1);
	Elems = Elems.data;
      else
	Elems = readmatrix(filename);
      end
    end
    
    [~,j] = size(Elems); 
    if( j == 7 )
      %% we have only triangles, hurray
      Elem3 = [Elem3;Elems];
    else
      %% so we have a combination of quadrilaterals and triangles, or
      %% just quadrilaterals
      i = find(isnan(Elems(:,end)));
      if(isempty(i))
	%% only quadrilaterals
	Elem4 = [Elem4;Elems];
      else
	Elem3 = [Elem3;Elems(i,1:7)];
	Elems(i,:) = [];
	Elem4 = [Elem4;Elems];
      end
    end
    filename = strtrim(fgets(ifp));
  end

  disp("Read Geometry");

%% find out if there are evalelements
  if( isempty( Elem3 ) )
    EElem3 = [];
  else
    i = find( Elem3(:,5) == 2 );
    if(~isempty(i))
      EElem3 = Elem3(i,:);
      Elem3(i,:) = [];
    else
      EElem3 = [];
    end
  end

  if( isempty( Elem4 ) )
  EElem4 = [];
  else
    i = find( Elem4(:,6) == 2 );
    if(~isempty(i))
      EElem4 = Elem4(i,:);
      Elem4(i,:) = [];
    else
      EElem4 = [];
    end
  end

  
  
% check consistency between datatype and evalname
%if( isempty(evaldir) && (datatype > 1) )
%    datatype = 1;
%end

%% clean up the elements and the mesh
%% problem is, that NumCalc allows jumps in the element numbers
%% vtk does not


  


  if( isempty(EElem3) && isempty(EElem4) && (datatype > 1) )
    datatype = 1;
  end
  
  if( isempty(Elem3) && isempty(Elem4) && (datatype < 2 || datatype > 3) )
    datatype = 3;
  end

  %% lets see what data we have in the first place
  nentries = length( dir( sprintf('%s/%s',rootdir,datadir) )); 
  if(  nentries > 3)
    %% this works if we have more then one entry
    [fnames] = dir( sprintf('%s/%s/be.*',rootdir,datadir) );
    nf = length( fnames );
    fstep0 = sscanf([fnames.name],"be.%d");
  else
    if( nentries > 2 )
      [fnames] = dir( sprintf('%s/%s',rootdir,datadir) );
      fstep0 = sscanf(fnames(3).name,"be.%d");
    else
      fstep0 = [];
    end
  end

  
    
  if(fstep == 0)
    %% works under linux,
    %% windows user may have to adapt this part
    fstep = fstep0;
  else
    fstep = intersect(fstep,fstep0);
  end

  
  
  printpb = 0;
  printpe = 0;
  printvb = 0;
  printve = 0;
  pB = [];
  pE = [];
  vB = [];
  vE = [];
  %%
  %% checkout which nodes are eval nodes and which aren't
  %%
  
  nnr = nodes(:,1);
  evalnodes=[];

  if(~isempty(EElem3))
    i = ismember(nnr,EElem3(:,2:4));
  else
    i = [];
  end
  if(~isempty(i))
    evalnodes = nodes(i,:);
    nodes(i,:) = [];
  end
  if( ~isempty(EElem4))
    i = ismember(nnr,EElem4(:,2:5));
  else
    i = [];
  end
  if(~isempty(i))
    evalnodes = [evalnodes;nodes(i,:)];
    nodes(i,:) = [];
  end
  % nodes should contain only bemnodes now

  Elem = [];
  EElem = [];
  % we will use export_vtk thus triangles are assumed to have -1 in
  % the last column
  if(~isempty(Elem3))
    Elem = Elem3(:,1:4);
  end
  if(~isempty(Elem4))
    if(~isempty(Elem))
      Elem(:,5) = -1;
    end
    Elem = [Elem;Elem4(:,1:5)];
  end
  Elem(:,1) = []; % get rid of the elemnumber

  if(~isempty(EElem3))
    EElem = EElem3(:,1:4);
  end
  if(~isempty(EElem4))
    if(~isempty(EElem))
      EElem(:,5) = -1;
    end
    EElem = [EElem;EElem4(:,1:5)];
  end
  EElem(:,1) = [];


  %% in theory there could be jumps in the node and element numbers
  %% for the results of numcalc this is not a problem, however, for
  %% the vtk geometry this could become a problem
  nn = find(diff(nodes(:,1)) - 1);

  if ~isempty(nn)
    for n = nn,
      delta = nodes(n+1,1) - n+1;
      [i,j] = find(Elem3 > nodes(n,1)); % n is coorect
      Elem3(i,j) = Elem3(i,j) - delta;
      [i,j] = find(Elem4 > nodes(n,1));
      Elem4(i,j) = Elem4(i,j) - delta;
    end
  end
  %% the same for the eval
  nn = find(diff(evalnodes(:,1)) - 1);
  if ~isempty(nn)
    for n = nn,
      delta = evalnodes(n+1,1) - n+1;
      [i,j] = find(EElem3 > evalnodes(n,1)); % n is coorect
      EElem3(i,j) = EElem3(i,j) - delta;
      [i,j] = find(EElem4 > evalnodes(n,1));
      EElem4(i,j) = EElem4(i,j) - delta;
    end
  end

  % we do not need the nodenumbers anymore now
  evalnodes(:,1) = [];
  nodes(:,1) = [];


  
%
% we look at the data for all frequency steps so one can make an
% animation out of the data, if you want things just for one
% specific frequency you'll need to change nf below
%
  disp("Starting Frequency Loop");
  for istep0 = 1:length(fstep),
    istep = fstep(istep0);
    if( ismember(datatype,[0,4,6]) )  
      %% read the pressure data at the boundary
      filename = sprintf('%s/%s/be.%d/pBoundary',rootdir, ...
                         datadir, istep);
      printpb = 1;
      %% should work in both matlab and octave
      % throw away the first 3 lines
      % the first column should be the element number
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
                         datadir, istep);
      
      D = importdata(filename, ' ', 3);
      pE = D.data;
      pE(:,2) = pE(:,2) + pE(:,3)*1.0i;
      pE(:,3) = [];
      printpe = 1;
    end
    if( ismember(datatype, [1,5,6]) )
      %% read the velocity at the boundary
      filename = sprintf('%s/%s/be.%d/vBoundary',rootdir, ...
                         datadir, istep);
      D = importdata(filename, ' ', 3);
      vB = D.data;
      vB(:,2) = vB(:,2) + vB(:,3)*1.0i;
      printvb = 1;

    end
    if( ismember(datatype, [5,6]) )
        %% read the velocity at the eval grid
        filename = sprintf('%s/%s/be.%d/vEvalGrid',rootdir, ...
                            datadir, istep);
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
    %  if(printpb)
    %      pBound = zeros(nnodes,1);
    %  end
    %  if(printvb) 
    %      vBound = zeros(nnodes,1);
    %  end
    %   
    %nnodes = rows(nodes);
    %% this is some old code snipped, that takes some time, however
    %% if you want to have a smooth graph, that works with nodal
    %% values instead of constant elements, the lines should work
    %% however, currently we use CELL structures to display the
    %% constant elements
   %    disp("Get the values for the BEM nodes")
   %    % for paraview the nodes need to have the value not the elements
   %    for j = 1:nnodes,
   %      if( ~isempty(Elem3) )
   %            % find the elements that contain the nodenumber
   %        [n31,n32] = find( nodes(j,1)  == Elem3(:,2:4) );
   %     	    % n31 contains the rownumbers of Elem3 that contain node j
   %      else
   %        n31 = 0;
   %      end
   %      
   %      if( ~isempty(Elem4) )
   %        [n41,n42] = find( nodes(j,1) == Elem4(:,2:5) );
   %      else
   %        n41 = 0;
   %      end
   %      n = length(n41) + length(n31);
   %      
   %      if( n31 )
   %        if(printpb)
   %                    % take the mean value of all elements ==
   %                    % collocation nodes containing the vertex
   %                    % 
   %                    % the index juggling may seem a bit unnecessary
   %                    % right now, but remember, it is possible to have
   %                    % jumps in the node numbers and element numbers so
   %                    % this can be helpful if the script will be
   %                    % extended for that case
   %          for j1 = 1:length(n31),
   %            pBound(j) = pBound(j) + pB( find( pB(:,1) == Elem3(n31(j1),1) ...
   %                                            ), 2);
   %          end
   %        end
   %        if(printvb)
   %          for j1 = 1:length(n31),
   %            vBound(j) = vBound(j) + vB( find( vB(:,1) == Elem3(n31(j1),1) ...
   %                                            ), 2);
   %          end
   %        end 
   %      end
   %      if( n41 )
   %        if(printpb)
   %          for j1 = 1: length(n41),
   %            pBound(j) = pBound(j) + pB( find( pB(:,1) == Elem4(n41(j1),1) ...
   %                                            ), 2);
   %          end
   %        end
   %        if(printvb)
   %          for j1 = 1:length(n41),
   %            vBound(j) = vBound(j) + vB( find( vB(:,1) == Elem4(n41(j1),1) ...
   %                                            ), 2);
   %          end
   %        end 
   %        
   %      end
   %      if(printpb)
   %        pBound(j) = pBound(j) / n;
   %      end
   %      if(printvb)
   %        vBound(j) = vBound(j) / n;
   %      end
   %    end
    
        
%%% put everything together in one file for now, think about the
%%% rest later

   %   if(!isempty(pE))
   %     pBound(pE(:,1)) = pE(:,2);
   %   end

    disp("Data is read")
    %% just a temporary solution, sooner or later we may want
    %% to distinguish between eval and bem
    %Elem3 = [Elem3;EElem3];
    %Elem4 = [Elem4;EElem4];
    
    %% velocity won't work this way, because it may be a vector
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %          write the vtk-files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %% if you are lazy like me, there is a chance that nodes also
    %% contains evaluation nodes, just saying
    filename = sprintf('%s_bound_%d.vtk', name, istep);
    if(printpb)
      pB = pB(:,2);
      i = find( pB == 0 );
      pB = 20 * log10(abs(pB)/2.0e-5);
      pB(i) = -500;
      export_vtk(nodes,Elem,filename,pB);
    end
    if(printvb)
      vB = vB(:,2);
      i = find( vB == 0);
      vB = 20 * log10(abs(vB)/5.0e-8);
      vB(i) = -500;
      export_vtk(nodes,Elem,filename,vB);
    end
    filename = sprintf('%s_eval_%d.vtk', name, istep);
    

    if(printpe)
      pE = pE(:,2);
      i = find( pE == 0 );
      pE = 20 * log10(abs(pE)/2.0e-5);
      pE(i) = -500;
      export_vtk(evalnodes,EElem,filename,pE);
    end
%      if(printpb + printvb == 2)
%          writefile(filename, nodes, Elem3, Elem4, pBound, vBound);
%      else
%          if(printpb)
%              writefile(filename, nodes, Elem3, Elem4, pBound);
%          end
%          if(printvb)
%              writefile(filename, nodes, Elem3, Elem4, vBound);
%          end
%      end
%      
%      if(~isempty(pE))
%          pE(:,1) = [];
%      end
%      if(~isempty(vE))
%          vE(:,1) = [];
%      end
%      filename = sprintf('%s_eval_%d.vtk', name, i);
    
    % now there is a slight problem, because we have not defined the
    % enodes yet, this is really lazy, we do not distinguish between
    % eval nodes and bem nodes
 %    enodes = nodes;
 %    if(printpe + printve == 2)
 %      writefile(filename, enodes, EElem3, EElem4, pE, vE);
 %    else
 %      if(printpe)
 %        writefile(filename, enodes, EElem3, EElem4, pE);
 %      end
 %      if(printve)
 %        writefile(filename, enodes, EElem3, EElem4, vE);
 %      end
 %    end
    disp("Data is written  ");
  end  % loop over freqsteps
end


%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    %
%    %                 Function writefile
%    %
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    function writefile(filename, N, E3, E4, data1, data2)    
%    %
%    % write the mesh data and the values at the mesh nodes to a vtk file
%    % 
%    % input: filename  : name of the output file
%    %        N         : list of nodes (matrix nx3)
%    %        E3        : list of all triangular elements (matrix mx3)
%    %        E4        : list of all quadrilateral elements (matrix mx4)
%    %        data1     : vector containing the comlex pressure or the
%    %                    complex velocity
%    %        data2     : optional data vector if both pressure and
%    %                    velocity shall be displayed
%    if(nargin < 6)
%        data2 = [];
%    end
%    % E3 and E4 still contain the indices of the elements
%    % get rid of them
%     
%    if(~isempty(E3))
%        E3(:,1) = [];
%    end
%    if(~isempty(E4))
%        E4(:,1) = [];
%    end
%     
%    nnodes = size(N,1);
%     
%     
%     
%    %% lets do some correction of element numbers and vtk format
%    %% N contains nodenumber and coordinates
%    nrmin = min( N(:,1) );
%    if(!isempty(E3))
%      E3 = E3(:,1:3);
%    end
%    if(!isempty(E4))
%      E4 = E4(:,1:4);
%    end
%     
%     
%    %% clean up E3 and E4
%    %% this only works if nodes does not contain more nodes then the
%    %% elements are using,
%    %% however there still could be jumps, clean that,
%    %% if you want to use that feature uncomment the next lines
%    %   for i = 1:nnodes,
%    %     [i1,i2] = find( E3 == N(i,1) );
%    %     for j = 1:length(i1),
%    %       E3(i1(j),i2(j)) = i - 1;
%    %     end
%    %     [i1,i2] = find( E4 == N(i,1) );
%    %     for j = 1:length(i1),
%    %       E4(i1(j),i2(j)) = i - 1;
%    %     end
%    %   end
%     
%    % lets assume that there are not jumps in the nodes, however nodes may
%    % contain more nodenumbers than necessary for the elems
%     
%     
%    dlist = [];
%    for i = 1:nnodes,
%      j = [];
%      if( !isempty(E4) )
%        j = find( N(i,1) == E4);
%      end
%     
%      if( !isempty(E3) )
%        j = find( N(i,1) == E3);
%      end
%     
%      if( isempty(j) )
%        dlist = [dlist,i];
%      end
%    end
%     
%    % get rid of the nodes that are not used
%    N(dlist,:) = [];
%    nmin = min(N(:,1));
%     
%    nnodes = rows(N);
%    % as I said assume, that there are no jumps
%    if(!isempty(E3))
%      E3 = E3 - nmin;
%    end
%     
%    if(!isempty(E4))
%      E4 = E4 - nmin;
%    end
%     	   
%     
%    N = N(:,2:4);
%     
%     
%     
%    fh = fopen(filename,'w');
%    if (fh == -1) 
%        error('Sorry cound not open the output file.')
%    end
%    %
%    %  write the header
%    %
%    fprintf(fh, '# vtk DataFile Version 3.0\n');
%    fprintf(fh, 'Generated by data2vtk\n');
%    fprintf(fh, 'ASCII\n');
%    %    fprintf(fh, '\n');
%    fprintf(fh, 'DATASET POLYDATA\n');
%    fprintf(fh, 'POINTS %d double\n',nnodes);
%    fprintf(fh, '%e %e %e\n', N');
%    %
%    %  write the triangle geo data
%    %
%     
%    fprintf(fh, 'POLYGONS %d %d\n',size(E3,1)+size(E4,1), ...
%            size(E3,1)*4+size(E4,1)*5);
%     
%    if(size(E3,1) > 0)
%        fprintf(fh, '3 %d %d %d\n', E3');
%    end
%    %
%    % write 4-sided geometry
%    %
%    if(size(E4,1) > 0)
%        fprintf(fh, '4 %d %d %d %d\n', E4');
%    end
%     
%    fprintf(fh, 'POINT_DATA %d\n', size(N,1) );
%    fprintf(fh, 'SCALARS Amplitude1(dB) double\n');
%    fprintf(fh, 'LOOKUP_TABLE default\n');
%     
%    data = 20*log10(abs(data1)/2e-5);
%    % this may happen if you use weird inputfiles
%    i = find(data == -Inf);
%    data(i) = -1000;
%    fprintf(fh, '%e\n', data);
%     
%    fprintf(fh, 'SCALARS Phase1(rad) double\n');
%    fprintf(fh, 'LOOKUP_TABLE default\n');
%    data = angle(data1);
%    fprintf(fh, '%e\n', data);
%     
%    if(~isempty(data2))
%        fprintf(fh, 'SCALARS Amplitude2(dB) double\n');
%        fprintf(fh, 'LOOKUP_TABLE default\n');
%        i = find(abs(data2) < 1e-10);
%        data2(i) = 1e-10;
%        data = 20*log10(abs(data2)/2e-5);
%        fprintf(fh, '%e\n', data);
%     
%        fprintf(fh, 'SCALARS Phase2(rad) double\n');
%        fprintf(fh, 'LOOKUP_TABLE default\n');
%        data = angle(data2);
%        fprintf(fh, '%e\n', data);
%    end
%    fclose(fh);
%     				% function end
%     
%     



%   function [Elem3,Elem4] = read_elements(fh)
%   % reads element number and element nodes from filehandle
%   % assumes mesh2hrtf format
%   % would be easier to use dlmread(filename,"emptyvalue",-1)
%   nelems = str2num(fgets(fh));
%   Elem3 = [];
%   Elem4 = [];
%   for j = 1:nelems
%     str = fgets(fh);
%     el = sscanf(str,'%f');
%     n = length(el);
%     if(n == 7)
%    				% we have a triangle
%       el = el(1:4);
%       Elem3 = [Elem3;el'];  
%     else
%    				% quadrilateral
%       el = el(1:5);
%       Elem4 = [Elem4;el'];  
%     end
%   end

function export_vtk(Nodes,Elems,filename,cvals)

%%% Without reduced, Nodes is a matrix number nodes times 3
% E is a matrix number of elements times 3 or 4
% we can mix triangles and quadrilaterals, but since we need
% a matrix, if we mix trianalges and quadrilaterals, the 4th
% vertex for triangles is set to -1
  if(nargin<4)
    cvals = [];
  end
  
  
  [~,~,ngroups] = size(Nodes);
  
  [rN,cN] = size(Nodes(:,:,1));
  if( rN < cN )
    %% probably wrong format
    for i = 1:ngroups,
      Nodes(:,:,i) = Nodes(:,:,i)';
    end
  end


  
  %% paraviews nodes start with 0, this should fix that
  minE = min(Elems(find(Elems(:)>-1)));
  
  Elems = Elems - minE;
  
  
  fid = fopen(filename, 'wt');
  if (fid < 0)
    error(sprintf("Sorry could not open file %s\n",filename));
  end
  fprintf(fid, "# vtk DataFile Version 3.0\n");
  fprintf(fid, "Comment goes here\n");
  fprintf(fid, "ASCII\n");
  fprintf(fid, "\n");
  fprintf(fid, "DATASET POLYDATA\n");
  
  [n1,n2,n3]=size(Nodes);
  if( n1 == 3 )
    
    for i = 1:n3,
	       % its a little bit tougher if we have multidim matrices
      DNodes(:,:,i) = Nodes(:,:,i)';
    end
    Nodes = DNodes;
    [n1,n2,n3]=size(Nodes);
  end
  
  fprintf(fid, "POINTS   %d double\n", n3*n1);
  
%  for i=1:n3, 
%    for j=1:n2,
%      fprintf(fid,"%f %f %f\n",Nodes(1,j,i),Nodes(2,j,i),Nodes(3,j,i));
%    end
%  end

  % fprintf is really slow in MATLAB so let's use dlmwrite, but
  % surprise that again is not that easy in MATLAB, and additionally
  % the internet says do not dlmwrite in matlab but writematrix
  if is_octave()
    dlmwrite(fid,Nodes,' ');
  else
    fclose(fid);
    %works only for matlab > 2019
    %writematrix(Nodes,filename,'Delimiter',' ','WriteMode','append');
    % should work with older version, however it is not recommended
    dlmwrite(filename,Nodes,'-append','Delimiter',' ');
  end
  
      % see if there are possibly triangles in E that are marked by -1
  trielems = 0;
  [e1,e2,e3]=size(Elems);
  
		     % it is assumed that the nodenumber starts with 0
		     % in nodes
  
  for i = 1:e3,   % maybe we have different faces
    		  % do the triangles first
	% there could be some hidden triangles in there where the last
	% index is -1, remember the elem number is alread removed
    
    
    [~,cE] = size(Elems(:,:,i));
    if(cE == 4) 
      trielems = find(Elems(:,4,i) < 0);
    else
      [rE,~] = size(Elems(:,:,i));
      trielems = [1:rE];
    end
    
    if(isempty(trielems))
      trielems = 0;
      Quads = Elems(:,:,i);
      [quadelems,~] = size(Elems);
    else
      Quads = Elems(:,:,i);
      Quads(trielems,:) = [];
      Etri = Elems(trielems,:,i);
      [trielems,~] = size(Etri);
      if(isempty(Quads))
	quadelems = 0;
      else
	[quadelems,~] = size(Quads);
      end
      [~,cE] = size(Etri);
      if(cE == 4)
				% get rid of the -1
	Etri(:,4) = [];
      end
      
    end

    if ~is_octave()
      fid = fopen(filename,'at');
      if fid < 0
	error( sprintf('Sorry could not open %s\n',filename) );
      end
    end
    fprintf(fid, "POLYGONS %d %d\n", trielems + quadelems, ...
	    4*trielems + 5*quadelems);
    if(trielems > 0)
      if is_octave()
	dlmwrite(fid,[ones(trielems,1) * 3, Etri],' ');
      else
	dlmwrite(filename,[ones(trielems,1) * 3, Etri],'-append',...
		 'Delimiter',' ');
	%writematrix(filename,[ones(trielems,1) * 3, Etri],'WriteMode',
	%	    'append');
      end
    end
%    for j = 1:trielems,
%      fprintf(fid,'3 %d %d %d \n', Etri(j,1),Etri(j,2),Etri(j,3));
%    end
      				%     keyboard
    
    
    
    
    if(~isempty(Quads))
  %	fprintf(fid, "POLYGONS %d %d\n", rows(Quads), 5*rows(Quads));
%      [rQ,~] = size(Quads);
%      for j = 1:rQ,
% 	fprintf(fid,'4 %d %d %d %d\n', Quads(j,1),Quads(j,2),Quads(j,3),Quads(j,4));
%      end
      if is_octave()
	dlmwrite(fid,[ones(quadelems,1)*4,Quads],' ');
      else
	dlmwrite(filename,[ones(quadelems,1)*4,Quads],'-append',...
		 'Delimiter',' ');
	%writematrix([ones(quadelems,1)*4,Quads],filename,'WriteMode',
	%	    'append')
      end
    end	
  end

  
  if(~isempty(cvals))
     % we switched the dimension of the matrices somewhere in the code
    if( ~is_octave() )
      fid = fopen(filename,'at');
      if fid < 0
	error( sprintf("Sorry could not open %s\n",filename) );
      end
    end
    [rN,~] = size(Nodes);
    if( length(cvals) == rN )
      %% we are good to go
      fprintf(fid, "POINT_DATA %d\n", length(cvals));
      text1 = sprintf('SCALARS dataset1 double\n');
      fprintf(fid, text1);
      fprintf(fid, 'LOOKUP_TABLE default\n');	
      fprintf(fid, '%e\n', cvals');
    end
    [rE,~] = size( Elems );
    if( length(cvals) == rE )
      fprintf(fid, "Cell_DATA %d\n", length(cvals));
      fprintf(fid, "SCALARS cell_scalars double 1\n");
      fprintf(fid, "LOOKUP_TABLE default\n");
      fprintf(fid, '%e\n', cvals');
    end
  end
  fclose(fid);
end

function n = is_octave()
  n = exist('OCTAVE_VERSION','builtin');
end
