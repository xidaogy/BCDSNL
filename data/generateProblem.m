%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is a component of the software package BCDSNL 
% Copyright (C) 2022 Mitsuhiro Nishijima and Kazuhide Nakata
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars

sDim = 2;
randSeed = 1;
noOfSensors = 1000;
noisyFac = 0;
radiorange = 0.1;
noOfAnchors = noOfSensors/10;
anchorType = 2;
[xMatrix0,distanceMatrix0,noOfAnchors] = generateProblemFunc(sDim,noisyFac,radiorange,noOfSensors,anchorType,noOfAnchors,randSeed);

if noisyFac == 0
    ubdForSenToAnchorEdge = sDim+1;
else
    ubdForSenToAnchorEdge = (sDim+1)*2;
end
minDegree = sDim + 2;
distanceMatrix = distanceMatrix0;
saDistMat = distanceMatrix(:,noOfSensors+1:noOfSensors+noOfAnchors);

[S,ssaDistance] = calcS(distanceMatrix, noOfSensors, noOfAnchors);
noOfssEdges = nnz(distanceMatrix(:,1:noOfSensors));
noOfsaEdges = nnz(distanceMatrix(:,noOfSensors+1:noOfSensors+noOfAnchors));

for i = 1:noOfSensors
    IncidenceIdx{i} = find(S(i,:));
end

filename = "/data/d"+sDim+"_s"+noOfSensors+"_a"+noOfAnchors+"_n"+sprintf('%1.0e',noisyFac)+"_r1e-1"+'.mat';
save(filename,'distanceMatrix0','IncidenceIdx','minDegree','noisyFac','noOfAnchors','noOfsaEdges','noOfSensors','noOfssEdges','radiorange','S','saDistMat','sDim','ssaDistance','ubdForSenToAnchorEdge','xMatrix0');