function [xMatrix0,distanceMatrix0,noOfAnchors] = generateProblemFunc(sDim,noisyFac,radiorange,noOfSensors,anchorType,noOfAnchors,randSeed)
% This program is a modification of 'generateProblem.m' in SFSDP software package made by Sunyoung Kim, Masakazu Kojima, Hayato Waki, and Makoto Yamashita.

% A MATLAB program for generating a sensor network localization problem
% Sunyoung Kim, Masakazu Kojima^* and Hayato Waki
% July 28, 2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sDim :    the dimension of the space in which sensors and anchors are
%           located; sDim is either 2 or 3. 
%           When sDim = 2, sensors and anchors are located in [0,1]^2.
%           When sDim = 3, sensors and anchors are located in [0,1]^3.
% noisyFac :    noisy factor 
%               = 0 --- no noise 
%               = \sigma > 0 --- noise with the distribution N(0,\sigma).
% radiorange :  radio range; 
%               If \|x_p - \x_q\| <= radio range, a distance (with noise) is given between x_p and x_q. 
%               If \|x_p - \x_q\| > radio range, no distance is given between x_p and x_q. 
%               Here all sensors are placed in [0,1]^{sDim} randomely.                
% noOfSensors : the number of sensors.
% anchorType :	parameters on how we locate anchors. 
%   anchorType = 0
%       all anchors are placed at grid points on the boundary and interior
%       of [0,1]^{sDim}
%   anchorType = 1
%       all anchors are placed at grid points in the interior of
%       [0,1]^{sDim}
% 	anchorType = 2
%       all anchors are placed randomly in [0,1]^{sDim}
%   anchorType = 3
%       sDim + 1 anchors on the origin and the coordinate axis
%   anchorType = 4
%       sDim + 1 anchors near the center
%   anchorType = 5
%       noOfAnchors are placed randomly + 
%       anchors are placed at the corners of [0,1]^2 or [0,1]^3
%   anchorType = 10
%       no anchor
% noOfAnchors : the number of anchors
% randSeed  :   a random seed number for a random distribution of sensors and
%               anchors when anchorType = 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xMatrix0 :    an sDim x (noOfSensors + noOfAnchors) matrix to represent 
%               locations of sensors and anchors, where noOfSensors sensors
%               are placed in the first noOfSensors columns, and
%               noOfAnchors are place in the last noOfAnchors columns. 
% distanceMatrix0 :  a noOfSensors x (noOfSensors + noOfAnchors) upper
%                   triangular matrix to represent distances from sensors
%                   to sensors and anchors, where the (p,q)th element
%                   d_{pq} denotes the distance from the pth sensor to the qth
%                   sensor or anchor. When noisyFac = 0,
%                       d_{pq} = norm(xMatrix0(:,p) - xMatrix0(:,q))
%                           if norm(xMatrix0(:,p) - xMatrix0(:,q)) <= radiorange,
%                       d_{pq} = 0 
%                           if norm(xMatrix0(:,p) - xMatrix0(:,q)) > radiorange,
%                   When \sigma = noisyFac > 0,
%                       d_{pq} = max(1+epsilon,0.1)*(norm(xMatrix0(:,p) - xMatrix0(:,q))
%                           if norm(xMatrix0(:,p) - xMatrix0(:,q)) <= radiorange,
%                       d_{pq} = 0 
%                           if norm(xMatrix0(:,p) - xMatrix0(:,q)) > radiorange,
%                   Here epsilon is a random number chosen from N(0,\sigma). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xbd = 1.0;

% Checking anchorType and noOfAnchors  --->
if ((anchorType < 0) || (anchorType > 6)) && (anchorType ~= 10)
    error('anchorType should be 0, 1, 2, 3, 4 and 10.');
elseif mod(anchorType,1) > 0
    error('anchorType should be 0, 1, 2, 3, 4 and 10.');
elseif (anchorType == 0) || (anchorType == 1)
    if noOfAnchors < power(2,sDim)
        error('When anchorType = 0 or 1, noOfAnchors has to be n^{sDim} for some integer n >= 2.'); 
    else
        gridForAnchors = power(noOfAnchors,1/sDim); 
        if mod(gridForAnchors,1) > 0
            error('When anchorType = 0 or 1, noOfAnchors has to be n^{sDim} for some integer n >= 2.'); 
        end
    end    
elseif (anchorType == 2)
    if noOfAnchors < sDim+1
       error('noOfAnchors should be not less than sDim+1 when anchorType = 2.');
    end
elseif (anchorType == 3) || (anchorType == 4)
    if noOfAnchors ~= sDim+1
       error('noOfAnchors should be sDim+1 when anchorType = 3 or 4.');
    end 
elseif anchorType == 5
    if noOfAnchors < 0
       error('noOfAnchors should be a nonnegative integer.');
    end 
elseif anchorType == 10
    if noOfAnchors > 0
        error('noOfAnchors should be 0 when anchorType = 10.');
    end
end
% <--- Checking anchorType and noOfAnchors

startingTime = tic; % cputime; 

rand('seed',randSeed);

% placing anchors on grids ---> 
if (anchorType == 0) || (anchorType == 1)
    gridForAnchors = power(noOfAnchors,1/sDim); 
    if noOfAnchors == 1
        noOfAnchors = 1;
        anchorMatrix = ones(sDim,1)*0.5;
        xMatrix0 = [xbd*rand(sDim,noOfSensors),anchorMatrix];
    else % noOfAnchors > 1
        xMatrix0 = xbd*rand(sDim,noOfSensors+noOfAnchors);
        anchorMatrix = sparse(sDim,noOfAnchors);
        h0 = xbd/(gridForAnchors+2);
        if anchorType == 0
            h = xbd / (gridForAnchors - 1);
            origin = zeros(sDim,1);
        else
            h = (xbd-h0) / (gridForAnchors-1);
            origin = (h0/2) * ones(sDim,1);
        end
        pointer = 0;
        if sDim == 2
            for i = 1:gridForAnchors
                for j=1:gridForAnchors
                    pointer = pointer+1;
                    anchorMatrix(1,pointer) = (i-1)*h + origin(1,1);
                    anchorMatrix(2,pointer) = (j-1)*h + origin(2,1);
                end
            end
        elseif sDim == 3
            for i = 1:gridForAnchors
                for j=1:gridForAnchors
                    for k=1:gridForAnchors
                        pointer = pointer+1;
                        anchorMatrix(1,pointer) = (i-1)*h + origin(1,1);
                        anchorMatrix(2,pointer) = (j-1)*h + origin(2,1);
                        anchorMatrix(3,pointer) = (k-1)*h + origin(3,1);
                    end
                end
            end
        end
        xMatrix0 = [xMatrix0(:,1:noOfSensors), anchorMatrix];
    end
% <--- placing anchors on grids
% placing anchors randomly --->
elseif anchorType == 2
    xMatrix0 = xbd*rand(sDim,noOfSensors+noOfAnchors);
% <--- placing anchors randomly
% placing anchors on the coordinate axes ---> 
elseif anchorType == 3
    noOfAnchors = sDim + 1;
    xMatrix0 = xbd*rand(sDim,noOfSensors+noOfAnchors);
    anchorMatrix = sparse(sDim,noOfAnchors); 
    pointer = 1; 
    anchorMatrix(:,pointer) = zeros(sDim,1); 
    for i=1:sDim
        pointer = pointer+1; 
        anchorMatrix(:,pointer) = zeros(sDim,1);
        anchorMatrix(i,pointer) = xbd*0.5;
    end
    xMatrix0 = [xMatrix0(:,1:noOfSensors), anchorMatrix];
% <--- placing anchors on the coordinate axes
% placing anchors in a small box near the center ---> 
elseif anchorType == 4
    noOfAnchors = sDim + 1;
    xMatrix0 = xbd*rand(sDim,noOfSensors+noOfAnchors);
    pointer = 1;
    anchorMatrix(:,pointer) = 0.5*xbd*ones(sDim,1); 
    for i=1:sDim
        pointer = pointer+1; 
        anchorMatrix(:,pointer) = 0.5*xbd*ones(sDim,1);
        anchorMatrix(i,pointer) = 0.6*xbd;
    end
    xMatrix0 = [xMatrix0(:,1:noOfSensors), anchorMatrix];
elseif anchorType == 5
    xMatrix0 = xbd*rand(sDim,noOfSensors);
    if sDim == 2
        for i=0:1
            for j=0:1
                xMatrix0 = [xMatrix0,xbd*[i;j]];
            end
        end
    elseif sDim == 3
        for i=0:1
            for j=0:1
                for k=0:1
                    xMatrix0 = [xMatrix0,xbd*[i;j;k]];
                end
            end
        end
    end
    xMatrix0 = [xMatrix0,xbd*rand(sDim,noOfAnchors)]; 
    if sDim == 2
        noOfAnchors = noOfAnchors + 4; 
    elseif sDim == 3
        noOfAnchors = noOfAnchors + 8;         
    end
elseif anchorType == 10
    noOfAnchors = 0;    
    xMatrix0 = xbd*rand(sDim,noOfSensors);    
end
% <--- placing anchors in a small box near the center

randSeed = randSeed + 13;

[distanceMatrix0] = computeDistance2(noOfSensors,xMatrix0,radiorange,noisyFac,randSeed); 

% fprintf('## elapsed time for generating a sensor network problem = %8.2f\n',toc(startingTime)); 

if 0 == 1
    sensorDistMat = distanceMatrix0(:,1:noOfSensors)+distanceMatrix0(:,1:noOfSensors)'+(noOfSensors+1)*speye(noOfSensors,noOfSensors);
    permutation = symamd(sensorDistMat);
    UMat = chol(sensorDistMat(permutation,permutation));
    figure(1);
    spy(UMat + UMat');
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [distanceMatrix0] = computeDistance(noOfSensors,xMatrix0,radiorange,noisyFac,randSeed)
% 
% sDim = size(xMatrix0,1);
% randn('seed',randSeed);
% noOfAnchors = size(xMatrix0,2) - noOfSensors;
% distanceMatrix0 = sparse(noOfSensors,noOfSensors+noOfAnchors);
% 
% startingTime = tic; 
% 
% for r = noOfSensors+1:noOfSensors+noOfAnchors
%     for p=1:noOfSensors
%         d0 = norm(xMatrix0(:,p)-xMatrix0(:,r));
%         if d0 <= radiorange
%             rate = max([1+randn(1,1)*noisyFac,0.1]);
%             distanceMatrix0(p,r) = d0*rate;
%         end
%     end
% end
% 
% fprintf('##0 %6.1f\n',toc(startingTime));
% startingTime = tic; 
% 
% idx = [1:noOfSensors];
% zeroPone = repmat(0.1,1,noOfSensors+noOfAnchors);
% for q = idx
%     tmpIdx = repmat(q,1,q-1);
%     tmpMat = xMatrix0(:,1:q-1)-xMatrix0(:,tmpIdx);
%     tmpMat = tmpMat.*tmpMat;
%     tmpMat = sum(tmpMat,1);
%     tmpMat = sqrt(tmpMat);
%     I = find(tmpMat <= radiorange);
%     if ~isempty(I)
%         s = length(I);
%         tmp = [1+randn(1,s)*noisyFac;zeroPone(1,1:s)];
%         tmpMat(1,I) = tmpMat(1,I).* max(tmp,[],1);
%         distanceMatrix0(I,q) = tmpMat(1,I)';
%     end
% end
% 
% fprintf('##1 %6.1f\n',toc(startingTime));
% 
% return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [distanceMatrix0] ...
    = computeDistance2(noOfSensors,xMatrix0,radiorange,noisyFac,randSeed)
noOfAnchors = size(xMatrix0,2) - noOfSensors;

if noisyFac == 0 
%    startingTime = tic; 
    rowIdx = [];
    colIdx = [];
    value = [];
    for r=noOfSensors+1:noOfSensors+noOfAnchors
        tmpIdx = repmat(r,1,noOfSensors);
        tmpMat = xMatrix0(:,1:noOfSensors)-xMatrix0(:,tmpIdx);
        normVector = sqrt((sum(tmpMat .* tmpMat,1)));
        rowAdd = find(normVector < radiorange);
        lenRowAdd = length(rowAdd);
        if lenRowAdd > 0
            rowIdx = [rowIdx,rowAdd];
            colIdx = [colIdx,repmat(r-noOfSensors,1,lenRowAdd)];
            %        valueAdd = normVector(rowAdd) .* (ones(1,lenRowAdd)+randn(1,s)*noisyFac);
            value = [value,normVector(rowAdd)];
        end
    end
    senToAnchDistMat = sparse(rowIdx,colIdx,value,noOfSensors,noOfAnchors);
%     fprintf('##0 %6.1f\n',toc(startingTime));
%     startingTime = tic;
    distanceMatrix0 = sparse(noOfSensors,noOfSensors);
    idx = [1:noOfSensors];
    for q = idx
        tmpIdx = repmat(q,1,q-1);
        tmpMat = xMatrix0(:,1:q-1)-xMatrix0(:,tmpIdx);
        tmpMat = tmpMat.*tmpMat;
        tmpMat = sum(tmpMat,1);
        tmpMat = sqrt(tmpMat);
        I = find(tmpMat <= radiorange);
        if ~isempty(I)
%            s = length(I);
%            tmp = [1+randn(1,s)*noisyFac;zeroPone(1,1:s)];
%            tmpMat(1,I) = tmpMat(1,I).* max(tmp,[],1);
            distanceMatrix0(I,q) = tmpMat(1,I)';
        end
    end
%    fprintf('##1 %6.1f\n',toc(startingTime));
else
    randn('seed',randSeed);
    zeroPone = repmat(0.1,1,noOfSensors+noOfAnchors);
%     startingTime = tic; 
    rowIdx = [];
    colIdx = [];
    value = [];
    for r=noOfSensors+1:noOfSensors+noOfAnchors
        tmpIdx = repmat(r,1,noOfSensors);
        tmpMat = xMatrix0(:,1:noOfSensors)-xMatrix0(:,tmpIdx);
        normVector = sqrt((sum(tmpMat .* tmpMat,1)));
        rowAdd = find(normVector < radiorange);
        lenRowAdd = length(rowAdd);
        if lenRowAdd > 0
            rowIdx = [rowIdx,rowAdd];
            colIdx = [colIdx,repmat(r-noOfSensors,1,lenRowAdd)];
            rateVector = max([ones(1,lenRowAdd)+randn(1,lenRowAdd)*noisyFac;zeroPone(1,1:lenRowAdd)],[],1);
            valueAdd = normVector(rowAdd) .* rateVector;
            value = [value,valueAdd];
        end
    end
    senToAnchDistMat = sparse(rowIdx,colIdx,value,noOfSensors,noOfAnchors);
%     fprintf('##0 %6.1f\n',toc(startingTime));
%     startingTime = tic;
    distanceMatrix0 = sparse(noOfSensors,noOfSensors);
    idx = [1:noOfSensors];
    for q = idx
        tmpIdx = repmat(q,1,q-1);
        tmpMat = xMatrix0(:,1:q-1)-xMatrix0(:,tmpIdx);
        tmpMat = tmpMat.*tmpMat;
        tmpMat = sum(tmpMat,1);
        tmpMat = sqrt(tmpMat);
        I = find(tmpMat <= radiorange);
        if ~isempty(I)
            s = length(I);
            tmp = [1+randn(1,s)*noisyFac;zeroPone(1,1:s)];
            tmpMat(1,I) = tmpMat(1,I).* max(tmp,[],1);
            distanceMatrix0(I,q) = tmpMat(1,I)';
        end
    end
%     fprintf('##1 %6.1f\n',toc(startingTime));
end

distanceMatrix0 = [distanceMatrix0, senToAnchDistMat];


% rowIdx = [];
% colIdx = [];
% value = [];
% 
% for p=1:noOfSensors-1
%     tmpIdx = repmat(p,1,noOfSensors-p);
%     tmpMat = xMatrix0(:,p+1:noOfSensors)-xMatrix0(:,tmpIdx);
%     normVector = sqrt((sum(tmpMat .* tmpMat,1))); 
%     rowAdd = find(normVector < radiorange); 
%     if ~isempty(rowAdd) 
%         rowIdx = [rowIdx,rowAdd+p];
%         colIdx = [colIdx,repmat(p,1,length(rowAdd))];
%         value = [value,normVector(rowAdd)];
%     end        
% end
% distanceMatrix0= [sparse(rowIdx,colIdx,value,noOfSensors,noOfSensors)',senToAnchDistMat];

return

