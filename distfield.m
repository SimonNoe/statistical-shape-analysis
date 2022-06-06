function [distances] = distfield(target,reference)
%DISTFIELD Computes signed distances for each point of target to reference cloud
%   Detailed explanation goes here

% Search for nearest neighbors of each target point
Mdl = KDTreeSearcher(reference);
Idx = knnsearch(Mdl,target);

% Check which reference points are inside target
tri = alphaShape(target);
IsInside = inShape(tri,reference(:,1),reference(:,2),reference(:,3));
% tri = delaunayn(target); % Partition into tetrahedrons
% IsInside = tsearchn(target,tri,reference); % Find closest triangle point if inside element
% IsInside = ~isnan(IsInside); % Convert to logical vector
IsInside = IsInside(Idx); % Correspond to distance vector

% Compute distances
coords = reference(Idx,:);
distances = zeros(length(coords),1);
for i = 1:length(coords)
    if IsInside(i) % Sign distances based on location
        distances(i) = -(norm(coords(i)-target(i)));
    else
        distances(i) = norm(coords(i)-target(i));
    end
end
end

