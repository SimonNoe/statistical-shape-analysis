function [Vout] = Distortion(Vin,Tin,point,method,norm,sigma,strength)
%DISTORTION Introduces a spherical distortion in Vin at input point.
%   Introduces a spehrical distortion in Vin at input point. Size is
%   determined by strength and distance to input point, range by sigma. It
%   is possible to choose between the distance vector (from point to
%   vertex) and the vertex normal for point displacement.
%   Currently requires user correction for normal vectors pointing inside
%   the shape.

% Gaussian
mu = 0;
% sigma = 10;
max_dist = icdf('normal',0.95,mu,sigma);
% strength = 100;

% Point selection
Mdl = KDTreeSearcher(Vin);
[Idx,dist] = rangesearch(Mdl,point,max_dist);
Idx = Idx{1};
dist = dist{1};
tri = triangulation(Tin,Vin);

% Calculate and apply distortion values
dist_val = pdf('normal',dist,mu,sigma)*strength; % Distortion value
Vout = Vin;
for i = 1:length(Idx)
    if method
        if (Vin(Idx(i),:) == point)
            Vout(Idx(i),:) = Vin(Idx(i),:); % Do nothing
        else
            % Distance vector method
            dist_vect = ((Vin(Idx(i),:)-point)/dist(i))*dist_val(i);
            Vout(Idx(i),:) = Vin(Idx(i),:) + dist_vect;
        end
    else
        if (Vin(Idx(i),:) == point)
            Vout(Idx(i),:) = Vin(Idx(i),:); % Do nothing
        else
            % Vertex normal method
            norm_vect = vertexNormal(tri,Idx(i));
            if norm
                Vout(Idx(i),:) = Vin(Idx(i),:) - norm_vect*dist_val(i);
            else
                Vout(Idx(i),:) = Vin(Idx(i),:) + norm_vect*dist_val(i);
            end
        end
    end
end
    

