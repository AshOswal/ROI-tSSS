function [R,EX,EY,EZ,R2,EX2,EY2,EZ2] = ft_getpos(grad, fid, channels)
% Compute the inputs to TSSS code from Fieldtrip representation of MEG
% sensors and headshape.
%_______________________________________________________________________
% Copyright (C) 2017  Vladimir Litvak based on MNE code by Matti Hamalainen
% and Eric Larson
%
% Redistribution and use of the Software in source and binary forms, with or without 
% modification, are permitted for non-commercial use.
% 
% The Software is provided "as is" without warranties of any kind, either express or
% implied including, without limitation, warranties that the Software is free of defects,
% merchantable, fit for a particular purpose. Developer/user agrees to bear the entire risk 
% in connection with its use and distribution of any and all parts of the Software under this license.
% 

fid =  ft_convert_units(fid, 'm');
grad = ft_datatype_sens(grad, 'amplitude', 'T', 'distance', 'm');
     
M = headcoordinates_mne(fid.fid.pnt(1, :), fid.fid.pnt(2, :), fid.fid.pnt(3, :));

fid  = ft_transform_headshape(M, fid);
grad = ft_transform_sens(M, grad);

[dum, meg_chan] = spm_match_str(channels, grad.label);

pos = [];
for i = 1:length(meg_chan)
    pos(i).r0 = grad.coilpos(find(grad.tra(meg_chan(i), :), 1, 'first'), :);
    pos(i).ez = grad.coilori(find(grad.tra(meg_chan(i), :), 1, 'first'), :);
    
    pos(i).r2 = grad.coilpos(find(grad.tra(meg_chan(i), :), 1,'last'), :);
    pos(i).ez2= grad.coilori(find(grad.tra(meg_chan(i), :), 1, 'last'), :);

    [pos(i).ex, pos(i).ey, pos(i).ez] = get_plane_vectors(pos(i).ez);
    [pos(i).ex2, pos(i).ey2, pos(i).ez2] = get_plane_vectors(pos(i).ez2);
    
end

R   =  cat(1, pos(:).r0)';
EX  =  cat(1, pos(:).ex)';
EY  =  cat(1, pos(:).ey)';
EZ  =  cat(1, pos(:).ez)';
R2  =  cat(1, pos(:).r2)';
EX2 =  cat(1, pos(:).ex2)';
EY2 =  cat(1, pos(:).ey2)';
EZ2 =  cat(1, pos(:).ez2)';

function trans = headcoordinates_mne(nas, lpa, rpa)
% Returns the homogenous coordinate transformation matrix
% that converts the specified fiducials in any coordinate system (e.g. MRI)
% into the rotated and translated headccordinate system.
%
% M1 = headcoordinates(nas, lpa, rpa)
%
% The headcoordinate system in Elekta is defined as follows:
% X-axis from the origin towards the RPA point (exactly through)
% Y-axis from the origin towards the nasion (exactly through)
% Z-axis from the origin upwards orthogonal to the XY-plane
% Origin: Intersection of the line through LPA and RPA and a line orthogonal to L passing through the nasion.
%_______________________________________________________________________
% Copyright (C) 2017  Vladimir Litvak based on MNE code by Matti Hamalainen
% and Eric Larson

% $Id: ft_getpos.m 153 2017-04-12 13:19:28Z vladimir $

% ensure that they are row vectors
lpa = lpa(:)';
rpa = rpa(:)';
nas = nas(:)';

diff_1 = nas-lpa;
ex = rpa-lpa;
alpha = dot(diff_1, ex) / dot(ex, ex);

ex = ex/norm(ex);
trans = eye(4);
move = (1-alpha)* lpa + alpha * rpa;
trans(1:3, 4) = move;
trans(1:3, 1) = ex;
ey = nas - move;
ey = ey/norm(ey);
trans(1:3, 2) = ey;
trans(1:3, 3) = cross(ex, ey);

trans = inv(trans);

function [ex, ey, ez] = get_plane_vectors(ez)
% Get two orthogonal vectors orthogonal to ez (ez will be modified).
%_______________________________________________________________________
% Copyright (C) 2017  Vladimir Litvak based on MNE code by Matti Hamalainen
% and Eric Larson

ez = ez(:)';

ez_len = norm(ez);
if ez_len == 0
    error('Zero length normal. Cannot proceed.')
end

if abs(ez_len - abs(ez(3))) < 1e-5
    ex = [1, 0, 0];
else
    ex = zeros(1, 3);
    if ez(2) < ez(3)
        if ez(1) < ez(2)
            ex(1) = 1;
        else
            ex(2) = 1;
        end
    else
        if ez(1) < ez(3)
            ex(1) = 1;
        else
            ex(3) = 1;
        end
    end
end

ez = ez/ez_len;
ex = ex -(ez*ex') * ez;
ex = ex/norm(ex);
ey = cross(ez, ex);