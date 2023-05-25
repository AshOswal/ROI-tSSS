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
