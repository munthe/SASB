function [td, offset] = calc_trm_delays(ele_pos, VS_pos, c,ref_point)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%   ele_pos     N_ele x 3. Position of the elements creating the VS
%   VS_pos      1 x 3. Position of the VS
%   c           Speed of sound.
% OUTPUT
%   td          1 x N_ele. Positive time delays for each trm element creating the VS.
% 	offset      Offset (>=0) due to forcing 'td' to be positive.
% DESCRIPTION
%   Calculates the trm time delays (positive only) based on 1) position of the elements, 2) position of
%   the VS.
% HOW TO USE
%   Just Do It 
% VERSION		
%       v1  06/10/2006
%       v2  10/06/2010 - modified infront calculation to be relative to the
%       closest element and not 0
% AUTHOR    Jacob kortbek
%           Martin Christian Hemmsen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of active elements
N = size(ele_pos,1);

% Initialize 
td      = zeros(1,N);
offset  = 0;

% Vectors from transmitting elements to VS
vect = ones(N,1)*VS_pos - ele_pos;

% Calc distance (positive)
dist = sqrt(sum(vect.^2,2))';


% Convert to range [0:max] (distances relative to minimum distance)
[min_dist index_min] = min(dist);
% dist = dist-min_dist;

% added to remove small numbers
% dist = dist.*dist > eps;

% VS in front or behind transducer
infront = VS_pos(3) > ele_pos(index_min,3);

% new
offset_ref = sqrt(sum((VS_pos-ref_point).^2))/c;

% Convert to time-delays considering whether the VS is in front or behind
% the transducer.

if infront
    dist = -dist;
    offset_ref = -offset_ref;
end
td = dist/c;

% Find delay offset
% min(td) : negative when 'infront' and zero when 'behind'
offset = -min(td); %

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Offset is removed 18/3-2011
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% offset = 0;
% Compensate for delay offset
td = td - offset_ref;
% Find delay offset
% min(td) : negative when 'infront' and zero when 'behind'
offset = -min(td); %
offset = 0;