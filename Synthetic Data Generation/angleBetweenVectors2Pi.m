function [v4_rot] = angleBetweenVectors2Pi(a,b,b0)


v1 = b    ;      % random vector
v1 = v1 / norm(v1)  ;    % normalize

% Generate another vector orthogonal to v1
                         % another random vector
v2 = a    ;       % remove component along v1
v2 = v2 / norm(v2) ;

                         % another random vector
v4 = b0 ;          % remove component along v1
v4 = v4 / norm(v4) ;
% Assume v1 and v2 are orthogonal vectors in 3D
v1 = v1 / norm(v1);  % normalize
v2 = v2 / norm(v2);  % normalize

v3 = cross(v1, v2);  % v3 is orthogonal to both v1 and v2
v3 = v3 / norm(v3);  % ensure orthonormality

% STEP 3: Construct rotation matrix to map [v1 v2 v3] → [x z y]
% Target axes (columns of identity matrix permuted)
target1 = [1; 0; 0];  % x-axis
target2 = [0; 0; 1];  % z-axis
target3 = -cross(target1, target2);  % y-axis (right-hand rule)

% Original basis formed by v1, v2, v3
original_basis = [v1, v2, v3];       % 3x3
target_basis   = [target1, target2, target3];  % 3x3

% STEP 4: Compute rotation matrix
R = target_basis * original_basis';  % rotates original basis to target

% STEP 5: Apply rotation
v1_rot = R * v1
v2_rot = R * v2
v4_rot = R * v4
   
end
