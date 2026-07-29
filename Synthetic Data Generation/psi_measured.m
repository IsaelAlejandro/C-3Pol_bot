function[angle]=psi_measured(a,b)
% Normalize input vectors
a = a / norm(a);
b = b / norm(b);

% Define a reference normal vector (plane normal)
% Option 1: use the cross product of a and b (good for defining plane)
n = cross(a, b);  % This is perpendicular to the plane of a and b

% Optionally normalize the normal vector
n = n / norm(n);

% Compute angle using atan2
dot_ab = dot(a, b);
cross_ab = cross(a, b);
angle = mod(atan2(norm(cross_ab), dot_ab),2*pi);  % Gives angle in [0, pi]

% Determine direction (sign) using the reference normal
sign_dir = sign(dot(n, cross_ab)) % +1 or -1 depending on orientation
if sign_dir < 0
    angle = mod(2*pi - angle,2*pi);
end