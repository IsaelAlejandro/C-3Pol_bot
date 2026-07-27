function plot_ellipse3D(major_axis,minor_axis)

 center = [0,0,0];
t = linspace(0, 2*pi, 100);  % Parameter t from 0 to 2π

% Parametric equation of the ellipse
ellipse_points = center' + major_axis' * cos(t) + minor_axis' * sin(t);

% Plot the ellipse
figure;
plot3(ellipse_points(1,:), ellipse_points(2,:), ellipse_points(3,:), 'g', 'LineWidth', 2);
hold on;

% Plot the center
scatter3(center(1), center(2), center(3), 100, 'ro', 'filled');

% Plot the major and minor axis vectors
quiver3(center(1), center(2), center(3), major_axis(1), major_axis(2), major_axis(3), 'm', 'LineWidth', 2);
quiver3(center(1), center(2), center(3), minor_axis(1), minor_axis(2), minor_axis(3), 'b', 'LineWidth', 2);
spin3d=2*[major_axis(2)*minor_axis(3)-minor_axis(2)*major_axis(3), major_axis(3)*minor_axis(1)-minor_axis(3)*major_axis(1),major_axis(1)*minor_axis(2)-minor_axis(1)*major_axis(2)];
quiver3(center(1), center(2), center(3), spin3d(1), spin3d(2), spin3d(3), 'r', 'LineWidth', 2);


% Labels and formatting
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D Ellipse');
grid on; axis equal;
legend('Ellipse', 'Center', 'Major Axis', 'Minor Axis');
axis equal
xlim([-1, 1]); ylim([-1, 1]); zlim([-1, 1])
view(110,27);% Set 3D view
end