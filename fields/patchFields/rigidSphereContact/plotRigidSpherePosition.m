a = -225/((4.5e-3)^2);
b = -(9e-3)*a;
acceleration = @(t) a*t.^2 +b*t;
t = (0:0.5:9)*1e-3;
figure(1)
plot(t,acceleration(t))
xlabel('Time [s]')
ylabel('Acceleration [g]');

figure(2)
clf
v0 = -2.5;
radii = (94/2) * 1e-3;
heightBrick = 34.23 * 1e-3;
sphereCenter = radii + heightBrick;
position = @(t,v0, x0) a/12*t.^4 + b/6*t.^3 + v0*t + x0;
plot(t,position(t, v0, sphereCenter));
xlabel('Time [s]')
ylabel('Sphere''s center position [m]');
hold on
plot(t,position(t, v0, sphereCenter)-radii);

% Minimum height brick 34.23/2*1e-3
top = @(t) t*0 + heightBrick;
plot(t,top(t));
axis([0 9e-3 0 0.09])
