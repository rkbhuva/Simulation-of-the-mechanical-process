% Hand in 3 : Collisions of particles, comparison hard-sphere and DEM
% Name : Ridham Kishorbhai Bhuva
% Matr. No. : 245921
% Course : CEE

%% PLEASE DO COPY SUB ANSWERS IN NEW FILE AND RUN IT!!! DO COPY FUNCTIONS TOO!!!

%% ANSWER 2 : Hard-sphere algorithm - Function can be found at end of the code

clear all;
clc;

% Mass of particles in kg
m1 = 0.05;
m2 = 0.05;

% Diameter of particles in m
d1 = 1;
d2 = 1;

% Coefficient of restitution
e = 1; % change e according e = 1 or 0.75

% Initialize position and velocity
x1 = [0, 0, 0];
x2 = [1.1, 1.3, 0];
v1i = [0, 0, 0];
v2i = [-1, -1, 0];


% Determine collision time of given two particles
[tMin] = collisionTime(d1, d2, x1, x2, v1i, v2i);

% New position of particles
x1New = x1 + v1i*tMin;
x2New = x2 + v2i*tMin;

% Determine post collision velocity from new position
[v1f, v2f] = collision(m1, m2, x1New, x2New, v1i, v2i, e);

tMin
v1f
v2f

%% ANSWER 3 : Matlab function to calculate normal vectore and overlap can be find at the end of the code 

%% ANSWER 4 : To use Leapfrog algorithm, normal force need to be include in the code to calculate acceleration


clear all;
clc;

% Mass of particles in kg
m = 0.05; % m1 = m2 = m

% Diameter of particles in m
d1 = 1;  
d2 = 1;

% Initial position of particles in m
x1i = [0 0 0];
x2i = [1.1 1.3 0];

% Initial velocity of particles in m/s
v1i = [0 0 0];
v2i = [-1 -1 0];

% Allocating time and time step
tstart = 0;
tend = 2; % should be minimum 1.9 for compelete collision
dt = 1.0e-4;
t = tstart:dt:tend;

% Number of steps
N = length(t);

% Spring constant in N/m
kn = 750;

% Preallocate arrays
x2 = zeros(N+1, 3);
v2 = zeros(N+1, 3);
n = zeros(N, 3);
overlap = zeros(N, 1);
force = zeros(N, 3);
acceleration = zeros(N, 3);
vp2 = zeros(N, 3);

% Initialize position and velocity 
x1(1, :) = x1i;
x2(1, :) = x2i;
v1(1, :) = v1i;
v2(1, :) = v2i;


% As mentioned in question, moves particle 2 with fixed time step dt
% 
% However, acceleration sholud be zero till the collision
% 
% After collision particle 1 will be still fixed 
% as mentioned particle 2 moves with fixed time step dt
%
% e = 1 ideal considered, therefore kn will be constant (knLoad = knUnload)
%
% NOTE : If particle 1 should also moves after collision then please consider 
% code of ANSWER 5 with e = 1!!! 
% (Then ANSWER 4 & 5 will be in one code ANSWER 5)

% v1 & v2 are corrector velocities
% vp1 & vp2 are predictor velocities

for i = 1:N
    
    % Determine normal vector and overlap
    [n(i, :), overlap(i)] = funA(d1, d2, x1, x2(i, :));
    
    % Force and acceleration
    force(i, :) = -kn * overlap(i) .* n(i, :);
    acceleration(i, :) = force(i, :) / m;
    
    % Predictor velocity for particle 2    
    vp2(i, :) = v2(i, :) + dt * acceleration(i, :);

    % New position of particle 2 after time step    
    x2(i+1, :) = x2(i, :) + dt * vp2(i, :) + 0.5 * dt * dt * acceleration(i, :);
        
    % Determine new velocity of particles using implicite acceleration
    [n_new, overlap_new] = funA(d1, d2, x1, x2(i+1, :));
    accelerationNew = -kn * overlap_new .* n_new / m;

    % Corrector velocity for particle 2
    v2(i+1, :) = vp2(i, :) + 0.5 * dt * (acceleration(i, :) + accelerationNew);
    
    
end

%% ANSWER 5 : As mentioned in the question, " apply normal force on both particles"
% In this sub question force will be apply on both particles
%
% Subscript 2 will be use for this sub task as x21,x22,v21 & v22

% v21 & v22 are corrector velocities
% vp21 & vp22 are predictor velocities

clear all;
clc;

% Mass of particles in kg
m = 0.05; % m1 = m2 = m

% Diameter of particles in m
d1 = 1;  
d2 = 1;

% Initial position of particles in m
x1i = [0 0 0];
x2i = [1.1 1.3 0];

% Initial velocity of particles in m/s
v1i = [0 0 0];
v2i = [-1 -1 0];

% Allocating time and time step
tstart = 0;
tend = 2; % should be minimum 1.9 for complete collision
dt = 1.0e-4;
t = tstart:dt:tend;

% Number of steps
N = length(t);

% Coefficient of restitution
e = 0.75; % change e according e = 1 or 0.75

% Spring constant for loading and unloading
knLoad = 750;
knUnload = knLoad * (e^2);

% Preallocate arrays
x21 = zeros(N+1, 3);
x22 = zeros(N+1, 3);
v21 = zeros(N+1, 3);
v22 = zeros(N+1, 3);
n2 = zeros(N, 3);
overlap2 = zeros(N, 1);
v12 = zeros(N, 3);
force2 = zeros(N, 3);
acceleration2 = zeros(N, 3);
vp21 = zeros(N, 3);
vp22 = zeros(N, 3);

% Initialize position and velocity 
x21(1, :) = x1i;
x22(1, :) = x2i;
xp21(1, :) = x1i;
xp22(1, :) = x2i;
v21(1, :) = v1i;
v22(1, :) = v2i;
vp21(1, :) = v1i;
vp22(1, :) = v2i;



for i = 1:N
    
    % Determine normal vector and overlap
    [n2(i, :), overlap2(i)] = funA(d1, d2, x21(i, :), x22(i, :));

    % Determine loading and unloading 'kn' based on overlap
    v12(i, :) = v21(i, :) - v22(i, :);

          if dot(v12(i, 1:3),n2(i, 1:3)) < 0 
            kn = knLoad;
          else
            kn = knUnload;
          end
                
    % Force and acceleration
    force2(i, :) = -kn * overlap2(i) .* n2(i, :);
    acceleration2(i, :) = force2(i, :) / m;
    
    % Predictor velocity for particles    
    vp21(i+1, :) = v21(i, :) + dt * acceleration2(i, :);
    vp22(i+1, :) = v22(i, :) + dt * acceleration2(i, :);
    
    % New position of particles after time step    
    xp21(i+1, :) = x21(i, :) + dt * vp21(i, :) + 0.5 * dt * dt * acceleration2(i, :);
    xp22(i+1, :) = x22(i, :) + dt * vp22(i, :) + 0.5 * dt * dt * acceleration2(i, :);
    
    % Determine new velocity of particles using implicite acceleration
    [n2_new, overlap2_new] = funA(d1, d2, xp21(i+1, :), xp22(i+1, :));
    accelerationNew2 = -kn * overlap2_new .* n2_new / m;
    
    % Corrector velocity for particles
    v21(i+1, :) = v21(i, :) + 0.5 * dt * (acceleration2(i, :) + accelerationNew2);
    v22(i+1, :) = v22(i, :) + 0.5 * dt * (acceleration2(i, :) + accelerationNew2);

    x21(i+1, :) = x21(i, :) + dt * v21(i, :) + 0.5 * dt * dt * acceleration2(i, :);
    x22(i+1, :) = x22(i, :) + dt * v22(i, :) + 0.5 * dt * dt * acceleration2(i, :);
   
end

%% ANSWER 5 (a) : Find the contact point coordinates before the collision

contactIndices = find(overlap2 > 0); % Indices where overlap > 0

% Choose the first contact time step for demonstration
contactIndex = contactIndices(1); % First contact time step
    
% Positions of the particles at the contact time step
x1_contact = x21(contactIndex, :); 
x2_contact = x22(contactIndex, :); 
    
% Normal vector at the contact time step
n_contact = n2(contactIndex, :); 
    
% Calculate the contact point
Contact_point_coordinates_before_collision = x1_contact + (d1 / 2) * n_contact % Print in command window

%% ANSWER 5 (b & c) : The post-collision velocities of particles

% Find collision start and end times
collisionStartIndex = find(overlap2 > 0, 1, 'first'); % First index where overlap > 0
collisionEndIndex = find(overlap2 > 0, 1, 'last'); % Last index where overlap > 0

v1_postCollision = v21(collisionEndIndex + 1, :) % Velocity of particle 1
v2_postCollision = v22(collisionEndIndex + 1, :) % Velocity of particle 2
% Print in command window

%% ANSWER 5 (d) : The duration of the collision

% Calculate collision duration
collisionStartTime = collisionStartIndex * dt;
collisionEndTime = collisionEndIndex * dt;

collisionDuration = collisionEndTime - collisionStartTime % Print in command window


%% ANSWER 3 : Function to Calculate normal vector and overlap
% This function needed for ANSWER 4 & 5

function [n, overlap] = funA(d1, d2, x1, x2)

    % Relative position of particles
    x12 = x1 - x2;
    % Normal vector
    n = x12 / sqrt(dot(x12,x12));
    % Mean diameter (distace between two particle while in contact)
    dMean = (d1 + d2)/2;
    % Distance between two particle in space
    x12Dist = sqrt( sum(x12.^2) );
    % Determine overlap
    if x12Dist < dMean
        overlap = dMean - x12Dist;
    else
        overlap = 0;
    end
    
end
%% ANSWER 2 functions : Function to Calculate collision time and post collision states
% This functions needed for ANSWER 2

function [tCol] = collisionTime(d1, d2, x1, x2, v1i, v2i)

    % Relative position and velocity of particles
    x12 = x1 - x2;
    v12 = v1i - v2i;

    % Time calculation
    t1 = ( -dot(x12, v12) - sqrt( dot(x12,v12)^2 - dot(v12,v12)*(dot(x12,x12) - (0.5*d1 + 0.5*d2)^2)) ) / dot(v12,v12);
    t2 = ( -dot(x12, v12) + sqrt( dot(x12,v12)^2 - dot(v12,v12)*(dot(x12,x12) - (0.5*d1 + 0.5*d2)^2)) ) / dot(v12,v12);
    
    % Min time / collision time
    if (t1 > 0 && t2 > 0)
        tCol = min(t1,t2);
    elseif (t1 < 0 && t2 < 0)
        tCol = -1;
    else
        if t1 > 0
            tCol = t1;
        else
            tCol = t2;
        end
    end

end

function [v1f, v2f] = collision(m1,m2, x1,x2, v1i,v2i, e)

    % Pre calculations
    mStar = m1*m2/(m1+m2);
    x12 = x1 - x2;
    n = x12 / sqrt(dot(x12,x12));
    v12i = v1i - v2i;
    
    % Velocity after collision
    if e == 1  
        v1f = (m1-m2)/(m1+m2) .* v1i + 2*m2/(m1+m2) .* v2i;
        v2f = 2*m1/(m1+m2) .* v1i + (m2-m1)/(m1+m2) .* v2i;
    else
        J = -(1+e) .* dot(v12i,n) * mStar .* n;    
        v1f = v1i + J/m1;
        v2f = v2i - J/m2;
    end
end
