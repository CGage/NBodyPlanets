
function [] = solarSystemOfDeath()
close all;
numPlanets = 200;
timeSteps = 100;
planets = init(numPlanets);

for i = 1: timeSteps
plotPlanets(planets, numPlanets);
% [planets, numPlanets] = addPlanet(planets, numPlanets);
% [planets, numPlanets] = planetCollision(planets, numPlanets);
planets.xpos = planets.xpos + randi([-10, 10],1,numPlanets);
    planets.ypos = planets.ypos + randi([-10, 10],1,numPlanets);


end
end

%% This method sets the initial world with numPlanets many planets and 
%returns the planets as a struct.
function planets = init(numPlanets)
    mapHeight = 1000;
    mapWidth = 1000;
    maxPlanetRadius = 25;
    minPlanetRadius = 10;
    maxPlanetDensity = 10;
    minPlanetDensity = 1;
    maxSpeed = 2;
    cmap = hsv(numPlanets);
    planets = struct('radius', 0, 'density', 0,'speed', 0,'direction', 0, ...
        'xpos' , 0, 'ypos', 0, 'colour', [0 0 0]);
    for i = 1: numPlanets
        planets.radius(i) = randi([minPlanetRadius maxPlanetRadius]);
        planets.density(i) = randi([minPlanetDensity, maxPlanetDensity]);
        planets.speed(i) = randi([0, maxSpeed]);
        planets.direction(i) = rand()* (2*pi);
        planets.xpos(i) = randi([10, mapWidth - 10]);
        planets.ypos(i) = randi([10, mapHeight - 10]);
        planets.colour(i, :) = cmap(i, :);
    end
end

%% This method plots the planets stored in the planets struct
function [] = plotPlanets(planets, numPlanets)
scatter(planets.xpos, planets.ypos, planets.radius, planets.colour, 'filled');    
axis([0 1000 0 1000])    
pause(0.01);  
end

%% This function removes planets from the struct when they collide (has bugs)
function [planets, numPlanets] = planetCollision(planets, numPlanets)
planetsToDestroy = zeros(1, numPlanets);
for i = 1: numPlanets
    for j = 1: numPlanets
        dist = sqrt((planets.xpos(i) - planets.xpos(j))^2 + (planets.xpos(i) - planets.xpos(j))^2);
        if ((dist - planets.radius(i) - planets.radius(j) < 0)&& i ~= j)
            planetsToDestroy(i) = 1;
            break;           
        end
    end
end
for i = numPlanets:-1:1
   if (planetsToDestroy == 1)
   planets(i) = []; 
   end

end
numPlanets = numPlanets - sum(planetsToDestroy);
end

%% This function calculates the next time step of the simulation
function [planets, numPlanets] = nextStep(planets, numPlanet)
 
end

%% Adds a new random planet to the planets struct
function [planets, numPlanets] = addPlanet(planets, numPlanets)
    
    i = numPlanets + 1;
    numPlanets = i;
mapHeight = 1000;
    mapWidth = 1000;
    maxPlanetRadius = 25;
    minPlanetRadius = 10;
    maxPlanetDensity = 10;
    minPlanetDensity = 1;
    maxSpeed = 2;
    cmap = hsv(numPlanets);
	planets.radius(i) = randi([minPlanetRadius maxPlanetRadius]);
	planets.density(i) = randi([minPlanetDensity, maxPlanetDensity]);
	planets.speed(i) = randi([0, maxSpeed]);
	planets.direction(i) = rand()* (2*pi);
	planets.xpos(i) = randi([10, mapWidth - 10]);
	planets.ypos(i) = randi([10, mapHeight - 10]);
	planets.colour(i, :) = cmap(i, :);
end

%% This method returns the gravitational force exerted on planet 1 by planet 2
function g = gravity(planet1, planet2)

end