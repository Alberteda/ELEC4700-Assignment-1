% Course: ELEC4700A, Assignment 1
% Name: Oritseserundede Eda (Albert)
% CUID: 100993421
clc
clear all 
clearvars 
set(0,'DefaultFigureWindowStyle','docked')
% Declaring constants to be used in our modelling system. 
electron_mass = 9.10938215e-31; % electron mass
effective_mass = 0.26 * electron_mass; % electron effective mass
temperature = 300; % assumed temperature for thermal velocity in Kelvin
electron_charge  = 1.60217653e-19; % electron charge 
dirac_constant = 1.054571596e-34; % dirac constant 
plank_constant = dirac_constant * 2 * pi; % plank constant 
boltzmann_constant = 1.3806504e-23; % boltzmann constant  
permittivity = 8.854187817e-12; % vacuum permittivity
u_permittivity = 1.2566370614e-6; % vacuum permittivity
c = 299792458; % speed of light         
g = 9.80665; % metres (32.1740 ft) per s²       
am = 1.66053892e-27;
thermal_vel = sqrt((boltzmann_constant * 300) / effective_mass);

% standard deviation for the x and y velocities 
Standdev = thermal_vel/(sqrt(2));
delta_t = 7.5 * 10 ^ -15;
arrayformovsion = zeros(1, 1000);
tmpvel1 = (1:1:1000);
show1 = 300;
wid = 200 * 10 ^ -9; % x-coordinate boundaries
len = 100 * 10 ^ -9; % y-coordinate boundaries

size = 50; % number of electrons in the system 
collisiontime = 0.2 * 10 ^ -12; % time between electron collisions
xpos = rand(1, size) .* wid;
ypos = rand(1, size) .* len;
dft = rand(1,size)*2*pi;
xvel = randn(1,size) .*Standdev;
yvel = randn(1,size) .*Standdev;
nopbox = true;
% construct a while to avoid electrons from starting inside the boundary
while nopbox == true
   inbox = ((xpos <= (1.30 * wid/2) & (xpos >= (0.90 * wid/2))) & ((ypos < (len/3)) | ypos >= (2*len/3)));
   if (sum(inbox) > 0)
       xpos(inbox) = rand(1, sum(inbox)) .* wid;
       ypos(inbox) = rand(1, sum(inbox)) .* len;
   else
       nopbox = false;
   end
       
end

for i = 1:1000
    xopos = xpos;
    xpos = xopos + delta_t*xvel;
    yopos = ypos;
    ypos = yopos + delta_t*yvel;
    side1 = (ypos >= len);
    side2 = (ypos < 0);
    yvel(side1) = -yvel(side1);
    yvel(side2) = -yvel(side2);
    side3 = ((xpos) >= wid);
    side4 = ((xpos) < 0);
    xopos(side3) = xopos(side3) - wid;
    xpos(side3) = xpos(side3) - wid;
    xopos(side4) = xopos(side4) + wid;
    xpos(side4) = xpos(side4) + wid;
    inbox = ((xpos <= (1.30 * wid/2) & (xpos >= (0.90 * wid/2))) & ((ypos < (len/3)) | ypos >= (2*len/3)));
    middle = ((xopos<= (1.30 * wid/2) & (xopos>= (0.90 * wid/2))) & ((yopos > (len/3)) & yopos <= (2*len/3)));
   
    plot (xpos,ypos,'r.');
    pause(0.05)
    hold on
   
    % dimension for the bottle neck boundary
    line([0.90*wid/2 0.90*wid/2], [len 2*len/3]);
    line([1.30*wid/2 1.30*wid/2], [len 2*len/3]);
    line([0.90*wid/2 1.30*wid/2], [len len]);
    line([0.90*wid/2 1.30*wid/2], [2*len/3 2*len/3]);
    line([0.90*wid/2 0.90*wid/2], [0 len/3]);
    line([1.30*wid/2 1.30*wid/2], [0 len/3]);
    line([0.90*wid/2 1.30*wid/2], [0 0]);
    line([0.90*wid/2 1.30*wid/2], [len/3 len/3]);
    side5 = (xpos > 0.90*wid/2);
    side6 = (xpos < 1.30*wid/2);
    side7 = (ypos > 2*len/3);
    side8 = (ypos < len/3);
    side9 = (xpos > 0.90*wid/2);
    side10 = (xpos < 1.30*wid/2);
    side11 = (ypos > 2*len/3);
    side12 = ((side5) & (side6));
    side13 = ((side7)| (side8));
    side14 = ((side9)& (side10));
    side15 = ((side11)& (side12));
    vrms = sqrt ((xvel.^2)+(yvel.^2));
    xvel(inbox&(~middle)) = -xvel(inbox&(~middle));
    yvel(inbox&middle) = -yvel(inbox&middle);
    mft = (thermal_vel * electron_charge)/200;
    mfp = mean(vrms) * mft;


    vrms = sqrt ((xvel.^2)+(yvel.^2));
    show1 = (sqrt (2) * (mean(vrms)^2) * effective_mass )/boltzmann_constant;
    arrayformovsion (1,i)=  show1;
   % generation of values for mappings   
    [xgraph, ygraph] = meshgrid(0:(wid/10):wid, 0:(len/10):len);
    electronmate = zeros(11, 11);
    temperturemate = zeros(11, 11);
    numelectrol = 0;
    totalvel = 0;
  % temperature calculation per mapping 
   for  j = 1:10
    xmin = xgraph(1, j);
    xmax = xgraph(1,j+1);
    for k = 1:10
        ymin = ygraph(k, 1);
        ymax = ygraph(k+1, 1);
        for m = 1:size
            if((xpos(m) > xmin) && (xpos(m) < xmax) && ((ypos(m) > ymin) && ypos(m) < ymax))
                numelectrol = numelectrol + 1;
                electronmate(j, k) = electronmate(j, j) + 1;
                totalvel = totalvel + sqrt((xvel(m) .^ 2) + (xvel(m) .^ 2));
                temperturemate(j, k) = ((sqrt(2)*(totalvel/numelectrol) ^ 2) * effective_mass) / boltzmann_constant;
            end
        end
        totalvel = 0;
        numelectrol = 0;
    end
   end
   
end



fprintf("Mean Free Time is = %f\n", mft);
fprintf("Mean Free Patth is = %f\n", mfp);


figure (1)
title(["Average temperature value = " num2str(show1)]);
xlabel("particle on x axis");
ylabel("particle on y axis");
 

figure (2)
plot(tmpvel1, arrayformovsion);
title('Temperature across Time');
ylabel('Temperature');
xlabel('Time');
plank_constantold on
   

figure(3); 
histogram(vrms, 15);
title('Histogram of Thermal Velocities');
xlabel("Xlabel");
ylabel("Ylabel");

figure(4); 
surf(electronmate);
title('Density Mapping');
xlabel("Xlabel");
ylabel("Ylabel");

figure(5); 
surf(temperturemate);
title('Temperature Mapping');
xlabel("Xlabel");
ylabel("Ylabel");

