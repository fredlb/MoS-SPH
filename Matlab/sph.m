%% SPH
%% Partiles
clear all
close all

n = 20;
xbegin = 0.5;
ybegin = 0;
[X,Y] = meshgrid(xbegin:1/n:(xbegin+1),ybegin:1/n:(ybegin+1));
N = size(meshgrid(0:1/n:1,0:1/n:1));
Nx = N(2);
Ny = N(1);
n = Nx*Ny;
X = reshape(X,[n 1]);
Y = reshape(Y,[n 1]);

poss = zeros(n,2);
poss(:,1) = X;
poss(:,2) = Y;

velocity = zeros(n,2);
Acc = zeros(n,2);
density = zeros(n,1);
pressure = zeros(n,1);

%% Constants 
pm = 1;                     %Particle mass (changes later)

x = 20;                     %Average number of particles in kernel
A = 1*1;                    %Area
h = sqrt((A*x)/(n*pi));     %kernel size

dt = 0.005;
k = 30;

g = -9.8;                   
md0 = 988;
surfaceTension = 0.0728;
surfaceLimit = sqrt(md0/x);
cr = 0.2;                   %Dampening on border cases
viscosityConstant = 3.5;
runTime = 10;


%% Video recording
vidObj = VideoWriter('fluid6','Motion JPEG AVI');
vidObj.FrameRate = round(1/dt);
open(vidObj);

%Write particles
figure
plot(poss(:,1),poss(:,2),'*')
hold on
plot(poss(round(n/2):n,1),poss(round(n/2):n,2),'r*')
hold off
axis tight
set(gca,'nextplot','replacechildren');
xlim([0 2])
ylim([0 2])

%% Adjust mass after desierd density
md = 0;
for i = 1:n
    for j = 1:n
        r = poss(i,:) - poss(j,:);
        if((r*r') < h^2)
            md = md + pm*Wkernel(r,h,1);
        end
    end
end
amd = md/n;
pm = (amd*md0)/(amd*amd);

%% Begin time loop
for t = 0:dt:runTime
%% Density and pressure
for i = 1:n
    md = 0;
    for j = 1:n
        r = poss(i,:) - poss(j,:);
        if((r*r') < h^2)
            md = md + pm*Wkernel(r,h,1);
        end
    end
    density(i) = md;
    pressure(i) = k*(md-md0);
end

%% Forces
for i = 1:n
    fp = 0;             %Pressure force
    fv = 0;             %Viscosity force
    ni = 0;             %Normal
    lci = 0;            %laplacian of ci (gradient of normal)
    fs = 0;             %Surface tension force
   
    fg = g*density(i);  %Gravity
    
    for j = 1:n
        r = poss(i,:) - poss(j,:);
        if((r*r')< h^2)
            ni = ni + (pm/density(j))*Wkernel(r,h,4);
            lci = lci + (pm/density(j))*Wkernel(r,h,5);
        if r ~= 0
            fp = fp +(((pressure(i)/(density(i)^2)))+((pressure(j)/(density(j)^2))))*pm*Wkernel(r,h,2);
            u = velocity(j,:) - velocity(i,:);
            fv = fv + u*(pm/density(j))*Wkernel(r,h,3);
        end
        end
    end
    
    if(sqrt(dot(ni,ni)) > surfaceLimit)
        fs = -surfaceTension*lci*(ni/sqrt(dot(ni,ni)));
    end
    
    fp = -density(i)*fp;
    fv = (viscosityConstant)*fv;
    
    %If you want to play with starting force
    fa = 0;
    %{
    if(t < 0.6)
        if(i > n/2)
            fa = [120 0];
        end
    end
    %}
    Acc(i,:) = (fp + fv + fs + [0 fg] + fa);
end
%% Time integration
for i = 1:n
    a = Acc(i,:)/density(i);
    if(t == 0)
        velocity(i,:) = velocity(i,:) + 0.5*dt*a;
        poss(i,:) = poss(i,:) + dt*velocity(i,:);
    else
        velocity(i,:) = velocity(i,:) + a*dt;
        poss(i,:) = poss(i,:) + velocity(i,:)*dt;
    end

%% Colision detection and response
    if(poss(i,2) < 0)
        x = poss(i,:);
        cp = [poss(i,1) 0];
        d = sqrt(dot(cp-x,cp-x));
        normal = [0 1];
       
        ui = velocity(i,:);
        poss(i,:) =cp + d*normal;
        velocity(i,:) =ui - (1 + cr*(d/(dt*sqrt(dot(ui,ui)))))*(ui*normal')*normal;
    end
    
    if(poss(i,1) > 2)
        x = poss(i,:);
        cp = [2 poss(i,2)];
        d = sqrt(dot(cp-x,cp-x));
        normal = [-1 0];
       
        ui = velocity(i,:);
        poss(i,:) =cp +  d*normal;
        velocity(i,:) =ui - (1 + cr*(d/(dt*sqrt(dot(ui,ui)))))*(ui*normal')*normal;
    end
    
    if(poss(i,1) < 0)
        x = poss(i,:);
        cp = [0 poss(i,2)];
        d = sqrt(dot(cp-x,cp-x));
        normal = [1 0];
       
        ui = velocity(i,:);
        poss(i,:) =cp + d*normal;
        velocity(i,:) =ui - (1 + cr*(d/(dt*sqrt(dot(ui,ui)))))*(ui*normal')*normal;
    end
    
end

%% Make pretty movie
currFrame = getframe;
writeVideo(vidObj,currFrame);

plot(poss(:,1),poss(:,2),'*')
xlim([0 2])
ylim([0 2])
hold on
plot(poss(round(n/2):n,1),poss(round(n/2):n,2),'r*')
hold off
xlim([0 2])
ylim([0 2])
end
close(vidObj)