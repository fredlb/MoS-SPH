clear all
close all

%particle = [mass, possx, possy, velox, veloy,massdensity, pressure, forcex, forcey ];

[X,Y] = meshgrid(0.5:1/20:1.5,0:1/20:1);
N = size(meshgrid(0:1/20:1,0:1/20:1));
Nx = N(2);
Ny = N(1);
N = Nx*Ny;
X = reshape(X,[N 1]);
Y = reshape(Y,[N 1]);

%{
h = 0.05;
[X,Y] = meshgrid(0:h/1.3:0.5,0:h/1.3:0.5);
N = size(meshgrid(0:h/1.3:0.5,0:h/1.3:0.5));
Nx = N(2);
Ny = N(1);
N = Nx*Ny;
X = reshape(X,[N 1]);
Y = reshape(Y,[N 1]);
%}

poss = zeros(N,2);
poss(:,1) = X;
poss(:,2) = Y;

velocity = zeros(N,2);
Acc = zeros(N,2);
density = zeros(N,1);
pressure = zeros(N,1);

pm = 1; %Particle mass

x = 20;
n = N; 
A = 1*1;
h = sqrt((A*x)/(n*pi));

g = -9.8; 
md0 = 988;
k = 30;
sigma = 0.0728;
dt = 0.005;
cr = 0.2;
vc = 3.5;
l = sqrt(md0/x);

vidObj = VideoWriter('fluid6','Motion JPEG AVI');
vidObj.FrameRate = round(1/dt);
open(vidObj);

figure
yl = 0;
plot(poss(:,1),poss(:,2),'*')
hold on
plot(poss(round(N/2):N,1),poss(round(N/2):N,2),'r*')
hold off
axis tight
set(gca,'nextplot','replacechildren');
%hold on
%figure(2)
%plot(poss(50,1),poss(50,2),'r*')
%hold on 
xlim([0 2])
ylim([yl 2])

time = 1;
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


for t = 0:dt:1
%Density and pressure
Apast = Acc;
%Vpast = velocity;
%Ppast = position;

num = 0;
for i = 1:n
    md = 0;
    for j = 1:n
        r = poss(i,:) - poss(j,:);
        if((r*r') < h^2)
            md = md + pm*Wkernel(r,h,1);
            num = num +  1;
        end
    end
    %num
    density(i) = md;
    pressure(i) = k*(md-md0);
end
%num/(N)

%Force
for i = 1:n
    fp = 0;
    fv = 0;
    ni = 0;
    lci = 0;
    
    fg = g*density(i);
    %fg = 0;
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
    fs = 0;
    if(sqrt(dot(ni,ni)) > l)
        fs = -sigma*lci*(ni/sqrt(dot(ni,ni)));
    end
    fa = 0;
    %{
    if(t < 0.6)
        if(i > n/2)
            fa = [120 0];
        end
    end
    %}
    fp = -density(i)*fp;
    fv = (vc)*fv;
    Acc(i,:) = (fp + fv + fs + [0 fg] + fa);
end
%-------------------Time integration-----------------------------------
for i = 1:n
    a = Acc(i,:)/density(i);
    if(t == 0)
        velocity(i,:) = velocity(i,:) + 0.5*dt*a;
        poss(i,:) = poss(i,:) + dt*velocity(i,:);
    else
        velocity(i,:) = velocity(i,:) + a*dt;
        poss(i,:) = poss(i,:) + velocity(i,:)*dt;
    end
    
    if(poss(i,2) < yl)
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
        poss(i,:) = cp +  d*normal;
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



%figure(1)

currFrame = getframe;
writeVideo(vidObj,currFrame);

plot(poss(:,1),poss(:,2),'*')
xlim([0 2])
ylim([yl 2])
%figure(2)
hold on
plot(poss(round(N/2):N,1),poss(round(N/2):N,2),'r*')
hold off
xlim([0 2])
ylim([yl 2])


%{
Pposs = poss(50,:)
Pvelovity = velovity(50,:)
PAcc = Acc(50,:)
Pdensity = density(50)
Ppressure = pressure(50)
%}
end
close(vidObj)