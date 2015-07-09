function record=fdtd2ds(source);
% fdtd2ds: finitite differences time domain 2ds

% Define velocity model
nx = 400;
ny = 400;
nx2 = nx+1;
ny2 = ny+1;
vp0 = 2000;
vp = vp0*ones(nx,ny);
            
vp = vp .* vp;
rho = 1500;

% Define source position and receiver position
isx = 200;
isy = 100;
irx = 200;
iry = 300;

% source and fd parameters fd= finite difference
nt = 3000;
fc = 1000;
dx = vp0/(30*fc);
dy = dx;
dt = dx/vp0*1/(sqrt(2)*2*pi);

% Define source-time function
tsour=1/fc;
t = (0:nt-1)*dt;
lsour=tsour/dt;
t0=tsour*1.5;
T0=tsour*1.5;
tau=pi*(t-t0)/T0;
a=4;
if source==0
    fs=(1-a*tau.*tau).*exp(-2*tau.*tau);
else
    fs=source;
end;

% Define PML (perfectly matched layer boundaries) properties
npml = 30;
pmlfac = 50;
pmlexp = 2;
qx = zeros(nx2,ny2);% inutile
qy = zeros(nx2,ny2);% inutile
qx = zeros(nx,ny);
qy = zeros(nx,ny);
for a = 1:npml
    qx(a,:) = pmlfac*(npml-a)^pmlexp; % left explained by the transposition of the matrices when plotted  x is vertical axis in the code and y is horizontal, the trickt to change it is to flip the matrices when plotting
    qx(nx-a+1,:) = pmlfac*(npml-a)^pmlexp; % right
    qy(:,a) = pmlfac*(npml-a)^pmlexp; % top
    qy(:,ny-a+1) = pmlfac*(npml-a)^pmlexp; % bottom
end;
qx = [qx(:,1) qx]; % pour que ça fasse la taille nx2, ny2 
qx = [qx(1,:);qx];
qy = [qy(:,1) qy];
qy = [qy(1,:);qy];



% Initialize fields
px = zeros(nx2,ny2);
py = zeros(nx2,ny2);
ux = zeros(nx2,ny2);
uy = zeros(nx2,ny2);
sfd = zeros(nt,1);

% Main loop
for a = 2:nt

    % Inject source funtion
    px(isx,isy) = px(isx,isy) + dt*0.5*fs(a);
    py(isx,isy) = py(isx,isy) + dt*0.5*fs(a);
    
    % Update px
    diffop = (ux(2:nx2,1:ny) - ux(1:nx,1:ny))/dx;
    pmlop = qx(2:nx2,2:ny2).*px(2:nx2,2:ny2);
    px(2:nx2,2:ny2) = px(2:nx2,2:ny2) - dt*(pmlop + rho*vp.*diffop);
    
    % Update py
    diffop = (uy(1:nx,2:ny2) - uy(1:nx,1:ny))/dy;
    pmlop = qy(2:nx2,2:ny2).*py(2:nx2,2:ny2);
    py(2:nx2,2:ny2) = py(2:nx2,2:ny2) - dt*(pmlop + rho*vp.*diffop);
    
    % Update ux
    diffop = (px(2:nx2,2:ny2) - px(1:nx,2:ny2) + py(2:nx2,2:ny2) - py(1:nx,2:ny2))/dx;
    pmlop = 0.5*(qx(2:nx2,2:ny2)+qx(1:nx,2:ny2)).*ux(1:nx,1:ny);
    ux(1:nx,1:ny) = ux(1:nx,1:ny) - dt/rho*(pmlop + diffop);
    
    % Update uy
    diffop = (px(2:nx2,2:ny2) - px(2:nx2,1:ny) + py(2:nx2,2:ny2) - py(2:nx2,1:ny))/dy;
    pmlop = 0.5*(qy(2:nx2,2:ny2)+qy(2:nx2,1:ny)).*uy(1:nx,1:ny);
    uy(1:nx,1:ny) = uy(1:nx,1:ny) - dt/rho*(pmlop + diffop);
        
    % Display fields or save seismograms
%     if (mod(a,10)==5)
%         clf;
%         imagesc(px' + py'); axis equal; 
%         caxis([-5e-7 5e-7]); 
%         colorbar; 
%         title(sprintf('Time %f ms',a*dt*1000));
%         pause(0.01)
%    end;
   sfd(a) = px(irx-1,iry-1) + py(irx-1,iry-1);
end;
% tfd = t*1000;
% clf;
% plot(tfd,sfd);
% hold on 
record=sfd;



    
