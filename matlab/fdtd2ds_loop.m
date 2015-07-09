function record=fdtd2ds_loop(source);
% fdtd2ds: finitite differences time domain 2ds

% Define velocity model
nx = 400;
ny = 400;
nx2 = nx+1;
ny2 = ny+1;
vp0 = 2000;
vp = vp0*ones(nx2,ny2);
            
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


dx = 1/dx;
dy = 1/dy;
oneOverRho= 1/rho;

% Main loop
for a = 2:nt

    % Inject source funtion
    px(isx,isy) = px(isx,isy) + dt*0.5*fs(a);
    py(isx,isy) = py(isx,isy) + dt*0.5*fs(a);
    
    for j=1:ny
      for i=1:ny
    
        % Update px
        diffop = (ux(i+1,j) - ux(i,j))*dx;
        pmlop = qx(i+1,j+1)*px(i+1,j+1);
        px(i+1,j+1) = px(i+1,j+1) - dt*(pmlop + rho*vp(i+1,j+1)*diffop);

        % Update py
        diffop = (uy(i,j+1) - uy(i,j))*dy;
        pmlop = qy(i+1,j+1)*py(i+1,j+1);
        py(i+1,j+1) = py(i+1,j+1) - dt*(pmlop + rho*vp(i+1,j+1)*diffop);

        % Update ux
        diffop = (px(i+1,j+1) - px(i,j+1) + py(i+1,j+1) - py(i,j+1))*dx;
        pmlop = 0.5*(qx(i+1,j+1)+qx(i,j+1))*ux(i,j);
        ux(i,j) = ux(i,j) - dt*oneOverRho*(pmlop + diffop);

        % Update uy
        diffop = (px(i+1,j+1) - px(i+1,j) + py(i+1,j+1) - py(i+1,j))*dy;
        pmlop = 0.5*(qy(i+1,j+1)+qy(i+1,j))*uy(i,j);
        uy(i,j) = uy(i,j) - dt*oneOverRho*(pmlop + diffop);
        

        end;
    end;
    sfd(a) = px(irx-1,iry-1) + py(irx-1,iry-1);
end;
% tfd = t*1000;
% clf;
% plot(tfd,sfd);
% hold on 
record=sfd;



    
