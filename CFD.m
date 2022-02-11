% READ BEFORE RUNNING CODE:
% This code has been set up so no inputs are needed to run simulation.
% To run code type "CFD()" (without quotation marks) into Command Window and press ENTER.
%
% After running code 3 plots will appear.
% Figure 1: Streamline plot of cavity-driven flow
% Figure 2: Surface plot of the divergence from continuity (VERY small - on the order of 1E-13).
% Figure 3: Surface plot of U-velocity, V-velocity, and Pressure in cavity.
%
% Author: Nicholas Yearout

function[]=CFD()

L = 1;
h = 4.567;
beta = 1.345;
mu  = 1;
rho = 1;
N = 32;
INCR = 1e2;
omega = 4*pi;
IC = 1;
PLOT=1;
BCFLAG=2;
SAVEMOVIE=0;
Ue = 150*ones(N+2,N+1);

dt = 1e-3;
Nt = 1.5e3;

 [SFCN,Output_Real,Output,U1,V1,P1] = Main(N,N,L,L,Nt,dt,rho,mu,beta,h,omega,IC,INCR,PLOT,Ue,SAVEMOVIE);
 
if (SAVEMOVIE==1)
filename='MOVDATA';
save (filename,'SFCN')
end
  

function[SFCN,Output_Real,Output,U1,V1,P1] = Main(Nx,Ny,Lx,Ly,Nt,dt,rho,mu,beta,h,omega,IC,INCR,PLOT,Ue,SAVEMOVIE)

% ----- Build grid ---------
XFLAG = 0;
YFLAG = 0;
PLOTGRID = 0;
[xu,yu,xv,yv,xp,yp] = Grid(Lx,Ly,Nx,Ny,XFLAG,YFLAG,PLOTGRID);

[Xu,Yu] = meshgrid(xu,yu);
[Xv,Yv] = meshgrid(xv,yv);
[Xp,Yp] = meshgrid(xp,yp);

%---- Assign boundary condition coefficients -------

% left BCs
au1 = 1.0*ones(Ny+2,1);
bu1 = 0.0*ones(Ny+2,1);
av1 = 1.0*ones(Ny+1,1);
bv1 = 0.0*ones(Ny+1,1);

% right BCs
au2 = 1.0*ones(Ny+2,1);
bu2 = 0.0*ones(Ny+2,1);
av2 = 1.0*ones(Ny+1,1);
bv2 = 0.0*ones(Ny+1,1);

% lower BCs
au3 = 1.0*ones(1,Nx+1);
bu3 = 0.0*ones(1,Nx+1);
av3 = 1.0*ones(1,Nx+2);
bv3 = 0.0*ones(1,Nx+2);

% upper BCs
au4 = 1.0*ones(1,Nx+1);
bu4 = 0.0*ones(1,Nx+1);
av4 = 1.0*ones(1,Nx+2);
bv4 = 0.0*ones(1,Nx+2);

% ----- Build LHS Matrices ------
Au = BuildAu(Nx,Ny,dt,rho,mu,xu,yu,yv,au1,bu1,au2,bu2,au3,bu3,au4,bu4);
Av = BuildAv(Nx,Ny,dt,rho,mu,xv,yv,xu,av1,bv1,av2,bv2,av3,bv3,av4,bv4);
Ap = BuildAp(Nx,Ny,xp,yp,xu,yv);
Au = decomposition(Au);
Av = decomposition(Av);
Ap = decomposition(Ap);

% ------- Initial conditions ------------
if (IC == 1)
    
 U0 = zeros*Xu;
 V0 = zeros*Xv;
 U1 = zeros*Xu;
 V1 = zeros*Xv;
 P0 = zeros*Xp;
 P1 = zeros*Xp;
 
elseif (IC == 2)
    
 U0 = -cos(beta*Xu-h).*sin(beta*Yu-h);
 V0 = sin(beta*Xv-h).*cos(beta*Yv-h);
 U1 = -cos(beta*Xu-h).*sin(beta*Yu-h);
 V1 = sin(beta*Xv-h).*cos(beta*Yv-h);
 P0 = sin(beta*Xp-h).*sin(beta*Yp-h);
 P1 = sin(beta*Xp-h).*sin(beta*Yp-h);
 
end

[NLU0,NLV0] = NonLin_new(U0,V0,Nx,Ny,rho,xu,yv);

figure

t = 0;

NSX=10;
NSY=10;
ys = linspace(min(yp),max(yp),NSX);
xs = linspace(min(xp),max(xp),NSY);
[XS,YS]= meshgrid(xs,ys);
stfig=figure;

if (SAVEMOVIE==1)
SFCN=[];
end

for k=1:Nt
    
 t = t + dt;
 
 % ---- Boundary Forcing --------
 Ve = 0.0*(sin(beta*Xv-h).*cos(beta*Yv-h)*cos(omega*t));
 Pe = 0.0*(sin(beta*Xp-h).*sin(beta*Yp-h)*cos(omega*t));
 
 gu1 = 0.0*Ue(:,1);
 gu2 = 0.0*Ue(:,Nx+1);
 gu3 = 0.0*Ue(1,:);
 gu4 = Ue(Ny+2,:);
 
 gv1 = Ve(:,1);
 gv2 = Ve(:,Nx+2);
 gv3 = Ve(1,:);
 gv4 = Ve(Ny+1,:);

 % ---- Body Forcing --------
 Bx = 0.0*(cos(beta*Xu-h).*sin(beta*Yu-h)*( rho*omega*sin(omega*t) - 2*beta^2*mu*cos(omega*t)) - ...
 0.5*rho*beta*sin(2.0*beta*Xu-2*h)*cos(omega*t)^2 + ...
 beta*cos(beta*Xu-h).*sin(beta*Yu-h)*cos(omega*t));

 By = 0.0*(sin(beta*Xv-h).*cos(beta*Yv-h)*(-rho*omega*sin(omega*t) + 2.0*mu*beta^2*cos(omega*t) ) -...
 0.5*rho*beta*sin(2*beta*Yv-2*h)*cos(omega*t)^2 + ...
 beta*sin(beta*Xv-h).*cos(beta*Yv-h)*cos(omega*t));

 %------- Build right-hand-side vectors ------
 Fu = RHSu(U0,NLU0,P0,Bx,Nx,Ny,dt,rho,xu,yv,gu1,gu2,gu3,gu4);
 Fv = RHSv(V0,NLV0,P0,By,Nx,Ny,dt,rho,xu,yv,gv1,gv2,gv3,gv4);
 
 %------ Solve matrix problems for U and V -------
 TMPU = Au\Fu;
 TMPV = Av\Fv;
 
 for i=1:Ny+2
 U1(i,:) = TMPU( (i-1)*(Nx+1) + 1 : i*(Nx+1));
 end
 
 for i=1:Ny+1
 V1(i,:) = TMPV( (i-1)*(Nx+2) + 1 : i*(Nx+2));
 end
 
 Fp = RHSp(U1,V1,Nx,Ny,rho,dt,xu,yv);
 
 TMPP = Ap\Fp;
 for i=1:Ny
 PHI(i,:) = TMPP( (i-1)*(Nx) + 1 : i*(Nx));
 end
 
 P0=P1;
 [U1,V1,P1] = Correction(U1,V1,P1,PHI,Nx,Ny,rho,dt,xp,yp,gu1,gu2,gu3,gu4,gv1,gv2,gv3,gv4);
 
[div] = Div(Nx,Ny,U1,V1,xu,yv);

[psi] = Psi(Nx,Ny,U1,V1,xp,yp);

if (SAVEMOVIE==1)
SFCN=[SFCN; psi];
end

Xq=0.5*ones(17,1);
Yq=[1.00; 0.9766; 0.9688; 0.9609; 0.9531; 0.8516; 0.7344; 0.6172; 0.500; 0.4531; 0.2813; 0.1719; 0.1016; 0.0703; 0.0625; 0.0547; 0.00];
Output=interp2(Xu,Yu,U1,Xq,Yq);
Output_Real=Output/100;
 
 if ( mod(k,INCR) == 0 )
     
 k;
 Du = max(max(abs(U1-U0)))/dt
 Dv = max(max(abs(V1-V0)))/dt
 Dp = max(max(abs(P1-P0)))/dt
 
 if (PLOT == 1)
 
 Svals = interp2(Xp,Yp,psi,XS,YS);
 Svals = Svals(:)';
 figure(1)
 hold off
%  contour(Xp,Yp,psi,[Svals],'-k');
 contour(Xp,Yp,psi,'-k');
 title('Streamline Plot')
     
 figure (2)
 surf(Xp,Yp,div), shading interp
 xlabel('x')
 ylabel('y')
 zlabel('Div')
 title('Divergence Plot')
     
 figure(3)
 subplot(1,3,1)
 surf(Xu,Yu,U1), shading interp
 xlabel('x')
 ylabel('y')
 zlabel('Un')
 subplot(2,3,2)

 subplot(1,3,2)
 surf(Xv,Yv,V1), shading interp
 xlabel('x')
 ylabel('y')
 zlabel('Vn')

 subplot(1,3,3)
 surf(Xp,Yp,P1), shading interp
 xlabel('x')
 ylabel('y')
 zlabel('Pn')
 
 pause(0.1)
 
 end
 end
 
 % update the velocity fields and nonlinear terms
 U0 = U1;
 V0 = V1;
 [NLU0,NLV0] = NonLin_new(U0,V0,Nx,Ny,rho,xu,yv);
 
end

Eu = max(max(abs(U1-Ue)));
Ev = max(max(abs(V1-Ve)));
end

function[div] = Div(Nx,Ny,U1,V1,xu,yv)

for j = 1:Nx
for i = 1:Ny
    dx=(xu(j+1)-xu(j));
    dy=(yv(i+1)-yv(i));
    div(i,j)=(U1(i+1,j+1)-U1(i+1,j))/dx + (V1(i+1,j+1)-V1(i,j+1))/dy;
    
    CheckU=(U1(i+1,j+1)-U1(i+1,j))/dx;
    CheckV=(V1(i+1,j+1)-V1(i,j+1))/dy;
end
end

end

function[psi] = Psi(Nx,Ny,U1,V1,xp,yp)

% Initialize PSI matrix
% This is the trapezoidal rule but dy^2 & dx^2 changed to dy & dx
psi=zeros(Ny,Nx);

for i=1:Ny-1
    Us=(U1(i+1,1)+U1(i+1,2))/2;
    Un=(U1(i+2,1)+U1(i+2,2))/2;
    psi(i+1,1)=psi(i,1)+((yp(i+1)-yp(i)))*((Us+Un)/2);
end

for j = 1:Nx-1
for i = 1:Ny
    Vw=(V1(i,j+1)+V1(i+1,j+1))/2;
    Ve=(V1(i,j+2)+V1(i+1,j+2))/2;
    psi(i,j+1)=psi(i,j)-((xp(j+1)-xp(j)))*((Ve+Vw)/2);
end
end

end

function[Ap] = BuildAp(Nx,Ny,xp,yp,xu,yv)


ctr = 0;

% Internal Nodes
for j = 1:Nx
for i = 1:Ny   

  xj = xp(j);
  yi = yp(i);
  dx = xu(j+1) - xu(j);
  dy = yv(i+1) - yv(i);
  
  if (j < Nx)
    xe = xp(j+1);
    ae = 1/(xe-xj)/dx;   
  else
    ae = 0;
  end
  
  if (j > 1)
    xw = xp(j-1);
    aw = 1/(xj-xw)/dx;
  else
    aw = 0;
  end
  
  if (i < Ny)
    yn = yp(i+1);
    an = 1/(yn-yi)/dy;
  else
    an = 0;
  end
  
  if (i > 1)
    ys = yp(i-1);
    as = 1/(yi-ys)/dy;  
  else 
    as = 0;
  end
  
  a =  -(ae + aw  + an + as);  
  
  ctr = ctr+1;
  row(ctr) = (i-1)*Nx + j;
  col(ctr) = (i-1)*Nx + j;
  A(ctr)   = a;
  
  if (j > 1)
    ctr = ctr+1;
    row(ctr) = (i-1)*Nx + j;
    col(ctr) = (i-1)*Nx + j - 1;
    A(ctr)   = aw;
  end  
  
  if (j < Nx)
    ctr = ctr+1;
    row(ctr) = (i-1)*Nx + j;
    col(ctr) = (i-1)*Nx + j + 1;
    A(ctr)   = ae;
  end
  
  if (i < Ny)
    ctr = ctr+1;
    row(ctr) = (i-1)*Nx + j;
    col(ctr) = (i)*Nx   + j;
    A(ctr)   = an;
  end
  
  if (i > 1)
    ctr = ctr+1;
    row(ctr) = (i-1)*Nx + j;
    col(ctr) = (i-2)*Nx + j;
    A(ctr)   = as;
  end  
  
  
end
end

% Add Lagrange Multiplier
for i = 1:Nx*Ny
  
  ctr = ctr+1;
  row(ctr) = i;
  col(ctr) = Nx*Ny + 1;
  A(ctr)   = 1.0;
  
end

ctr = ctr+1;
row(ctr) = Nx*Ny + 1;
col(ctr) = 1;
A(ctr)   = 1.0;

Ap = sparse(row,col,A);

end

function[Au] = BuildAu(Nx,Ny,dt,rho,mu,xu,yu,yv,a1,b1,a2,b2,a3,b3,a4,b4)


ctr=0;

%lower wall boundary conditions
dy = yu(2)-yu(1);
for j=2:Nx
  
  ctr = ctr+1;
  row(ctr) = j;
  col(ctr) = row(ctr);
  A(ctr)   = a3(j) - b3(j)/dy;
  
  ctr = ctr+1;
  row(ctr) = j;
  col(ctr) = row(ctr) + (Nx+1);
  A(ctr)   = b3(j)/dy;
  
end

%upper wall boundary conditions
dy = yu(Ny+2)-yu(Ny+1);
for j=2:Nx
  
  ctr = ctr+1;
  row(ctr) = (Ny+1)*(Nx+1) + j;
  col(ctr) = row(ctr);
  A(ctr)   = a4(j) + b4(j)/dy;
  
  ctr = ctr+1;
  row(ctr) = (Ny+1)*(Nx+1) + j;
  col(ctr) = row(ctr) - (Nx+1);
  A(ctr)   = -b4(j)/dy;
end

% left boundary condition
dx  = xu(2)-xu(1);
for i = 1:Ny+2
  
  ctr = ctr+1;
  row(ctr) = (i-1)*(Nx+1) + 1;
  col(ctr) = row(ctr);
  A(ctr)   = a1(i) - b1(i)/dx;
  
  ctr = ctr+1;
  row(ctr) = (i-1)*(Nx+1) + 1;
  col(ctr) = row(ctr) + 1;
  A(ctr)   = b1(i)/dx;
   
end

% right boundary condition
dx  = xu(Nx+1)-xu(Nx);
for i = 1:Ny+2
  
  ctr = ctr+1;
  row(ctr) = i*(Nx+1);
  col(ctr) = row(ctr);
  A(ctr)   = a2(i) + b2(i)/dx;
  
  ctr = ctr+1;
  row(ctr) = i*(Nx+1);
  col(ctr) = row(ctr) - 1;
  A(ctr)   = -b2(i)/dx;
  
end

% Internal Nodes
for j = 2:Nx
for i = 2:Ny+1   

  dx = 0.5*( xu(j+1) - xu(j-1) );  
  dy = yv(i) - yv(i-1);
  xe = xu(j+1);
  xj = xu(j);
  xw = xu(j-1);
  yn = yu(i+1);
  yi = yu(i);
  ys = yu(i-1);
  
  aw = -mu*dy/(xj-xw);
  ae = -mu*dy/(xe-xj);
  as = -mu*dx/(yi-ys);
  an = -mu*dx/(yn-yi);    
  ap =  rho*dx*dy/dt - (aw + ae + as + an);
  
  ctr = ctr+1;
  row(ctr) = (i-1)*(Nx+1)+j;
  col(ctr) = (i-1)*(Nx+1)+j;
  A(ctr)   = ap;
  
  ctr = ctr+1;
  row(ctr) = (i-1)*(Nx+1)+j;
  col(ctr) = (i-1)*(Nx+1)+j-1;
  A(ctr)   = aw;
  
  ctr = ctr+1;
  row(ctr) = (i-1)*(Nx+1)+j;
  col(ctr) = (i-1)*(Nx+1)+j+1;
  A(ctr)   = ae;
  
  ctr = ctr+1;
  row(ctr) = (i-1)*(Nx+1)+j;
  col(ctr) = (i)*(Nx+1)+j;
  A(ctr)   = an;
  
  ctr = ctr+1;
  row(ctr) = (i-1)*(Nx+1)+j;
  col(ctr) = (i-2)*(Nx+1)+j;
  A(ctr)   = as;
    
      
end
end

Au = sparse(row,col,A);

end

function[Av] = BuildAv(Nx,Ny,dt,rho,mu,xv,yv,xu,a1,b1,a2,b2,a3,b3,a4,b4)


ctr=0;

%lower wall boundary conditions
dy = yv(2)-yv(1);
for j=2:Nx+1
  
  ctr = ctr+1;
  row(ctr) = j;
  col(ctr) = j;
  A(ctr)   = a3(j)-b3(j)/dy;

  ctr = ctr+1;
  row(ctr) = j;
  col(ctr) = j + (Nx+2);
  A(ctr)   = b3(j)/dy;
  
end

%upper wall boundary conditions
dy = yv(Ny+1)-yv(Ny);
for j=2:Nx+1
  
  ctr = ctr+1;
  row(ctr) = Ny*(Nx+2) + j;
  col(ctr) = row(ctr);
  A(ctr)   = a4(j) + b4(j)/dy;
  
  ctr = ctr+1;
  row(ctr) = Ny*(Nx+2) + j;
  col(ctr) = row(ctr) - (Nx+2);
  A(ctr)   = -b4(j)/dy;
  
  
end

% left boundary condition
dx = xv(2)-xv(1);
for i = 1:Ny+1
  
  ctr = ctr+1;
  row(ctr) = (i-1)*(Nx+2) + 1;
  col(ctr) = row(ctr);
  A(ctr)   = a1(i) - b1(i)/dx;
  
  ctr = ctr+1;
  row(ctr) = (i-1)*(Nx+2) + 1;
  col(ctr) = row(ctr)+1;
  A(ctr)   = b1(i)/dx;
  
end

% right boundary condition
dx = xv(Nx+2) - xv(Nx+1);
for i = 1:Ny+1
  
  ctr = ctr+1;
  row(ctr) = i*(Nx+2);
  col(ctr) = row(ctr);
  A(ctr)   = a2(i) + b2(i)/dx;
  
  ctr = ctr+1;
  row(ctr) = i*(Nx+2);
  col(ctr) = row(ctr) -1;
  A(ctr)   = -b2(i)/dx;
  
end

% Internal Nodes
for j = 2:Nx+1
for i = 2:Ny   

  dx = xu(j) - xu(j-1);  
  dy = 0.5*( yv(i+1) - yv(i-1) );
  xe = xv(j+1);
  xj = xv(j);
  xw = xv(j-1);
  yn = yv(i+1);
  yi = yv(i);
  ys = yv(i-1);
  
  aw = -mu*dy/(xj-xw);
  ae = -mu*dy/(xe-xj);
  as = -mu*dx/(yi-ys);
  an = -mu*dx/(yn-yi);  
  ap = rho*dx*dy/dt - (aw + ae + as + an);
  
  ctr = ctr+1;
  row(ctr) = (i-1)*(Nx+2)+j;
  col(ctr) = (i-1)*(Nx+2)+j;
  A(ctr)   = ap;
  
  ctr = ctr+1;
  row(ctr) = (i-1)*(Nx+2)+j;
  col(ctr) = (i-1)*(Nx+2)+j-1;
  A(ctr)   = aw;
  
  ctr = ctr+1;
  row(ctr) = (i-1)*(Nx+2)+j;
  col(ctr) = (i-1)*(Nx+2)+j+1;
  A(ctr)   = ae;
  
  ctr = ctr+1;
  row(ctr) = (i-1)*(Nx+2)+j;
  col(ctr) = (i)*(Nx+2)+j;
  A(ctr)   = an;
  
  ctr = ctr+1;
  row(ctr) = (i-1)*(Nx+2)+j;
  col(ctr) = (i-2)*(Nx+2)+j;
  A(ctr)   = as;
    
end
end

Av = sparse(row,col,A);

end

function[U1,V1,P1] = Correction(U1,V1,P1,PHI,Nx,Ny,rho,dt,xp,yp,...
                     gu1,gu2,gu3,gu4,gv1,gv2,gv3,gv4)

  % update Pressure
  P1 = P1 + PHI;
  
  
  % update U
  for i = 2:Ny+1
  for j = 2:Nx
    
    U1(i,j) = U1(i,j) - dt/rho*( PHI(i-1,j) - PHI(i-1,j-1) )/( xp(j)-xp(j-1) );
    
  end
  end
  
  U1(:,1)    = gu1;
  U1(:,Nx+1) = gu2;
  U1(1,:)    = gu3;
  U1(Ny+2,:) = gu4;
 
  
  % update V
  for i = 2:Ny
  for j = 2:Nx+1
    
    V1(i,j) = V1(i,j) - dt/rho*( PHI(i,j-1) - PHI(i-1,j-1) )/( yp(i)-yp(i-1) );
    
  end
  end
  
  V1(:,1)    = gv1;
  V1(:,Nx+2) = gv2;
  V1(1,:)    = gv3;
  V1(Ny+1,:) = gv4;
  
end

function[xu,yu,xv,yv,xp,yp] = Grid(Lx,Ly,Nx,Ny,XFLAG,YFLAG,PLOTGRID)

% This routine builds a staggered grid of Nx by Ny finite volume cells
% in a rectangular domain of size 0<x<Lx and 0<y<Ly. 
% Grid is uniform in x-direction if XFLAG = 0 and non-uniform if XFLAG = 1
% Grid is uniform in y-direction if YFLAG = 0 and non-uniform if YFLAG = 1


% Initialize the vectors
xp = zeros(Nx,1);
xu = zeros(Nx+1,1);
xv = zeros(Nx+2,1);

yp = zeros(Ny,1);
yu = zeros(Ny+2,1);
yv = zeros(Ny+1,1);


% ----- compute x-locations -----

if ( XFLAG == 0 ) % uniform
  
  dx = Lx/Nx;
  for i=1:Nx+1;
    xu(i) = (i-1)*dx;
  end
  
elseif ( XFLAG == 1 ) % non-uniform

  xu = 0.5*Lx*( 1 + cos(pi*[Nx:-1:0]'/Nx) );
  
end

xv(1) = 0;
xv(Nx+2) = Lx;
for i = 2:Nx+1
  xv(i) = 0.5*( xu(i-1) + xu(i) );
  xp(i-1) = xv(i);
end


%----- compute y-locations ---------

if ( YFLAG == 0 ) % uniform
  
  dy = Ly/Ny;
  for i=1:Ny+1;
    yv(i) = (i-1)*dy;
  end

elseif ( YFLAG == 1 ) % non-uniform

  yv = 0.5*Ly*( 1 + cos(pi*[Ny:-1:0]'/Ny) );
  
end

yu(1) = 0;
yu(Ny+2) = Ly;
for i = 2:Ny+1
  yu(i) = 0.5*( yv(i-1) + yv(i) );
  yp(i-1) = yu(i); 
end

% Plot the grid
if ( PLOTGRID == 1 )
  
  figure
  hold on
  
  % Plot horizontal lines
  for i=1:Ny+1
    plot([0 Lx],[yv(i) yv(i)],'-k','LineWidth',2.0)  
  end
  
  % Plot vertical lines
  for i=1:Nx+1
    plot([xu(i) xu(i)],[0 Ly],'-k','LineWidth',2.0)
  end
  
  for i=1:Ny
  for j=1:Nx
    plot(xp(j),yp(i),'.k','MarkerSize',24)
  end
  end
  
  for i=1:Ny+2
  for j=1:Nx+1
    plot(xu(j),yu(i),'sk','MarkerSize',14,'LineWidth',2)
  end
  end
  
  for i=1:Ny+1
  for j=1:Nx+2
    plot(xv(j),yv(i),'^k','MarkerSize',10,'LineWidth',2)
  end
  end
  
end

end

function[] = Movie(NSX,NSY)

% inputs: NSX and NSY are the number of streamlines in the x and y directions. 
%         NSX=NSY=10, for a total of 100 lines.
load MOVDATA_psi
Ny=100;
Nx=100;
Lx = 1;
Ly = 1;
XFLAG = 0;
YFLAG = 0;
PLOTGRID = 0;
NSX=10;
NSY=10;

[xu,yu,xv,yv,xp,yp] = Grid(Lx,Ly,Nx,Ny,XFLAG,YFLAG,PLOTGRID);
[Xp,Yp] = meshgrid(xp,yp);

v = VideoWriter('Movie.avi');
v.FrameRate = 5;
v.Quality = 98;
open(v);

[rows,cols]=size(SFCN);
k=rows/(Ny)

ys = linspace(min(yp),max(yp),NSX);
xs = linspace(min(xp),max(xp),NSY);
[XS,YS]= meshgrid(xs,ys);

  
mygraph=figure
set(mygraph, 'PaperPositionMode','auto')  
set(gca,'FontSize',16,'FontName','Times')


k=10;
for i=1:k
  
  sf = SFCN( (i-1)*(Ny)+1:(i)*(Ny),:);
  
  hold off
  Svals = interp2(Xp,Yp,sf,XS,YS);
  Svals = Svals(:)';
  contour(Xp,Yp,sf,[Svals],'-k');
  hold on
   
  pause(.1)
  frame = getframe(gcf);
  writeVideo(v,frame);
end

end

function[NLU,NLV] = NonLin_new(U,V,Nx,Ny,rho,xu,yv)

NLU = zeros(Ny+2,Nx+1);
NLV = zeros(Ny+1,Nx+2);

% Compute NLU
for i = 2:Ny+1
for j=2:Nx
      
    dy = yv(i) - yv(i-1);
    dxe = xu(j+1) - xu(j);
    dxw = xu(j)   - xu(j-1);
    
    % Compute the mass fluxes so that mass is conserved on the 
    % x-momentum cell
    me = 0.5*rho*dy*( U(i,j)   + U(i,j+1)   );
    mw = 0.5*rho*dy*( U(i,j)   + U(i,j-1)   );
    mn = 0.5*rho*( V(i,j)*dxw   + V(i,j+1)*dxe );
    ms = 0.5*rho*( V(i-1,j)*dxw + V(i-1,j+1)*dxe );
    
    %---- Centered Differences-----------%
    ue = 0.5*( U(i,j) + U(i,j+1) );
    uw = 0.5*( U(i,j) + U(i,j-1) );
    
    if (i==Ny+1) % account for upper boundary where ctrd diff not necessary
      un = U(i+1,j);
    else
      un = 0.5*( U(i,j) + U(i+1,j) );
    end
    if (i==2) % account for lower boundary where ctrd diff not necessary
      us = U(i-1,j);
    else
      us = 0.5*( U(i,j) + U(i-1,j) );
    end  
    NLU(i,j) = me*ue + mn*un - mw*uw - ms*us;
    
end  
end

% Compute NLV
for i = 2:Ny
for j=2:Nx+1
        
    dx = xu(j) - xu(j-1);  
    dyn = yv(i+1) - yv(i);
    dys = yv(i)   - yv(i-1);
    
    
    %--------- Mass Fluxes for the Y-Momentum Cell ----------%
    me = 0.5*rho*( U(i,j)*dys   + U(i+1,j)*dyn );
    mw = 0.5*rho*( U(i,j-1)*dys   + U(i+1,j-1)*dyn );
    mn = 0.5*rho*dx*( V(i,j) +V(i+1,j) );   
    ms = 0.5*rho*dx*( V(i,j) +V(i-1,j) );   
    
    % --------- Centered Differences ----------%
    if (j==Nx+1) % account for right boundary where ctrd diff not necessary
      ve = V(i,j+1);
    else
      ve = 0.5*( V(i,j) + V(i,j+1) );
    end
    if (j==2) % account for left boundary where ctrd diff not necessary
      vw = V(i,j-1);
    else
      vw = 0.5*( V(i,j) + V(i,j-1) );
    end
    vn = 0.5*( V(i,j) + V(i+1,j) );
    vs = 0.5*( V(i,j) + V(i-1,j) );
    
    NLV(i,j) = me*ve + mn*vn - mw*vw - ms*vs;

    
end  
end

end

function[Fp] = RHSp(U1,V1,Nx,Ny,rho,dt,xu,yv)

Fp = zeros(Ny*Nx+1,1);

for i = 1:Ny
for j = 1:Nx
  
    dx = xu(j+1) - xu(j);  
    dy = yv(i+1) - yv(i);
   
    Fp( (i-1)*Nx + j ) = rho/dt*( ( U1(i+1,j+1) - U1(i+1,j) )/dx +  ( V1(i+1,j+1) - V1(i,j+1) )/dy );
    
end  
end

Fp(Ny*Nx+1) = 0;

end

function[Fu] = RHSu(U0,NLU,P,Bx,Nx,Ny,dt,rho,xu,yv,g1,g2,g3,g4)

Fu = zeros( (Ny+2)*(Nx+1),1 );

% lower/upper walls
for j = 2:Nx
  Fu(j) = g3(j); % Lower wall
  Fu((Nx+1)*(Ny+1)+j) = g4(j); % upper wall
end

% left/right boundaries
for i = 1:Ny+2
  Fu( (i-1)*(Nx+1)+1 ) = g1(i);   % left boundary
  Fu( i*(Nx+1) ) = g2(i);        %right boundary
end

% interior points
for i = 2:Ny+1
for j = 2:Nx
  
  dx = 0.5*( xu(j+1) - xu(j-1) );
  dy = yv(i) - yv(i-1);  
  
  % Add pressure term in line below.
  Fu( (i-1)*(Nx+1) + j ) = rho*U0(i,j)*dx*dy/dt - NLU(i,j) + Bx(i,j)*dx*dy - dy*(P(i-1,j)-P(i-1,j-1));
  
end  
end

end

function[Fv] = RHSv(V0,NLV,P,By,Nx,Ny,dt,rho,xu,yv,g1,g2,g3,g4)

Fv = zeros( (Ny+1)*(Nx+2),1 );

% lower/upper walls
for j = 2:Nx+1  
  Fv(j) = g3(j); % lower wall
  Fv((Nx+2)*Ny + j) = g4(j); % upper wall
end

% left/right boundaries
for i = 1:Ny+1
  Fv( (i-1)*(Nx+2)+1 ) = g1(i);   % left boundary
  Fv( i*(Nx+2) ) = g2(i);   %right boundary
end

% interior points
for i = 2:Ny
for j = 2:Nx+1
  
    dx = xu(j) - xu(j-1);  
    dy = 0.5*( yv(i+1) - yv(i-1) );
    
    % Add pressure term in line below.
    Fv( (i-1)*(Nx+2) + j ) = rho*V0(i,j)*dx*dy/dt - NLV(i,j) + By(i,j)*dy*dx - dx*(P(i,j-1)-P(i-1,j-1));
  
end
end

end

end