clear, close all
set(0,'defaultAxesFontSize',14), set(0,'DefaultLineLineWidth',1.5);
set(groot,'defaultAxesXGrid','on'), set(groot,'defaultAxesYGrid','on')
p_g=9.81; p_f=1.3e-4; p_N=1e-3; p_R0=1027;
iPlot=1; iSave=1;

dx=100; dy=100; dz=40;
Nx=240; Ny=200; Nz=50;
xc=dx/2:dx:Nx*dx-dx/2; yc=dy/2:dy:Ny*dy-dy/2;
Depth=-Nz*dz; 
zc=-dz/2:-dz:Depth+dz/2; zf=0:-dz:Depth;

[gridXY,gridYX]=meshgrid(xc,yc); gridXY=gridXY'; gridYX=gridYX';
[gridYZ,gridZY]=meshgrid(yc,zc); gridYZ=gridYZ'; gridZY=gridZY';

%% Topography
amp=5; wavelen=1200; k0=2*pi/wavelen; Umax0=0.18; %Umax=0.09; Umax=0.06; %Umax=0.03;
% fprintf('Intrinsic freq at generation = %.2f N = %.2f f \n',k0*Umax/p_N,k0*Umax/p_f)
H=Depth+amp*cos(k0*(gridXY-dx/2))+amp;
% H(1:12,:)=Depth+amp; H(12+12*5+2:end,:)=Depth+amp; 
% H([1:12*2  12*10+1:end],:)=Depth+amp;
H(:,[1 end])=0; % meridional walls
dH=2000+H(:,Ny/2); 

if iPlot
    figure, mesh(gridXY/1e3,gridYX/1e3,H)
    if amp~=0, zlim([Depth Depth+2*amp]), end
    xlabel('x (km)'), ylabel('y (km)'), zlabel('z (m)')
    title('Topographic depth')
end
fprintf('Amplitude = %d m, wavelength = %d m \n',amp,wavelen)
fprintf('Maximum bathy = %d m, average bathy = %.2f m \n',min(H,[],'all'), mean(H(:,2:end-1),'all'))
fprintf('Minimum partial cell thickness = %.2f m \n',min(dH(dH>0),[],'all'))
clear amp wavelen

%% Zonal velocity
zt=Depth+dz; zb=-100;
fprintf('Jet''s bottom at %d m and top at %d m \n', zt,zb)
factor=2*Umax0/(zt-zb)^2;
Uc=factor/2*(zc-zb).^2 .*(zc<zb);
if iPlot
    figure, plot(Uc,zc), ylim([zt 0])
    xlabel('U (m/s)'), ylabel('z (m)')
    title('Jet''s central velocity')
end
clear Depth factor Uc

%% Density
ym=(yc(1)+yc(end))/2; yf=15e3;
fprintf('Jet''s meridional axis at %.2f km \n',ym/1e3)
sigma=pi/(yf-ym);
thy=tanh(sigma*(gridYX-ym));
if 0
    figure, plot(yc,thy(1,:))
end

dUdzmax=2*Umax0/(zt-zb);
dR=p_f*dUdzmax*p_R0/(sigma*p_g);
dRdz=-p_N^2*p_R0/p_g;

Rref=ones(Nz,1)*p_R0;
R=ones(Nx,Ny,Nz)*p_R0;
for iz=1:Nz-1 % for each layer
    z=zc(iz);
    Rref(iz+1)=Rref(iz)-dRdz*dz;
    R(:,:,iz+1)=(z-zb)/(zt-zb)*dR*thy*(z<zb) + Rref(iz+1);
end
if iPlot
    figure, contour(gridYZ,gridZY,squeeze(R(1,:,:))), colorbar
    xlabel('y (km)'), ylabel('z (m)'), title('Density')
end
clear iz z Rref thy dUdzmax dRdz Umax

%% Temperature from density
alpha=2e-4;
T=(1-R/p_R0)/alpha +20;
if iPlot
    figure, contour(gridYZ/1e3,gridZY,squeeze(T(1,:,:))), colorbar
    xlabel('y (km)'), ylabel('z (m)'), title('Temperature')
end
clear alpha R

%% Reconstruct velocity from density
dRdy=(gridZY-zb)/(zt-zb).*(gridZY<zb) *dR .* cosh(sigma*(gridYZ-ym)).^(-2)*sigma;
dUdz=p_g*dRdy/(p_f*p_R0);

U_yz=zeros(Ny,Nz);    
for iz=1:Nz-1 % for each layer
    U_yz(:,iz+1)=U_yz(:,iz)-dUdz(:,iz)*dz;
end
if iPlot
    pcolor(gridYZ/1e3,gridZY,U_yz), shading interp, hold on
    contour(yc/1e3,zc,U_yz',[1e-10 .01:.01:.2],'k')
    colorbar, %colormap('bluewhitered')
    xlabel('y (km)'), ylabel('z (m)'), title('Jet''s zonal velocity')
end

U=zeros(Nx,Ny,Nz); V=zeros(Nx,Ny,Nz);
for ix=1:Nx, U(ix,:,:)=U_yz; end % repmat
Umax=max(U,[],'all');

fprintf('Jet''s max velocity = %.4f m/s \n',Umax)
fprintf('Intrinsic freq at generation = %.2f N = %.2f f \n',k0*Umax/p_N,k0*Umax/p_f)
clear ix iz dR dRdy dUdz sigma U_yz ym yf zb zt k0

%% Viscosity
Ah_min=0.05; Ah=ones(Nx,Ny,Nz)*Ah_min; Ns=10; % sponge layer thickness
for i=1:Ns
    Ah(:,[Ns+1-i Ny-Ns+i],:)=Ah_min*10*i; % amplified 10 times
    % Ah([Ns+1-i Nx-Ns+i],:,:)=Ah_min*10*i;
end
fprintf('Ah is between %.2e and %.2e \n',min(Ah,[],'all'),max(Ah,[],'all'))

if iPlot
    Ah_xy=squeeze(mean(Ah,3));
    figure, pcolor(xc/1e3,yc/1e3,Ah_xy'), shading flat, colorbar
    xlabel('x (km)'), ylabel('y (km)')
    % title('Biharmonic horizontal viscosity (m^4/s)')
    title('Laplacian horizontal viscosity (m^2/s)')
end
clear i Ah_xy

%% Save
if iSave
    % fid=fopen('BAT.bin'   ,'w','b'); fwrite(fid,H,'float32'); fclose(fid);

    fid=fopen('U.bin'     ,'w','b'); fwrite(fid,U,'float32'); fclose(fid);
    fid=fopen('V.bin'     ,'w','b'); fwrite(fid,V,'float32'); fclose(fid);
    fid=fopen('TEMP.bin'  ,'w','b'); fwrite(fid,T,'float32'); fclose(fid);

    fid=fopen('U_B.bin'   ,'w','b'); fwrite(fid,squeeze(U(1,:,:)),'float32'); fclose(fid);
    fid=fopen('V_B.bin'   ,'w','b'); fwrite(fid,squeeze(V(1,:,:)),'float32'); fclose(fid);
    fid=fopen('TEMP_B.bin','w','b'); fwrite(fid,squeeze(T(1,:,:)),'float32'); fclose(fid);
    
    fid=fopen('ViscAh.bin','w','b'); fwrite(fid,Ah,'float32'); fclose(fid);
end
clear fid iPlot iSave
clear gridXY gridYX gridYZ gridZY