
% EXAMPLE 1 demonstrates the use of BMD (Bispectral Mode Decomposition) on
% a single-variable 2-D flow field (the streamwise velocity component of
% the wake behind a cylinder at Re=500). No spectral estimation parameters
% are specified, so it is left to BMD to guess the parameters. This is not
% recommended, in general. Since the time step is not specified, BMD
% returns the frequency index in 'f'. BMD assumes that the first index of
% the data array is time.

clc, clear variables

%% Load data
load('wake_Re500.mat')
[nt,nx,ny]  = size(u);

%% Inspect wake data
figure
u_mean  = squeeze(mean(u,1)); v_mean  = squeeze(mean(v,1));
for ti=1:100
    subplot(2,1,1)
    pcolor(x,y,squeeze(u(ti,:,:))-u_mean); shading interp; axis equal tight
    xlabel('x'), ylabel('y'); title('u'''); caxis([-0.5 0.5]);
    subplot(2,1,2)
    pcolor(x,y,squeeze(v(ti,:,:))-v_mean); shading interp; axis equal tight
    xlabel('x'), ylabel('y'); title('v'''); caxis([-1 1]);
    drawnow
end

%%%%%%%%%
%% BMD %%
%%%%%%%%%
[B,P,f,idx] = bmd(u);
[f1,f2]     = ndgrid(f); % don't use meshgrid...

%% Plot mode magnitude bispectrum
figure
subplot(2,2,[1 3])
pcolor(f1-0.5,f2-0.5,log(abs(B))); axis equal tight, shading flat
xlim([0 f1(end)]); ylim([-f2(end) f2(end)/2]); 
xlabel('k'), ylabel('l')
title('Mode bispectrum')
hold on

%% Plot bispectral mode and interaction map for triad (k,l,k+l)
k   = 6;
l   = -3;
triadIdx    = find(idx==find(f1==k&f2==l)); % find the linear index of the triad (k,l,k+l)
plot(k,l,'r+'), text(k+1,l+1,['(' num2str(k) ',' num2str(l) ',' num2str(k+l) ')'])

subplot(2,2,2)
mode  = squeeze(P(1,triadIdx,:,:));
pcolor(x,y,real(mode)), axis equal tight, caxis(max(abs(mode(:)))*0.5*[-1 1]); shading interp
xlabel('x'), ylabel('y'), title(['\phi_{' num2str(k) num2str(l,'%+d') '} (bispectral mode)']);

subplot(2,2,4)
psi  = abs(squeeze(P(1,triadIdx,:,:).*P(2,triadIdx,:,:)));
pcolor(x,y,psi), axis equal tight, caxis(max(abs(psi(:)))*[0 1]); shading interp
xlabel('x'), ylabel('y'), title(['\psi_{' num2str(k) num2str(l,'%+d') '} (interaction map)']);

