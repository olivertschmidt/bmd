% EXAMPLE 3 demonstrates the use of CBMD (Cross-Bispectreal Mode
% Decomposition). Solely for the purpose of demonstration, we compute the
% CBMD in the special case where all three variables are the same, in which
% CBMD reduces to BMD. When computing the CBMD between different variables,
% consider computing more regions of the mode bispectrum as symetries are
% lost. All symmetries are lost if all three variables are both different
% and complex.

clc, clear variables

%% Load data
load('wake_Re500.mat')
[nt,nx,ny]  = size(u);

%% Assemble compound data array X containing all 3 variables in second index
X           = zeros(nt,3,nx,ny);
% Note: all three variables are set to 'u', and we should therefore recover
% the results from example 1.

% Q
X(:,1,:,:)  = u;
% R
X(:,2,:,:)  = u;
% S
X(:,3,:,:)  = u;
clear u v

%%%%%%%%%%
%% CBMD %%
%%%%%%%%%%
[B,P,f,idx] = cbmd(X);
[f1,f2]     = ndgrid(f);

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
triadIdx    = find(idx==find(f1==k&f2==l));
plot(k,l,'r+'), text(k+1,l+1,['(' num2str(k) ',' num2str(l) ',' num2str(k+l) ')'])

subplot(2,2,2)
mode  = squeeze(P(1,triadIdx,:,:));
pcolor(x,y,real(mode)), axis equal tight, caxis(max(abs(mode(:)))*0.5*[-1 1]); shading interp
xlabel('x'), ylabel('y'), title(['\phi_{' num2str(k) num2str(l,'%+d') '} (bispectral mode)']);

subplot(2,2,4)
psi  = abs(squeeze(P(1,triadIdx,:,:).*P(2,triadIdx,:,:)));
pcolor(x,y,psi), axis equal tight, caxis(max(abs(psi(:)))*[0 1]); shading interp
xlabel('x'), ylabel('y'), title(['\psi_{' num2str(k) num2str(l,'%+d') '} (interaction map)']);

