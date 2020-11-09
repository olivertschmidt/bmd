% EXAMPLE 2 demonstrates the manual specification of spectral estimation
% parameters and regions of the mode bispectrum, and its computation for
% data consisting of more than one variable (here the two velocity
% components of the wake data). The modes are plotted interactively by
% clicking on the mode bispectrum. Try setting opts.regions to 1:8 to
% compute all regions, or [1 2 7 8] to compute regions that are
% non-redundent with respect to complex conjugation.

clc, clear variables

%% Load data
load('wake_Re500.mat')
[nt,nx,ny]  = size(u);

nDFT        = 256;
nOvlp       = 128;
window      = hann(nDFT);
dV          = (x(2,1)-x(1,1))*(y(1,2)-y(1,1));
weight      = dV*ones(2*ny*nx,1);
opts.regions= 1; 
X           = zeros(nt,nx,ny,2);
X(:,:,:,1)  = u;
X(:,:,:,2)  = v; clear u v

%%%%%%%%%
%% BMD %%
%%%%%%%%%
[B,P,f,idx] = bmd(X,nDFT,weight,nOvlp,dt,opts);
[f1,f2]     = ndgrid(f);

%% Mode magnitude bispectrum
spec_fig = figure;
subplot(3,4,[1 2 5 6 9 10]);
contourf(f1,f2,log(abs(B)),100,'linecolor','none'); axis equal tight
c   = colorbar('Location','east'); c.Label.String = '|\lambda_1|';
xlabel('f_1'), ylabel('f_2'), zlabel('|\lambda_1|'); axis equal
xlim([min(f1(idx)) max(f1(idx))]); ylim([min(f2(idx)) max(f2(idx))])
title('Mode bispectrum')

%% Interactive visualization of modes
dcm             = datacursormode(spec_fig);
dcm.Enable      = 'on';
dcm.UpdateFcn   = @displayTriplet;
disp(' ');
disp('Click any point in the mode bispecrtum to plot modes, or ESC to exit visualization mode!');
while 1
    waitforbuttonpress;
    
    % exit while loop if ESC is pressed
    key = get(gcf,'currentcharacter');
    if key==27, break, end
       
    point       = getCursorInfo(dcm);
    triadIdx    = find(idx==find(f1==point.Position(1)&f2==point.Position(2)));
    
    subplot(3,4,3)
    mode  = squeeze(P(1,triadIdx,:,:,1));
    pcolor(x,y,real(mode)), axis equal tight, caxis(max(abs(mode(:)))*0.5*[-1 1])
    xlabel('x'), ylabel('y'), title('\phi^u_{k+l}')
    shading interp
    
    subplot(3,4,4)
    mode  = squeeze(P(1,triadIdx,:,:,2));
    pcolor(x,y,real(mode)), axis equal tight, caxis(max(abs(mode(:)))*0.5*[-1 1])
    xlabel('x'), ylabel('y'), title('\phi^v_{k+l}');
    shading interp
    
    subplot(3,4,7)
    mode  = squeeze(P(2,triadIdx,:,:,1));
    pcolor(x,y,real(mode)), axis equal tight, caxis(max(abs(mode(:)))*0.5*[-1 1])
    xlabel('x'), ylabel('y'), title('\phi^u_{k\circ l}');
    shading interp
    
    subplot(3,4,8)
    mode  = squeeze(P(2,triadIdx,:,:,2));
    pcolor(x,y,real(mode)), axis equal tight, caxis(max(abs(mode(:)))*0.5*[-1 1])
    xlabel('x'), ylabel('y'), title('\phi^v_{k\circ l}');
    shading interp   
    
    subplot(3,4,11)
    mode  = abs(squeeze(P(1,triadIdx,:,:,1)).*squeeze(P(2,triadIdx,:,:,1)));
    pcolor(x,y,real(mode)), axis equal tight
    xlabel('x'), ylabel('y'), title('\phi^u_{k\circ l}\circ\phi^u_{k+l}');
    shading interp
    
    subplot(3,4,12)
    mode  = abs(squeeze(P(1,triadIdx,:,:,2)).*squeeze(P(2,triadIdx,:,:,2)));
    pcolor(x,y,real(mode)), axis equal tight
    xlabel('x'), ylabel('y'), title('\phi^v_{k\circ l}\circ\phi^v_{+}');
    shading interp     
    
    drawnow   
    figure(spec_fig)
end

function txt = displayTriplet(~,info)
    x = info.Position(1);
    y = info.Position(2);
    txt = ['{' num2str(x,'%.2g') ',' num2str(y,'%.2g') ',' num2str(x+y,'%.2g') '}'];
end