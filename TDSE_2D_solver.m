%% Doruk Efe Gökmen -- 28/11/2017 -- 2D Schrödinger Equation Solver v0
clear;clf;clc;
% The time evolution of a given initial state is calculated by iteratively
% acting the time evolution operator on the state vector.
% The time evolution matrix is assumed to decomposable into kinetic energy 
% (K) and potential energy (V) parts with small error, given that the time 
% step is sufficiently small, although [K,V]=/=0.
T = 6660; %time length
N = 128; %number of spatial elements
hbar = 1;
m = 1; %mass of the particle
dt = 50e-2; %temporal integration step

% Position space
dx = 1; dy = 1;  %spatial step
x=-(N-1)/2:dx:(N-1)/2; y=-(N-1)/2:dy:(N-1)/2;
[X,Y]=meshgrid(x,y);

%%Reciprocal space
% kx = 2*pi/N * (0:N-1); ky = 2*pi/N * (0:N-1);
% [Kx,Ky]=meshgrid(kx,ky);

% % Initial state (old method)
% px = @(kx) 2*sin(kx*pi/N);  %initial x momentum (arbitrarily chosen)
% py = @(ky) 2*sin(ky*pi/N); %initial y momentum (arbitrarily chosen)
% TP = exp(-1i * (px(0)*X + py(5)*Y))'; %momentum translation operator
% psi = mvnpdf([X(:) Y(:)],[0 0],[200 1; 1 200]); %Gaussian
% psi = reshape(psi,length(y),length(x));

psi = zeros(N,N); %initialise
A = [1,1];
mu = [40 , 0 ; 0 , 40]; %location of the center of the two Gaussians
sigma = [90 , 90]; %width of two Gaussians
p = [0 , 0 ; 0 , 0];%momentum of the initial state
ppsi = @(x,y)  ...
    A(1) * exp(-1i * (p(1,1)*x + p(1,2)*y)... %momentum 
    -((x-mu(1,1))^2 + (y-mu(1,2))^2)/sigma(1))... %Gaussian
    +A(2) * exp(-1i * (p(2,1)*x + p(2,2)*y)... %symmetric part 
    -((x-mu(2,1))^2 + (y-mu(2,2))^2)/sigma(2));
for i=1:N
    for j=1:N
        psi(i,j)=ppsi(-N/2+i-1,-N/2+j-1)*10e-4;
    end
end
psi=psi/sqrt((sum(sum(conj(psi).*psi)))); %normalisation

% psi1 = mvnpdf([X(:) Y(:)],[30 0],[50 1; 1 250]);
% psi1 = reshape(psi1,length(y),length(x));
% psi2 = mvnpdf([X(:) Y(:)],[0 30],[50 1; 1 250]);
% psi2 = reshape(psi2,length(y),length(x));
% psi=psi1+psi2; %two Gaussian blobs

% %plot initial state
% figure (1)
% contour(-x,y,psi);



% Potential energy matrix in position space
V = zeros(N,N); %free particle (or "null potential")

% % Harmonic well potential
% k0=10^-4;
% v = @(xx,yy) 0.5*k0*(xx^2+yy^2);
% 
% for i=1:N
%     for j=1:N
%         V(i,j)=v(-N/2+i-1,-N/2+j-1);
%     end
% end

% Harmonic well with interaction barrier
V0=10; %interaction strength
k0=10^-4; %spring constant
v = @(xx,yy) 0.5*k0*(xx^2+yy^2);
for i=1:N
    for j=1:N
        V(i,j)=v(-N/2+i-1,-N/2+j-1);
        if abs(i-j)<5
            V(i,j)=v(-N/2+i-1,-N/2+j-1)+V0; %delta barrier
        end
    end
end

% % Square well potential
% for i=1:N
%     for j=1:N
%         if i<35 || j<35
%             V(i,j)=V0;
%         elseif i>N-35 || j>N-35
%             V(i,j)=V0;
%         else
%             V(i,j)=0;
%         end
%     end
% end

% % Square well potential with delta interaction
% for i=1:N
%     for j=1:N
%         if i<35 || j<35
%             V(i,j)=V0;
%         elseif i>N-35 || j>N-35
%             V(i,j)=V0;
%         elseif i==j
%             V(i,j)=V0; %delta barrier
%         else
%             V(i,j)=0;
%         end
%     end
% end

%plot the potential energy
% surf(x,y,V);

% Kinetic energy matrix in reciprocal space
cKE_ = @(kx,ky) hbar^2 * 4*(sin(kx*pi/N)^2+sin(ky*pi/N)^2); 
KE_=zeros(N,N);
for i=1:N
    for j=1:N
        KE_(i,j)=cKE_(i-1,j-1);
    end
end

% Time evolution operator
UV = exp(-1i*dt/hbar *V); %potential energy part (diagonal in position space)
UK_ = exp(-1i*dt/hbar *KE_); %kinetic energy part (diagonal in reciprocal space)

avg_x=zeros(1,T); %initialise the average values
avg_y=avg_x;

XgY=fliplr(tril(fliplr(X))); %"x for x>y" operator
YgX=fliplr(tril(fliplr(Y))); %"y for y>x" operator

XgY=triu(X);

figure(2)
for i = 1:T
    set(gcf, 'Position', get(0, 'Screensize'));
    subplot(1,2,1)
    contour(x,y,abs(psi)); %contour plot of the current state
    grid on
    title('Evolution of the probability distribution $$|\psi(x,y)|$$'...
        ,'interpreter','latex')
    xlabel('$$x$$','interpreter','latex')
    ylabel('$$y$$','interpreter','latex')
    pbaspect([1 1 1])
    subplot(1,2,2)
    s1=surf(x,y,V); hold on %plot the potential
    pbaspect([1 1 1])
    alpha(s1,0.1) 
    s2=surf(x,y,7*abs(psi)); %3D plot of the current state
    pbaspect([1 1 1])
    alpha(s2,1) 
    s3=surf(x,y,V); hold off
    pbaspect([1 1 1])
    alpha(s3,0.1)
    axis([-2*N/3 2*N/3 -2*N/3 2*N/3 0 0.6]);
    title('Evolution of the probability distribution $$|\psi(x,y)|$$'...
        ,'interpreter','latex')
    xlabel('$$x$$','interpreter','latex')
    ylabel('$$y$$','interpreter','latex')
    zlabel('$$|\psi(x,y)|$$','interpreter','latex')
    pause(0.001);
    
    phi = fft2(psi);
    phi = UK_.*phi;
%     %Plot the Fourier transform of the wavefunction
% %     subplot(1,3,3)
%     contour((-64:(N-65)),(-64:(N-65)),abs(phi))
% %     axis([0 N/8 -0 N/8]);
%     pbaspect([1 1 1])
%     title(...
%%'Evolution of the Fourier transform of the probability distribution $$|\phi(k_x,k_y)|$$',...
%     'interpreter','latex')
%     xlabel('$$k_x$$','interpreter','latex')
%     ylabel('$$k_y$$','interpreter','latex')
    psi = ifft2(phi);
    psi = UV.*psi;
  
    psi_xgy=fliplr(tril(fliplr(psi)));
    
%     P_xgy=0.5*abs(sum(triu(psi)));
%     plot(-(N-1)/2:(N-1)/2,P_xgy)
%     axis([-(N-1)/2 (N-1)/2 -0 1]);
%     pause(0.001);
    
    
    %position Hilbert space
    ket=psi;
    bra=conj(psi);

%   % Calculation of dynamical quantities
%     ket=X.*ket;
%     a=sum(sum(bra.*ket));
%     avg_x(i)=real(sum(sum(X.*conj(psi).*psi))); %average value of x
    
%    %average value of x given that x>y  
   avg_x(i)=real(sum(sum(XgY.*conj(psi).*psi))); %average value of x given that x>y
%     avg_x(i)=real(sum(sum(X.*conj(psi_xgy).*psi_xgy))); 
end

% Plot <dynamical variables>

%Average position of the first particle is represented by avg_x:
p_x1=avg_x;
%Average position of the second particle has a phase shift:
p_x2=-[avg_x(360:end),avg_x(1:359)];

figure(3)
plot(p_x1,1:T,'linewidth',3,'color',[0,0.2,0.8]), hold on
plot(p_x2,1:T,'linewidth',3,'color',[1,0.4,0.2])
grid on
title('Expectation value of positions of the two particles'...
,'interpreter','latex')
ll=legend('Position of the first particle', 'Position of the second particle');
set(ll,'interpreter','latex')
ylabel('$$t$$','interpreter','latex')
xlabel('Position','interpreter','latex')