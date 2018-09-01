%% Doruk Efe Gökmen -- 17/11/2017 -- 1D Schrödinger Equation Solver v1
% The time evolution of a given initial state is calculated by iteratively
% acting the time evolution operator on the state vector.
% The time evolution matrix is assumed to decomposable into kinetic energy 
% (K) and potential energy (V) parts with small error, given that the time 
% step is sufficiently small, although [K,V]=/=0.
% Kinetic energy matrix is diagonalised using the fast Fourier
% transform (fft) algorithm and the corresponding part of the time
% evolution operator is acted upon the state vector which is also
% transformed into the reciprocal space representation by taking fft. The
% modified state is transformed back into the position space representation
% by taking the inverse fast Fourier transform (ifft) and finally acted by
% the potential part of the time evolution operator to complete an
% iteration. Note that since all of these operations are unitary, the
% normalisation and the energy are conserved in this numerical scheme.

% In v1, instead of the matrix multiplications term by term vector 
% multiplications are used by taking the advantage that all matrices are 
% diagonal.
clf;
clc;
clear;
% Initialising parameters
T = 6000; %time length
N = 256; %number of spatial elements
dx = 1; %spatial integration step
hbar = 1;
m = 1; %mass of the particle
dt = 10e-2; %temporal integration step
k = 2*pi/N * (0:N-1); %reciprocal space

% Enter the initial wavefunction (state)
p0=0;%2*sin(k(50) / 2); %initial momentum (arbitrarily chosen)
TP = exp(-1i*p0*(0:(N-1)))'; %momentum translation operator
%psi = zeros(N,1);
psi = TP .* (exp((-((0:(N-1))-(N-4)/2).^2)/(N/20)^2))'; %Gaussian initial wavefunction
psi = psi / sqrt(dx*sum(psi'*psi)); %normalise the wavefunction

%Potential energy matrix in position space

% Potential well types -----------------------------------------------
disp('  ');
disp('  ');
disp('   POTENTIAL WELL TYPES - enter 1 or 2 or 3 ... or 6');
disp('  ');
disp('   1: Free particle');
disp(' ');
disp('   2: Harmonic well');
disp(' ');
disp('   3: Rectangular box');
disp(' ');
disp('   4: Rectangular box with a delta barrier in the middle');
disp(' ');
disp('   5: Periodic rectangular well lattice');
disp(' ');
disp('   6: Sinusoidal lattice');
disp(' ');
PotentialType = input('Specify well type: 1, 2, 3, 4, 5, 6 :');
disp(' ');
disp(' ');

switch PotentialType
    case 1
        v = zeros(N,1)'; %"null potential"
        % Free particle
        %V = diag(v);
    case 2
        k0 = 50e-5; % Harmonic potential well
        v = 0.5*k0*((0:N-1)-N/2).^2;
    case 3
        v = zeros(N,1)'; %"null potential"
        v(1:N/4) = 10e5; % Particle in a rectangular box
        v(N-N/4:N) = 10e5;
    case 4
        v = zeros(N,1)'; %"null potential"
        v(1:N/4) = 10e5; % Particle in a rectangular box
        v(N-N/4:N) = 10e5;
        v(N/2) = 10e5; %with a delta-barrier in the middle
    case 5
        v = 0.3 * square(0.04*pi*((0:(N-1))-N/2)); %periodic square well lattice
    case 6
        v = 0.4 * sin(0.1*(0:(N-1))); %sinusoidal lattice
end

% V = diag(v); %converts the potential vector into a diagonal matrix 
%(uncomment if working with matrices)

% Kinetic energy matrix in reciprocal space
KK_=hbar^2 * 4*sin(k / 2).^2/(2*m);
% K_ = diag(KK_); %(uncomment if working with matrices)

% Time evolution operator

%Time evolution operator can be decomposed into kinetic and potential parts
%with small error assuming that dt is sufficiently small.

%%Uncomment if you choose to work with matrices (intuitive)
% UV = expm(-1i*dt/hbar *V); %potential energy part (diagonal in position space)
% UK_ = expm(-1i*dt/hbar *K_); %kinetic energy part (reciprocal space representation)
%%BUT working with matrices is slow and unnecessary since they are diagonal
%%NB: Do not forget to switch the .*'s to *'s in the below for loop if you
%%choose to work with matrices.

% Let us work with vectors instead of matrices (efficient)
UV = exp(-1i*dt/hbar *v'); %take exponential of each element
UK_ = exp(-1i*dt/hbar *KK_'); 

avg_x=zeros(1,T);
avg_p=avg_x;
avg_KE=avg_x;

for i = 1:T
%     plot(v,"linewidth",4); hold on %plot the potential
%     plot(abs(psi),"linewidth",1); %plot the current state
%     plot(v,"linewidth",4); hold off
%     axis([0 N -0.05 1]);
%     pause(0.001);
    phi = fft(psi);
    phi = UK_.*phi; %term by term multiplication since vectors are used
    psi = ifft(phi);
    psi = UV.*psi; %same reasoning for .* as the previous line
%   psi = U*psi;

%     plot(abs(phi),"linewidth",1)
    
    %Calculate expectation values of dynamic variables
    
    %operators
    x_op = (0:(N-1))'; %x (position) operator in position space representation
    p_op = 2*sin(k)'; %p (momentum) operator in reciprocal space representation
    xx_op = x_op .* x_op; %x^2 operator in position space
    KE_op = p_op .* p_op /(2*m); % kinetic energy operator in reciprocal space
    
    %position Hilbert space
    ket=psi;
    bra=psi';
    ket=x_op .* ket;
    avg_x(i)=real(bra*ket); %Expectation value of x
    
    %reciprocal Hilbert space
    ket_=phi;
    bra_=phi';
    
%     ket_=p_op .* ket_;
%     avg_p(i)=real(bra_*ket_); %Expectation value of p
    
    ket_=KE_op .* ket_;
    avg_KE(i)=real(bra_*ket_); %Expectation value of kinetic energy
end

figure(2)
subplot(1,3,1)
plot(1:T,avg_x-N/2,1:T,-mean(avg_p)+avg_p,"linewidth",1)
title('Expectation value of position'...
,'interpreter','latex')
xlabel('$$t$$','interpreter','latex')
ylabel('Position $$\langle x \rangle$$','interpreter','latex')
pbaspect([2 1 1])

subplot(1,3,2)
plot(-mean(avg_p)+avg_p,"linewidth",1)
title('Expectation value of momentum'...
,'interpreter','latex')
xlabel('$$t$$','interpreter','latex')
ylabel('Momentum $$\langle p \rangle$$','interpreter','latex')
pbaspect([2 1 1])

subplot(1,3,3)
plot(avg_KE,"linewidth",1)
title('Expectation value of kinetic energy'...
,'interpreter','latex')
xlabel('$$x$$','interpreter','latex')
ylabel('Kinetic Energy $$\langle \frac{p^2}{2m} \rangle$$','interpreter','latex')
pbaspect([2 1 1])