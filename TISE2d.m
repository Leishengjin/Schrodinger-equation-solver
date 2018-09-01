%% 2D time independent Schrödinger equation solver. - 17/11/2017 
% The program is not original: mostly inspired by an entry in 
% https://math.stackexchange.com/

clc; clear all;

N = 64; % Size of the grid. Assume N > 3
L = 8;  % Size of space (L*L)
h = L/N; % Size of steps

% Construct d/dx and d/dy matrix which will have size N^2 * N^2
A=diag(ones(N-1,1),1)+diag(-1*ones(N-1,1),-1);
E = eye(N);

% Do the tensor product
del_x = kron(A, E)./(2*h);
del_y = kron(E, A)./(2*h);

%laplacian operator in 2D
K=-del_x'*del_x - del_y'*del_y;

% Construct the potential matrix
v = @(x,y) 10*(x^2+y^2); %given potential 
V = zeros(N*N,N*N);
for i = 1:N
   for j = 1:N
      pos = j + N*(i - 1);
      V(pos, pos) = v(-L/2 + h*(i-1), -L/2 + h*(j-1));
   end
end

% Construct the Hamiltonian matrix
H = K + V;

    [Psi, D] = eig(H);
    eval = diag(D);

    [useless, perm] = sort(eval);
    eval = eval(perm); Psi = Psi(:,perm);

    % Plot a choice of Psi
    fprintf('Preparing the plot... \n');
    k = 9;
    PsiFn = zeros(N,N);

    for j = 1:N
        for i = 1:N
            PsiFn(i,j) = Psi(i + N*(j - 1), k + 1);
        end
    end
    [y,z] = meshgrid(-L/2:h:L/2 - h);
    % Ploting |Psi|^2
    figure;
%     surf(y,z,PsiFn);
surf(y,z,abs(PsiFn),'EdgeColor','none','LineStyle','none','FaceLighting','phong');
