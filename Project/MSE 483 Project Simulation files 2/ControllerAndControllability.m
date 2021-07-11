% Variables
g = 9.81;
m = 0.03;       % wheel weight
R = 0.04;       % wheel radius 
Jw = (m*R^2)/2; % wheel inertia moment 
M = 0.6;        % body weight 
W = 0.14;       % body width 
D = 0.04;       % body depth 
H = 0.144;      % body height 
L = H/2;        % Distance of the center of mass from the wheel axle 
Jpsi = (M*L^2)/3; % Body pitch (about y-axis) inertia moment 
Jphi = (M*((W^2)+(D^2)))/12; % Body yaw (about z-axis) inertia moment
Jm = 1e-5;      % DC motor inertia moment 
Rm = 6.69;      % DC motor resistance 
Kb = 0.468;     % DC motor back EMF constant 
Kt = 0.317;     % DC motor back torque constant
n = 1;          % gear ratio
fm = 0.0022;    % Friction coefficient between body and DC motor 
fw = 0;         % Friction coefficient between wheel and floor

% alpha and beta
alpha = (n*Kt)/Rm;
beta = ((n*Kt*Kb)/Rm)+fm;

% E matrix components
E11 = (2*m+M)*(R^2)+2*Jw+2*(n^2)*Jm;
E12 = M*L*R-2*(n^2)*Jm;
E21 = E12;
E22 = M*(L^2)+Jpsi+2*(n^2)*Jm;

% E matrix
E = [E11 E12; E21 E22];

% F matrix components
F11 = beta+fw;
F12 = -beta;
F21 = -2*beta;
F22 = 2*beta;

% G matrix components
G11 = 0;
G12 = 0;
G21 = 0;
G22 = -M*g*L;

% H Matrix components
H11 = alpha/2;
H12 = H11;
H21 = -alpha;
H22 = H21;

% A matrix components
A32 = -(M*g*L*E12)/det(E);
A42 = (M*g*L*E11)/det(E);
A33 = -((beta+fw)*E22+2*beta*E12)/det(E);
A43 = ((beta+fw)*E12+2*beta*E11)/det(E);
A34 = beta*(E22+2*E12)/det(E);
A44 = -beta*(E12+2*E11)/det(E);

% B matrix components
B3 = alpha*((E22/2)+E12)/det(E);
B4 = -alpha*((E12/2)+E11)/det(E);

% A and B matrices
A = [0 0 1 0; 0 0 0 1; 0 A32 A33 A34; 0 A42 A43 A44]
B = [0;0;B3;B4]
C = [1 0 0 0];

% 5x5 A and B matrices
Anew = [0 0 1 0 0; 0 0 0 1 0; 0 A32 A33 A34 0; 0 A42 A43 A44 0; -1 0 0 0 0];
Bnew = [0;0;B3;B4;0];
C = [1 0 0 0];


% Pole placement using eigenvalues from criteria, A, and B
K = place(Anew,Bnew,[-1+0.805i -1-0.805i -10 -15 -100]);
Ac = Anew-Bnew*K;
Bc = [0;0;0;0;1000]; % Chosen as [0;0;0;0;1]r but simplified as such where r=1000 here
Cc = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0];
Dc = [0;0;0;0];
closedsys = ss(Ac,Bc,Cc,Dc);

% Plotting controlled step output (Out1=x1 | Out2=x2 | Out3=x3 | Out4=x4)
% x1 = theta
% x2 = phi
% x3 = theta dot
% x4 = phi dot
step(closedsys)

% Coefficients throughout, denoted by Z
Za = A32*A33+A42*A34;
Zb = (A33^2)+A43*A34;
Zc = A32+A33*A34+A34*A44;
Zd = A32*A43+A42*A44;
Ze = A33*A43+A43*A44;
Zf = A42+A43*A34+(A44^2);
Zg = Za;
Zh = Zb;
Zi = Zc;
Zj = Zd;
Zk = Ze;
Zl = Zf;
Zm = A32*Zb+A42*Zc;
Zn = A33*Zb+A43*Zc;
Zp = Za+A34*Zb+A44*Zc;
Zq = A32*Ze+A42*Zf;
Zr = A33*Ze+A43*Zf;
Zs = Zd+A34*Ze+A44*Zf;

% Controllability Matrix P
P11 = 0;
P12 = B3;
P13 = B3*A33+B4*A34;
P14 = B3*Zb+B4*Zc;
P21 = 0;
P22 = B4;
P23 = B3*A43+B4*A44;
P24 = B3*Ze+B4*Zf;
P31 = B3;
P32 = B3*A33+B4*A34;
P33 = B3*Zb+B4*Zc;
P34 = B3*Zn+B4*Zp;
P41 = B4;
P42 = B3*A43+B4*A44;
P43 = B3*Ze+B4*Zf;
P44 = B3*Zr+B4*Zs;

% Determining controllability
P = [P11 P12 P13 P14; P21 P22 P23 P24; P31 P32 P33 P34; P41 P42 P43 P44]
det_P = det(P)
rank_P = rank(P)
Co = ctrb(A,B)