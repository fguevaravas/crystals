% Example of ARP corresponding to lines (example 2.1)

% wavenumber
k = 1;
% wavelength
la = 2*pi/k;

% ARP parameters
a = 1; b = 1;

% Angle for tetragonal
theta_t = pi/4;

% wavevectors
c1 = [1,0]; c2=[0,1];

% primitive vectors
prim = la*inv([c1', c2'])';
a1 = prim(:,1);
a2 = prim(:,2);
A = [a1,a2]';

% Desired location of minimum
p = [0,0];

% Matrix for quadratic form
w1 = exp(1i*k*c1*p');
w2 = exp(1i*k*c2*p');
w1m = exp(-1i*k*c1*p');
w2m = exp(-1i*k*c2*p');
Mp = [w1, w2; 1i*k*c1(1)*w1, 1i*k*c2(1)*w2; 1i*k*c1(2)*w1, 1i*k*c2(2)*w2];
Mm = [w1m, w2m; -1i*k*c1(1)*w1m, -1i*k*c2(1)*w2m;...
    -1i*k*c1(2)*w1m, -1i*k*c2(2)*w2m];
Q = [Mp, Mm]'*[a, 0, 0;0, -b, 0; 0, 0, -b]*[Mp, Mm];

% Get minimum eigenvalue of Q
[V,D] = eig(Q);
[l_0,l] = min(diag(D));

K = eye(2); k1 = K(:,1); k2 = K(:,2);
[KV,KD] = eig(K*K');
sigma_1 = KD(1,1);
sigma_2 = KD(2,2);
[ZV,~] = eig(ones(2,2));
ZV = fliplr(ZV);
U = 1/sqrt(2)*[ZV,KV;ZV,-KV];
La = [4*a,0,-2*b*sigma_1,-2*b*sigma_2];

% amplitudes
alpha_1 = 1/sqrt(2);
alpha_2 = 0;
beta_1 = -1/sqrt(2);
beta_2 = 0;

% Calculate the field
xmin = -1*la; xmax = 2*la; nx=500; x = linspace(xmin,xmax,nx);
ymin = -1*la; ymax = 2*la; ny=500; y = linspace(ymin,ymax,ny);
[X,Y] = ndgrid(x,y);
u = alpha_1*exp(1i*(k1(1)*X + k1(2)*Y)) +...
    alpha_2*exp(1i*(k2(1)*X + k2(2)*Y)) +...
    beta_1*exp(-1i*(k1(1)*X + k1(2)*Y)) +...
    beta_2*exp(-1i*(k2(1)*X + k2(2)*Y));
ux = alpha_1*1i*k1(1)*exp(1i*(k1(1)*X + k1(2)*Y)) +...
    alpha_2*1i*k2(1)*exp(1i*(k2(1)*X + k2(2)*Y)) -...
    beta_1*1i*k1(1)*exp(-1i*(k1(1)*X + k1(2)*Y)) -...
    beta_2*1i*k2(1)*exp(-1i*(k2(1)*X + k2(2)*Y));
uy = alpha_1*1i*k1(2)*exp(1i*(k1(1)*X + k1(2)*Y)) +...
    alpha_2*1i*k2(2)*exp(1i*(k2(1)*X + k2(2)*Y)) -...
    beta_1*1i*k1(2)*exp(-1i*(k1(1)*X + k1(2)*Y)) -...
    beta_2*1i*k2(2)*exp(-1i*(k2(1)*X + k2(2)*Y));

% calculate acoustic radiation potential
rp = a*abs(u).^2 - b*( abs(ux).^2 + abs(uy).^2 );
imagesc(x/la,y/la,rp'); colorbar; axis equal; axis tight; axis xy;
hold on;
xlabel('x (in $$\ell$$)','Interpreter','latex'); 
ylabel('y (in $$\ell$$)','Interpreter','latex');
colorbar;
corner = a1+a2;
line([0,a1(1)]/la,[0,a1(2)]/la,'Color','w','LineWidth',2);
line([0,a2(1)]/la,[0,a2(2)]/la,'Color','w','LineWidth',2);
line([a1(1),corner(1)]/la,[a1(2),corner(2)]/la,'Color','w','LineWidth',2);
line([a2(1),corner(1)]/la,[a2(2),corner(2)]/la,'Color','w','LineWidth',2);

for k=-1:0.5:2
    line([k,k],[-1,2],'Color','r','LineWidth',2);
end
start = [1/2,1/2];
arrow(start,start+c1/2,'width',0.5,'length',10);
arrow(start,start+c2/2,'width',0.5,'length',10);

xlim([-1,2]); ylim([-1,2]);
set(gca,'FontSize',14);
filename = 'tetragonal_arp.png';
print('-dpng','-r300',filename);
system(['/usr/local/bin/mogrify -trim -define png:include-chunk=none ',...
    filename])