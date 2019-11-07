% Visualize achievable 2D Bravais lattices

function [] = plot2d(lattice)

if strcmp(lattice,'orthorhombic-centered')
    lattice_num = 1;
elseif strcmp(lattice,'hexagonal')
    lattice_num = 2;
elseif strcmp(lattice,'tetragonal')
    lattice_num = 3;
else
    disp('Unrecognized lattice');
    return
end

% wavenumber
k = 1;
% wavelength
la = 2*pi/k;

% ARP parameters
a = 1; b = 1;

% Angle for orthorhombic centered
g_oc = pi/4;

% Angle for tetragonal
theta_t = pi/4;

% Cell array: names and wavevectors
C = {'Orthorhombic-centered','Hexagonal','Tetragonal';...
    {[csc(g_oc)*sin(g_oc/2),csc(g_oc/2)/2],[csc(g_oc)*sin(g_oc/2),...
    -csc(g_oc/2)/2]},...
    {[1,1/sqrt(3)],[0,2/sqrt(3)]},...
    {[1,0],[0,1]}};

j = lattice_num;

% wavevectors
c1 = C{2,j}{1};
c2 = C{2,j}{2};
c1 = c1./norm(c1);
c2 = c2./norm(c2);
% primitive vectors
prim = la*inv([c1', c2'])';
a1 = prim(:,1);
a2 = prim(:,2);

% Desired location of minimum
p = [0,0];

% Matrix for quadratic form
w1 = exp(1i*k*c1*p');
w2 = exp(1i*k*c2*p');
MM = [w1, w2; 1i*k*c1(1)*w1, 1i*k*c2(1)*w2; 1i*k*c1(2)*w1, 1i*k*c2(2)*w2];
Q = [real(MM), -1*imag(MM)]'*[a, 0, 0;0, -b, 0; 0, 0, -b]*...
    [real(MM), -1*imag(MM)];

% Get minimum eigenvalue of Q
[V,D] = eig(Q);
[l_0,l] = min(diag(D));

% Get amplitudes
A = [V(1,l)+1i*V(3,l),V(2,l)+1i*V(4,l)];
if strcmp(C{1,j},'Tetragonal')
    A = [1i*sin(theta_t),1i*cos(theta_t)];
end

% Calculate the field
xmin = -3*la; xmax = 3*la; nx=2048; x = linspace(xmin,xmax,nx);
ymin = -3*la; ymax = 3*la; ny=2048; y = linspace(ymin,ymax,ny);
[X,Y] = ndgrid(x,y);
u = A(1)*exp(1i*k*(c1(1)*X + c1(2)*Y)) + A(2)*exp(1i*k*(c2(1)*X +...
    c2(2)*Y));
ux = A(1)*1i*k*c1(1)*exp(1i*k*(c1(1)*X + c1(2)*Y)) +...
    A(2)*1i*k*c2(1)*exp(1i*k*(c2(1)*X + c2(2)*Y));
uy = A(1)*1i*k*c1(2)*exp(1i*k*(c1(1)*X + c1(2)*Y)) +...
    A(2)*1i*k*c2(2)*exp(1i*k*(c2(1)*X + c2(2)*Y));

% calculate acoustic radiation potential
rp = a*real(u).^2 - b*( real(ux).^2 + real(uy).^2 );
figure();
% Make blobs around minima
thresh = l_0+0.1*(max(diag(D))-l_0);
ind1 = rp<thresh;
rp(~ind1)=0;
% Change the color of blobs that don't lie on vertices
ind2 = ind1;
for ii=-5:5
    for jj=-5:5
        ind2 = ind2 & (abs(X-(ii*a1(1)+jj*a2(1)))>la/5 |...
            abs(Y-(ii*a1(2)+jj*a2(2)))>la/8);
    end
end
rp(ind2)=2*l_0;
map = [1 0 0
    0.2 0.2 1
    1 1 1];
imagesc(x/la,y/la,rp'); axis equal; axis tight; axis xy; colormap(map);
xlabel('x (in $$\ell$$)','Interpreter','latex');
ylabel('y (in $$\ell$$)','Interpreter','latex');
% Plot the primitive vectors
hold on;
corner = a1+a2;
line([0,a1(1)]/la,[0,a1(2)]/la,'Color','k','LineWidth',2);
line([0,a2(1)]/la,[0,a2(2)]/la,'Color','k','LineWidth',2);
line([a1(1),corner(1)]/la,[a1(2),corner(2)]/la,'Color','k','LineWidth',2);
line([a2(1),corner(1)]/la,[a2(2),corner(2)]/la,'Color','k','LineWidth',2);
% Plot the unit cell vectors (points counterclockwise starting at origin)
if strcmp(lattice,'orthorhombic-centered')
    upt1 = [0,0];
    upt2 = -a1-a2;
    upt3 = -2*a1;
    upt4 = a2-a1;
    line([upt2(1),upt1(1)]/la,[upt2(2),upt1(2)]/la,'Color','k',...
        'LineWidth',2,'LineStyle','--');
    line([upt2(1),upt3(1)]/la,[upt2(2),upt3(2)]/la,'Color','k',...
        'LineWidth',2,'LineStyle','--');
    line([upt3(1),upt4(1)]/la,[upt3(2),upt4(2)]/la,'Color','k',...
        'LineWidth',2,'LineStyle','--');
    line([upt4(1),upt1(1)]/la,[upt4(2),upt1(2)]/la,'Color','k',...
        'LineWidth',2,'LineStyle','--');
    start=-a1/la;
    arrow(start,start+a1/la/1.5,'width',0.5,'length',10);
    arrow(start,start+a2/la/1.5,'width',0.5,'length',10);
elseif strcmp(lattice,'hexagonal')
    upt1 = a1;
    upt2 = a1+a2;
    upt3 = a2;
    upt4 = -a1;
    upt5 = -a1-a2;
    upt6 = -a2;
    line([upt1(1),upt2(1)]/la,[upt1(2),upt2(2)]/la,'Color','k',...
        'LineWidth',2,'LineStyle','--');
    line([upt2(1),upt3(1)]/la,[upt2(2),upt3(2)]/la,'Color','k',...
        'LineWidth',2,'LineStyle','--');
    line([upt3(1),upt4(1)]/la,[upt3(2),upt4(2)]/la,'Color','k',...
        'LineWidth',2,'LineStyle','--');
    line([upt4(1),upt5(1)]/la,[upt4(2),upt5(2)]/la,'Color','k',...
        'LineWidth',2,'LineStyle','--');
    line([upt5(1),upt6(1)]/la,[upt5(2),upt6(2)]/la,'Color','k',...
        'LineWidth',2,'LineStyle','--');
    line([upt6(1),upt1(1)]/la,[upt6(2),upt1(2)]/la,'Color','k',...
        'LineWidth',2,'LineStyle','--');
    start=-a1/la;
    arrow(start,start+a1/la/1.5,'width',0.5,'length',10);
    arrow(start,start+a2/la/1.5,'width',0.5,'length',10);
    xlim([-2*norm(a1/la),2*norm(a1/la)]);
    ylim([-2*norm(a1/la),2*norm(a1/la)]);
else
    start=(a1+a2)/la/2;
    arrow(start,start+a1/la/2.5,'width',0.5,'length',8);
    arrow(start,start+a2/la/2.5,'width',0.5,'length',8);
    xlim([-norm(a1/la),2*norm(a1/la)]);
    ylim([-norm(a1/la),2*norm(a1/la)]);
end
set(gca,'FontSize',14);
filename = [C{1,j},'.png'];
print('-dpng','-r300',filename);
system(['/usr/local/bin/mogrify -trim -define png:include-chunk=none ',...
    filename])
end