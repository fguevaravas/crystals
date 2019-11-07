% Another example of ARP corresponding to lines (example 2.2)

% ARP parameters
a = 1; b = 0;

% Directions
t1 = 0; t2 = pi/2;

% Location of mininmum
p = [0,0];

k1 = [cos(t1);sin(t1)];
k2 = [cos(t2);sin(t2)];
k1=k1./norm(k1);
k2=k2./norm(k2);
la = 2*pi/norm(k1);

% projectors
Vp = [eye(2); eye(2)]/sqrt(2);
Vm = [eye(2);-eye(2)]/sqrt(2);
Pp = Vp*Vp';
Pm = Vm*Vm';

dd = la*inv([k1,k2])'; % dual vectors

Mp = [exp(1i*p*k1),exp(1i*p*k2); 1i*k1*exp(1i*p*k1), 1i*k2*exp(1i*p*k2)];
Mn = [exp(-1i*p*k1),exp(-1i*p*k2);...
    -1i*k1*exp(-1i*p*k1), -1i*k2*exp(-1i*p*k2)];
M = [Mp, Mn];
Q = M'*[a,0,0;0,-b,0;0,0,-b]*M;

[V,D] = eig(Q);
[~,l] = min(diag(D));

K = [k1, k2];
[KV,KD] = eig(K*K');
sigma_1 = KD(1,1);
sigma_2 = KD(2,2);
[ZV,~] = eig(ones(2,2));
ZV = fliplr(ZV);
U = 1/sqrt(2)*[ZV,KV;ZV,-KV];
La = [4*a,0,-2*b*sigma_1,-2*b*sigma_2];
alpha_1 = U(1,2);
alpha_2 = U(2,2);
beta_1 = U(3,2);
beta_2 = U(4,2);

theta = -4:0.5:4;
num_n = 2;

A = inv(K);

xmin = -1*la; xmax = 2*la; nx=200; x = linspace(xmin,xmax,nx); 
ymin = -1*la; ymax = 2*la; ny=200; y = linspace(ymin,ymax,ny);
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

figure();
imagesc(x/la,y/la,rp'); axis equal; axis tight; axis xy;
xlabel('x (in $$\ell$$)','Interpreter','latex'); 
ylabel('y (in $$\ell$$)','Interpreter','latex'); 
colorbar;
hold on;
L = zeros(2,length(theta));
LL = zeros(2,length(theta));
Lnew = zeros(2,(2*num_n+1)*(2*num_n+1)*length(theta));
for n1=-num_n:1:num_n
    for n2=-num_n:1:num_n
        ind = 1;
        for i=1:length(theta)
            L(:,ind) = K'\(theta(i)*[1;-1]-2*pi*[n1;n2]);
            LL(:,ind) = K'\(theta(i)*[1;1]-2*pi*[n1;n2]);
            ind = ind + 1;
        end
        plot(L(1,:)/la,L(2,:)/la,'r-',LL(1,:)/la,LL(2,:)/la,'r-',...
            'LineWidth',2);
    end
end

a1 = A(:,1);
a2 = A(:,2);
corner = a1+a2;
line([0,a1(1)],[0,a1(2)],'Color','w','LineWidth',2);
line([0,a2(1)],[0,a2(2)],'Color','w','LineWidth',2);
line([a1(1),corner(1)],[a1(2),corner(2)],'Color','w','LineWidth',2);
line([a2(1),corner(1)],[a2(2),corner(2)],'Color','w','LineWidth',2);
xlim([-1,2]); ylim([-1,2]);
start = [1/2,1/2];
arrow(start,start+k1'/2,'width',0.5,'length',10);
arrow(start,start+k2'/2','width',0.5,'length',10);
set(gca,'FontSize',14);
filename = 'arp_lines.png';
print('-dpng','-r300',filename);
system(['/usr/local/bin/mogrify -trim -define png:include-chunk=none ',...
    filename])