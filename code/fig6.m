% Minima of ARP in 3D on lines (example 2.4)

% ARP parameters
a = 1; b = 0;
% wavelength
k = 1;
la = 2*pi/k;
% column to take amplitudes from
col = 3;

% color for patch
if col==3
    patch_color = 'red';
elseif col==2
    patch_color = 'none';
end

% Location of mininmum
p = [0,0,0];

% Reciprocal vectors
K = eye(3);
k1 = K(:,1);
k2 = K(:,2);
k3 = K(:,3);

Mp = [exp(1i*p*k1),exp(1i*p*k2),exp(1i*p*k3);...
    1i*k1*exp(1i*p*k1), 1i*k2*exp(1i*p*k2), 1i*k3*exp(1i*p*k3)];
Mn = [exp(-1i*p*k1),exp(-1i*p*k2),exp(-1i*p*k3);...
    -1i*k1*exp(-1i*p*k1), -1i*k2*exp(-1i*p*k2), -1i*k3*exp(-1i*p*k3)];
M = [Mp, Mn];
Q = M'*diag([a,-b,-b,-b])*M;

[V,D] = eig(Q);
[~,l] = min(diag(D));
[KV,KD] = eig(K*K');
sigma_1 = KD(1,1);
sigma_2 = KD(2,2);
sigma_3 = KD(3,3);
[ZV,~] = eig(ones(3,3));
ZV = fliplr(ZV);
U = 1/sqrt(2)*[ZV,KV;ZV,-KV];
La = [6*a,0,0,-2*b*sigma_1,-2*b*sigma_2,-2*b*sigma_3];
alpha_1 = U(1,col);
alpha_2 = U(2,col);
alpha_3 = U(3,col);
beta_1 = U(4,col);
beta_2 = U(5,col);
beta_3 = U(6,col);

A = la*inv(K);
a1 = A(:,1);
a2 = A(:,2);
a3 = A(:,3);

xmin = -1*la; xmax = 2*la; nx=50; x = linspace(xmin,xmax,nx); 
ymin = -1*la; ymax = 2*la; ny=50; y = linspace(ymin,ymax,ny);
zmin = -1*la; zmax = 2*la; nz=50; z = linspace(zmin,zmax,nz);
[X,Y,Z] = ndgrid(x,y,z);

u = alpha_1*exp(1i*(k1(1)*X + k1(2)*Y + k1(3)*Z)) +...
    alpha_2*exp(1i*(k2(1)*X + k2(2)*Y + k2(3)*Z)) +...
    alpha_3*exp(1i*(k3(1)*X + k3(2)*Y + k3(3)*Z)) +...
    beta_1*exp(-1i*(k1(1)*X + k1(2)*Y + k1(3)*Z)) +...
    beta_2*exp(-1i*(k2(1)*X + k2(2)*Y + k2(3)*Z)) +...
    beta_3*exp(-1i*(k3(1)*X + k3(2)*Y + k3(3)*Z));
ux = alpha_1*1i*k1(1)*exp(1i*(k1(1)*X + k1(2)*Y + k1(3)*Z)) +...
    alpha_2*1i*k2(1)*exp(1i*(k2(1)*X + k2(2)*Y + k2(3)*Z)) +...
    alpha_3*1i*k3(1)*exp(1i*(k3(1)*X + k3(2)*Y + k3(3)*Z)) -...
    beta_1*1i*k1(1)*exp(-1i*(k1(1)*X + k1(2)*Y + k1(3)*Z)) -...
    beta_2*1i*k2(1)*exp(-1i*(k2(1)*X + k2(2)*Y + k2(3)*Z)) -...
    beta_3*1i*k3(1)*exp(-1i*(k3(1)*X + k3(2)*Y + k3(3)*Z));
uy = alpha_1*1i*k1(2)*exp(1i*(k1(1)*X + k1(2)*Y + k1(3)*Z)) +...
    alpha_2*1i*k2(2)*exp(1i*(k2(1)*X + k2(2)*Y + k2(3)*Z)) +...
    alpha_3*1i*k3(2)*exp(1i*(k3(1)*X + k3(2)*Y + k3(3)*Z)) -...
    beta_1*1i*k1(2)*exp(-1i*(k1(1)*X + k1(2)*Y + k1(3)*Z)) -...
    beta_2*1i*k2(2)*exp(-1i*(k2(1)*X + k2(2)*Y + k2(3)*Z)) -...
    beta_3*1i*k3(2)*exp(-1i*(k3(1)*X + k3(2)*Y + k3(3)*Z));
uz = alpha_1*1i*k1(3)*exp(1i*(k1(1)*X + k1(2)*Y + k1(3)*Z)) +...
    alpha_2*1i*k2(3)*exp(1i*(k2(1)*X + k2(2)*Y + k2(3)*Z)) +...
    alpha_3*1i*k3(3)*exp(1i*(k3(1)*X + k3(2)*Y + k3(3)*Z)) -...
    beta_1*1i*k1(3)*exp(-1i*(k1(1)*X + k1(2)*Y + k1(3)*Z)) -...
    beta_2*1i*k2(3)*exp(-1i*(k2(1)*X + k2(2)*Y + k2(3)*Z)) -...
    beta_3*1i*k3(3)*exp(-1i*(k3(1)*X + k3(2)*Y + k3(3)*Z));

rp = a*abs(u).^2 - b*( abs(ux).^2 + abs(uy).^2 + abs(uz).^2);

min_ind = rp<(min(rp(:)) + 0.03*(max(rp(:))-min(rp(:))));

svec = [0 0 0 1 0 1 1 1
        0 0 1 0 1 0 1 1
        0 1 0 0 1 1 0 1];
    
rvec = svec;

uvec = U(:,col);
vvec = uvec(1:3);
tol = 10e-8;
i = 0;
for k = 1:8
    for j = 1:8
        s = svec(:,k);
        r = rvec(:,j);
        u1 = [(-1).^s; (-1).^s].*[vvec;-vvec];
        u2 = [(-1).^r; (-1).^r].*[vvec;-vvec];
        u3 = [(-1).^(s+r); (-1).^(s+r)].*[vvec;vvec];
        if (col==3)
            if((norm(Q*u1)<tol)&&...
                (norm(Q*u2)<tol)&&...
                (norm(Q*u3)<tol))
                i = i + 1;
                disp([s,r]);
                ss(:,i) = s;
                rr(:,i) = r;
            end
        elseif (col==2)
            if((norm(Q*u1)<tol)&&...
                (norm(Q*u2)<tol)&&...
                (norm(Q*u3)<tol)&&...
                (norm(cross((-1).^s,(-1).^r))>tol))
                i = i + 1;
                disp([s,r]);
                ss(:,i) = s;
                rr(:,i) = r;
            end
        end
    end
end

theta1 = -40; phi1 = -40;
theta2 = -40; phi2 = 40;
theta3 = 40; phi3 = 40;
theta4 = 40; phi4 = -40;

figure(); hold on;

for i = 1:size(ss,2)
    for n1 = -3:3
        for n2 = -3:3
            for n3 = -3:3
                s = ss(:,i);
                r = rr(:,i);
                n = [n1;n2;n3];
                p1 = K'\(theta1*(-1).^s + phi1*(-1).^r + 2*pi*n);
                p2 = K'\(theta2*(-1).^s + phi2*(-1).^r + 2*pi*n);
                p3 = K'\(theta3*(-1).^s + phi3*(-1).^r + 2*pi*n);
                p4 = K'\(theta4*(-1).^s + phi4*(-1).^r + 2*pi*n);
    
                xx = [p1(1),p2(1),p3(1),p4(1)];
                yy = [p1(2),p2(2),p3(2),p4(2)];
                zz = [p1(3),p2(3),p3(3),p4(3)];
                
                patch(xx/la,yy/la,zz/la,'r','edgecolor',patch_color,...
                    'facealpha',0.25,'edgealpha',1.0,'clipping','on');
            end
        end
    end
end

s1_x = linspace(0,a1(1),200); 
s1_y = linspace(0,a1(2),200); 
s1_z = linspace(0,a1(3),200);
plot3(s1_x/la,s1_y/la,s1_z/la,'k.');
s2_x = linspace(0,a2(1),200); 
s2_y = linspace(0,a2(2),200); 
s2_z = linspace(0,a2(3),200);
plot3(s2_x/la,s2_y/la,s2_z/la,'k.');
s3_x = linspace(0,a3(1),200); 
s3_y = linspace(0,a3(2),200); 
s3_z = linspace(0,a3(3),200);
plot3(s3_x/la,s3_y/la,s3_z/la,'k.');
s4_x = linspace(a1(1),a1(1)+a2(1),200); 
s4_y = linspace(a1(2),a1(2)+a2(2),200); 
s4_z = linspace(a1(3),a1(3)+a2(3),200);
plot3(s4_x/la,s4_y/la,s4_z/la,'k.');
s5_x = linspace(a1(1),a1(1)+a3(1),200); 
s5_y = linspace(a1(2),a1(2)+a3(2),200); 
s5_z = linspace(a1(3),a1(3)+a3(3),200);
plot3(s5_x/la,s5_y/la,s5_z/la,'k.');
s6_x = linspace(a2(1),a1(1)+a2(1),200); 
s6_y = linspace(a2(2),a1(2)+a2(2),200); 
s6_z = linspace(a2(3),a1(3)+a2(3),200);
plot3(s6_x/la,s6_y/la,s6_z/la,'k.');
s7_x = linspace(a2(1),a2(1)+a3(1),200); 
s7_y = linspace(a2(2),a2(2)+a3(2),200); 
s7_z = linspace(a2(3),a2(3)+a3(3),200);
plot3(s7_x/la,s7_y/la,s7_z/la,'k.');
s8_x = linspace(a3(1),a1(1)+a3(1),200); 
s8_y = linspace(a3(2),a1(2)+a3(2),200); 
s8_z = linspace(a3(3),a1(3)+a3(3),200);
plot3(s8_x/la,s8_y/la,s8_z/la,'k.');
s9_x = linspace(a3(1),a2(1)+a3(1),200); 
s9_y = linspace(a3(2),a2(2)+a3(2),200); 
s9_z = linspace(a3(3),a2(3)+a3(3),200);
plot3(s9_x/la,s9_y/la,s9_z/la,'k.');
s10_x = linspace(a1(1)+a2(1),a1(1)+a2(1)+a3(1),200); 
s10_y = linspace(a1(2)+a2(2),a1(2)+a2(2)+a3(2),200); 
s10_z = linspace(a1(3)+a2(3),a1(3)+a2(3)+a3(3),200);
plot3(s10_x/la,s10_y/la,s10_z/la,'k.');
s11_x = linspace(a1(1)+a3(1),a1(1)+a2(1)+a3(1),200); 
s11_y = linspace(a1(2)+a3(2),a1(2)+a2(2)+a3(2),200); 
s11_z = linspace(a1(3)+a3(3),a1(3)+a2(3)+a3(3),200);
plot3(s11_x/la,s11_y/la,s11_z/la,'k.');
s12_x = linspace(a2(1)+a3(1),a1(1)+a2(1)+a3(1),200); 
s12_y = linspace(a2(2)+a3(2),a1(2)+a2(2)+a3(2),200); 
s12_z = linspace(a2(3)+a3(3),a1(3)+a2(3)+a3(3),200);
plot3(s12_x/la,s12_y/la,s12_z/la,'k.');
xlabel('x (in $$\ell$$)','Interpreter','latex');
ylabel('y (in $$\ell$$)','Interpreter','latex');
zlabel('z (in $$\ell$$)','Interpreter','latex');

azim = -13; elev = 26;
axis equal;
axis([-2 2 -2 2 0 1])
light; view(azim, elev);

start = [-1,-1,0.1];
arrow3(start,start+k1'/1.2,'*',1,2);
arrow3(start,start+k2'/1.2,'*',1,2);
arrow3(start,start+k3'/1.2);

set(gca,'FontSize',14);
if col==2
    filename = 'arp_planes.png';
elseif col==3
    filename = 'arp_lines_3D.png';
end

if col==3
theta1 = -10; phi1 = -10;
theta2 = -10; phi2 = 10;
theta3 = 10; phi3 = 10;
theta4 = 10; phi4 = -10;
for i = 1:size(ss,2)
    for n1 = -1:1
        for n2 = -1:1
            for n3 = -1:1
                s = ss(:,i);
                r = rr(:,i);
                n = [n1;n2;n3];
                p1 = K'\(theta1*(-1).^s + phi1*(-1).^r + 2*pi*n);
                p2 = K'\(theta2*(-1).^s + phi2*(-1).^r + 2*pi*n);
                p3 = K'\(theta3*(-1).^s + phi3*(-1).^r + 2*pi*n);
                p4 = K'\(theta4*(-1).^s + phi4*(-1).^r + 2*pi*n);
    
                xx = [p1(1),p2(1),p3(1),p4(1)];
                yy = [p1(2),p2(2),p3(2),p4(2)];
                zz = [p1(3),p2(3),p3(3),p4(3)];
                
                patch(xx/la,yy/la,zz/la,'r','edgecolor',patch_color,...
                    'LineWidth',1,'facealpha',0.5,'clipping','on');
            end
        end
    end
end
end

print('-dpng','-r300',filename);
system(['/usr/local/bin/mogrify -trim -define png:include-chunk=none ',...
    filename])

% Restrict to unit cell
axis([0 1 0 1 0 1])
azim = -15; elev = 29; view(azim, elev);
if col==2
    filename = 'arp_planes_cell.png';
elseif col==3
    filename = 'arp_lines_3D_cell.png';
end
print('-dpng','-r300',filename);
system(['/usr/local/bin/mogrify -trim -define png:include-chunk=none ',...
    filename])