% Visualize achievable 3D Bravais lattices

function [] = plot3d(lattice)

if strcmp(lattice,'triclinic-primitive')
    lattice_num = 1;
    dist_pts = 1.2;
    nx = 300; ny = 300; nz = 300;
    thresh_const = 0.05;
    azim = -45; elev = -1;
    num_blue = 8;
    num_red = 1;
    start = [-1.5,2.5,0];
    arrow_scale = 1/2;
    arrow_width = 0.7;
    arrow_length = 1.4;
    r = 0.3;
    light_position = [-1,-1,1];
elseif strcmp(lattice,'orthorhombic-face-centered')
    lattice_num = 2;
    dist_pts = 1.2;
    nx = 250; ny = 250; nz = 250;
    thresh_const = 0.05;
    azim = 59; elev = 18;
    num_blue = 8;
    num_red = 1;
    start = [-1,-2,3.5];
    arrow_scale = 1.2;
    arrow_width = 0.8;
    arrow_length = 1.6;
    r = 0.2;
    light_position = [1,-1,1];
elseif strcmp(lattice,'tetragonal-body-centered')
    lattice_num = 3;
    dist_pts = 1.2;
    nx = 350; ny = 350; nz = 350;
    thresh_const = 0.05;
    azim = 62; elev = 12;
    num_blue = 8;
    num_red = 1;
    start = [-1/sqrt(2),-1/sqrt(2),-1/sqrt(2)];
    arrow_scale = 3/2;
    arrow_width = 1;
    arrow_length = 2;
    r = 0.1;
    light_position = [1,0,1];
elseif strcmp(lattice,'trigonal-primitive')
    lattice_num = 4;
    dist_pts = 1.2;
    nx = 250; ny = 250; nz = 250;
    thresh_const = 0.03;
    azim = 79; elev = 53;
    num_blue = 8;
    num_red = 1;
    start = [0.5,0.25,0];
    arrow_scale = 2;
    arrow_width = 1;
    arrow_length = 2;
    r = 0.15;
    light_position = [1,0,1];
elseif strcmp(lattice,'cubic-primitive')
    lattice_num = 5;
    dist_pts = 1.1;
    nx = 250; ny = 250; nz = 250;
    thresh_const = 0.02;
    azim = 58; elev = 8;
    num_blue = 8;
    num_red = 19;
    start = [-0.2,-0.2,-0.2];
    arrow_scale = 3;
    arrow_width = 1;
    arrow_length = 2;
    r = 0.075;
    light_position = [1,0,0.8];
elseif strcmp(lattice,'cubic-face-centered')
    lattice_num = 6;
    dist_pts = 1.2;
    nx = 250; ny = 250; nz = 250;
    thresh_const = 0.03;
    azim = 141; elev = 18;
    num_blue = 8;
    num_red = 1;
    start = [1.7,0,0];
    arrow_scale = 2;
    arrow_width = 0.9;
    arrow_length = 1.8;
    r = 0.1;
    light_position = [1,0,1];
elseif strcmp(lattice,'cubic-body-centered')
    lattice_num = 7;
    dist_pts = 1.2;
    nx = 250; ny = 250; nz = 250;
    thresh_const = 0.05;
    azim = 58; elev = 8;
    num_blue = 8;
    num_red = 1;
    start = [-1/sqrt(2),-1/sqrt(2),-1/sqrt(2)];
    arrow_scale = 1.8;
    arrow_width = 0.8;
    arrow_length = 1.6;
    r = 0.1;
    light_position = [1,0,1]; 
else
    disp('Unrecognized lattice');
    return
end

% wavenumber
k = 2;
% wavelength
la = 2*pi/k;

% ARP parameters
a = 1; b = 1;

% Triclinic primitive vectors
c1_tp = [1,2,7]; c2_tp = [8,3,5]; c3_tp = [1,3,5];
% Parameters for orthorhombic face-centered
a_ofc = 1; b_ofc = 2; c_ofc = 3;
% Parameters for trigonal primitive
a_tp = 1; c_tp = 2;

C = {'Triclinic-primitive','Orthorhombic-face-centered',...
    'Tetragonal-body-centered','Trigonal-primitive','Cubic-primitive',...
    'Cubic-face-centered','Cubic-body-centered';...
    {c1_tp,c2_tp,c3_tp},{[1/a_ofc,1/b_ofc,1/c_ofc],[-1/a_ofc,-1/b_ofc,1/c_ofc],...
    [1/a_ofc,-1/b_ofc,-1/c_ofc]},...
    {[0,1,1],[1,0,1],[1,1,0]},...
    {[0,-2/(3*a_tp),1/(3*c_tp)],[1/(sqrt(3)*a_tp),1/(3*a_tp),1/(3*c_tp)],...
    [-1/(sqrt(3)*a_tp),1/(3*a_tp),1/(3*c_tp)]},...
    {[1,0,0],[0,1,0],[0,0,1]},{[-1,1,1],[1,-1,1],[1,1,-1]},...
    {[0,1,1],[1,0,1],[1,1,0]}};

j = lattice_num;

% wavevectors
c1 = C{2,j}{1};
c2 = C{2,j}{2};
c3 = C{2,j}{3};
c1 = c1./norm(c1);
c2 = c2./norm(c2);
c3 = c3./norm(c3);

% primitive vectors
prim = la*inv([c1',c2',c3'])';
a1 = prim(:,1);
a2 = prim(:,2);
a3 = prim(:,3);

% Desired location of minimum
p = [0,0,0];

% Matrix for quadratic form
w1 = exp(1i*k*c1*p');
w2 = exp(1i*k*c2*p');
w3 = exp(1i*k*c3*p');
MM = [w1, w2, w3; 1i*k*c1(1)*w1, 1i*k*c2(1)*w2, 1i*k*c3(1)*w3;...
    1i*k*c1(2)*w1, 1i*k*c2(2)*w2, 1i*k*c3(2)*w3;...
    1i*k*c1(3)*w1, 1i*k*c2(3)*w2, 1i*k*c3(3)*w3];
Q = [real(MM), -1*imag(MM)]'*[a,0,0,0;0,-b,0,0;0,0,-b,0;0,0,0,-b]*...
    [real(MM), -1*imag(MM)];

% Get minimum eigenvalue of Q
[V,D] = eig(Q);
[l_0,l] = min(diag(D));

% Get amplitudes
A = [V(1,l)+1i*V(4,l),V(2,l)+1i*V(5,l),V(3,l)+1i*V(6,l)];
if strcmp(C{1,j},'Cubic-face-centered')
    A = [V(1,2)+1i*V(4,2),V(2,2)+1i*V(5,2),V(3,2)+1i*V(6,2)];
end
if strcmp(C{1,j},'Cubic-primitive')
    A = [1i,1i,1i]/sqrt(3);
end

% Calculate the field
xmin = -4.5*la; xmax = 4.5*la; x = linspace(xmin,xmax,nx);
ymin = -4.5*la; ymax = 4.5*la; y = linspace(ymin,ymax,ny);
zmin = -4.5*la; zmax = 4.5*la; z = linspace(zmin,zmax,nz);
[X,Y,Z] = ndgrid(x,y,z);
u = A(1)*exp(1i*k*(c1(1)*X + c1(2)*Y + c1(3)*Z)) +...
    A(2)*exp(1i*k*(c2(1)*X + c2(2)*Y + c2(3)*Z)) +...
    A(3)*exp(1i*k*(c3(1)*X + c3(2)*Y + c3(3)*Z));
ux = A(1)*1i*k*c1(1)*exp(1i*k*(c1(1)*X + c1(2)*Y + c1(3)*Z)) +...
    A(2)*1i*k*c2(1)*exp(1i*k*(c2(1)*X + c2(2)*Y + c2(3)*Z)) +...
    A(3)*1i*k*c3(1)*exp(1i*k*(c3(1)*X + c3(2)*Y + c3(3)*Z));
uy = A(1)*1i*k*c1(2)*exp(1i*k*(c1(1)*X + c1(2)*Y + c1(3)*Z)) +...
    A(2)*1i*k*c2(2)*exp(1i*k*(c2(1)*X + c2(2)*Y + c2(3)*Z)) +...
    A(3)*1i*k*c3(2)*exp(1i*k*(c3(1)*X + c3(2)*Y + c3(3)*Z));
uz = A(1)*1i*k*c1(3)*exp(1i*k*(c1(1)*X + c1(2)*Y + c1(3)*Z)) +...
    A(2)*1i*k*c2(3)*exp(1i*k*(c2(1)*X + c2(2)*Y + c2(3)*Z)) +...
    A(3)*1i*k*c3(3)*exp(1i*k*(c3(1)*X + c3(2)*Y + c3(3)*Z));

%calculate acoustic radiation potential
rp = a*real(u).^2 - b*( real(ux).^2 + real(uy).^2 + real(uz).^2);

figure; hold on;
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

if strcmp(lattice,'cubic-body-centered') ||...
        strcmp(lattice,'tetragonal-body-centered')
    line([-a2(1),a1(1)]/la,[-a2(2),a1(2)]/la,[-a2(3),a1(3)]/la,...
        'Color','k','LineStyle',':');
    line([a1(1),-a3(1)]/la,[a1(2),-a3(2)]/la,[a1(3),-a3(3)]/la,...
        'Color','k','LineStyle',':');
    line([-a3(1),-a1(1)-a2(1)-a3(1)]/la,[-a3(2),-a1(2)-a2(2)-a3(2)]/la,...
        [-a3(3),-a1(3)-a2(3)-a3(3)]/la,'Color','k','LineStyle',':');
    line([-a1(1)-a2(1)-a3(1),-a2(1)]/la,[-a1(2)-a2(2)-a3(2),-a2(2)]/la,...
        [-a1(3)-a2(3)-a3(3),-a2(3)]/la,'Color','k','LineStyle',':');
    line([-a2(1),a3(1)]/la,[-a2(2),a3(2)]/la,[-a2(3),a3(3)]/la,...
        'Color','k','LineStyle',':');
    line([a1(1),a1(1)+a2(1)+a3(1)]/la,[a1(2),a1(2)+a2(2)+a3(2)]/la,...
        [a1(3),a1(3)+a2(3)+a3(3)]/la,'Color','k','LineStyle',':');
    line([-a3(1),a2(1)]/la,[-a3(2),a2(2)]/la,[-a3(3),a2(3)]/la,...
        'Color','k','LineStyle',':');
    line([-a1(1)-a2(1)-a3(1),-a1(1)]/la,[-a1(2)-a2(2)-a3(2),-a1(2)]/la,...
        [-a1(3)-a2(3)-a3(3),-a1(3)]/la,'Color','k','LineStyle',':');
    line([a3(1),a1(1)+a2(1)+a3(1)]/la,[a3(2),a1(2)+a2(2)+a3(2)]/la,...
        [a3(3),a1(3)+a2(3)+a3(3)]/la,'Color','k','LineStyle',':');
    line([a1(1)+a2(1)+a3(1),a2(1)]/la,[a1(2)+a2(2)+a3(2),a2(2)]/la,...
        [a1(3)+a2(3)+a3(3),a2(3)]/la,'Color','k','LineStyle',':');
    line([a2(1),-a1(1)]/la,[a2(2),-a1(2)]/la,[a2(3),-a1(3)]/la,...
        'Color','k','LineStyle',':');
    line([-a1(1),a3(1)]/la,[-a1(2),a3(2)]/la,[-a1(3),a3(3)]/la,...
        'Color','k','LineStyle',':');
end

if strcmp(lattice,'orthorhombic-face-centered') ||...
        strcmp(lattice,'cubic-face-centered')
    line([a1(1)-a2(1)+a3(1),2*a1(1)]/la,[a1(2)-a2(2)+a3(2),2*a1(2)]/la,...
        [a1(3)-a2(3)+a3(3),2*a1(3)]/la,'Color','k','LineStyle',':');
    line([2*a1(1),a1(1)+a2(1)-a3(1)]/la,[2*a1(2),a1(2)+a2(2)-a3(2)]/la,...
        [2*a1(3),a1(3)+a2(3)-a3(3)]/la,'Color','k','LineStyle',':');
    line([a1(1)+a2(1)-a3(1),0]/la,[a1(2)+a2(2)-a3(2),0]/la,...
        [a1(3)+a2(3)-a3(3),0]/la,'Color','k','LineStyle',':');
    line([0,a1(1)-a2(1)+a3(1)]/la,[0,a1(2)-a2(2)+a3(2)]/la,...
        [0,a1(3)-a2(3)+a3(3)]/la,'Color','k','LineStyle',':');
    line([a1(1)-a2(1)+a3(1),2*a3(1)]/la,[a1(2)-a2(2)+a3(2),2*a3(2)]/la,...
        [a1(3)-a2(3)+a3(3),2*a3(3)]/la,'Color','k','LineStyle',':');
    line([2*a1(1),a1(1)+a2(1)+a3(1)]/la,[2*a1(2),a1(2)+a2(2)+a3(2)]/la,...
        [2*a1(3),a1(3)+a2(3)+a3(3)]/la,'Color','k','LineStyle',':');
    line([a1(1)+a2(1)-a3(1),2*a2(1)]/la,[a1(2)+a2(2)-a3(2),2*a2(2)]/la,...
        [a1(3)+a2(3)-a3(3),2*a2(3)]/la,'Color','k','LineStyle',':');
    line([0,-a1(1)+a2(1)+a3(1)]/la,[0,-a1(2)+a2(2)+a3(2)]/la,...
        [0,-a1(3)+a2(3)+a3(3)]/la,'Color','k','LineStyle',':');
    line([2*a3(1),a1(1)+a2(1)+a3(1)]/la,[2*a3(2),a1(2)+a2(2)+a3(2)]/la,...
        [2*a3(3),a1(3)+a2(3)+a3(3)]/la,'Color','k','LineStyle',':');
    line([a1(1)+a2(1)+a3(1),2*a2(1)]/la,[a1(2)+a2(2)+a3(2),2*a2(2)]/la,...
        [a1(3)+a2(3)+a3(3),2*a2(3)]/la,'Color','k','LineStyle',':');
    line([2*a2(1),-a1(1)+a2(1)+a3(1)]/la,...
        [2*a2(2),-a1(2)+a2(2)+a3(2)]/la,[2*a2(3),-a1(3)+a2(3)+a3(3)]/la,...
        'Color','k','LineStyle',':');
    line([-a1(1)+a2(1)+a3(1),2*a3(1)]/la,...
        [-a1(2)+a2(2)+a3(2),2*a3(2)]/la,[-a1(3)+a2(3)+a3(3),2*a3(3)]/la,...
        'Color','k','LineStyle',':');
end

% Find points in the unit cell
xyz = [0,0,0; a1'; a2'; a3'; a1'+a2'; a1'+a3'; a2'+a3'; a1'+a2'+a3'];
ctr_pt = sum(xyz)/8;
xyz = xyz - ones(8,1)*ctr_pt;
xyz = 1.3*xyz;
xyz = xyz + ones(8,1)*ctr_pt;
testpts = [X(:),Y(:),Z(:)];
in = inhull(testpts,xyz);
ind3 = reshape(in,size(X)); % indices of points in the unit cell

thresh = l_0+thresh_const*(max(diag(D))-l_0);
ind1 = rp<thresh; % indices of points around minima

ind2 = ind1; % indices of points near other minima
for ii=-5:5
    for jj=-5:5
        for kk=-5:5
            ind2 = ind2 & sqrt((X-(ii*a1(1)+jj*a2(1)+kk*a3(1))).^2 +...
                (Y-(ii*a1(2)+jj*a2(2)+kk*a3(2))).^2 + ...
                (Z-(ii*a1(3)+jj*a2(3)+kk*a3(3))).^2)>dist_pts;
        end
    end
end

light('Position',light_position)

% 
Points1 = [X(ind1&ind3&~ind2),Y(ind1&ind3&~ind2),Z(ind1&ind3&~ind2)];
[~,centroids1] = kmeans(Points1/la,num_blue);
[X1S,X2S,X3S] = sphere(20);
col = 'b';
for i = 1:num_blue
    surf(r*X1S+centroids1(i,1), r*X2S+centroids1(i,2),...
        r*X3S+centroids1(i,3),'facecolor',col,'edgecolor','none',...
        'facelighting','gouraud');
end

col = 'r';
Points2 = [X(ind2&ind3),Y(ind2&ind3),Z(ind2&ind3)];
[~,centroids2] = kmeans(Points2/la,num_red);
for i = 1:num_red
    surf(r*X1S+centroids2(i,1), r*X2S+centroids2(i,2),...
        r*X3S+centroids2(i,3),'facecolor',col,'edgecolor','none',...
        'facelighting','gouraud');
end

axis equal;

axis off;
view(azim,elev);
if strcmp(lattice,'cubic-face-centered')
    zlim([-0.5,2]);
end
if strcmp(lattice,'cubic-primitive')
    zlim([-0.25,1.15]);
end
if strcmp(lattice,'orthorhombic-face-centered')
    zlim([-3.5,4]);
end
arrow3(start,start+c1/arrow_scale,'*',arrow_width,arrow_length);
arrow3(start,start+c2/arrow_scale,'*',arrow_width,arrow_length);
arrow3(start,start+c3/arrow_scale,'*',arrow_width,arrow_length);
filename = [C{1,j},'.png'];
print('-dpng','-r300',filename);
system(['/usr/local/bin/mogrify -trim -define png:include-chunk=none ',...
    filename])

end