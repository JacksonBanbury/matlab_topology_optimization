% Topology Optimization script written by Jackson Banbury for AER1410 - Assignment 2
% November, 2017
% Based on Ole Sigmund's 99-line Topology Optimization code
% When running the code, choose 0 or 1 to plot without or with deformation
function top(deformation);
% Initializing
% loads provided data file
load ps_02_data.mat % ensure this is placed in the same working directory
E = 70e9; % Young's Modulus provided as 70 GPa
sigma = 325e6; % Yield strength provided as 325 MPa
nu = 0.3; % Poisson's Ratio
vfrac = 0.3; % volume fraction
penal = 3; % penalization factor
nelx = 50; % Number of elements in the x-direction
nely = 50; % Number of elements in the y-direction
volume = 1e-9; % Volume of a single element (1mm x 1mm x 1mm)

% Reforming the displacement matrix
dx_flat = reshape(dx.',1,[]); % flattens the displacement matrices into vectors
dy_flat = reshape(dy.',1,[]);
U = [dx_flat; dy_flat];
U = U(:)'; % This combines the flattened displacement matrices into the form seen in Sigmund's code (alternating x and y)
U = U.'.*(1/1000); % Converts the displacement from mm (as provided) into meters

if (deformation == 0)
%  Plotting Undeformed Mesh
    colormap(gray); imagesc(-rho); axis equal; axis tight; pause(1e-6);
    hold on
    quiver(dx,dy) % superimposes the nodal displacements as vectors over the undeformed mesh (scale = 1)
    hold off
else
% Separate function to plot deformed mesh
    plot_deformed(nelx,nely,U,rho);
    
end

[stresses_y, stresses_x, stresses_z]=stress(nelx,nely,U,rho,nu,E,sigma,penal);

[strain_energy_y,strain_energy_x,strain_energy_z] = strain_energy(stresses_y,stresses_x,stresses_z,E,volume);

compliance(nely,nelx,rho,U,E,nu,penal);


end

function plot_deformed(nelx,nely, U,rho);
% Plotting Deformed Mesh

colormap(gray); axis equal; % The below code is pulled from Sigmund's Topology course page: http://www.topopt.dtu.dk/DCAMM/
for ely = 1:nely
   for elx = 1:nelx
       n1 = (nely+1)*(elx-1)+ely;
       n2 = (nely+1)*elx +ely;
       Ue = 30000*U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1); 
       ly = ely-1; lx = elx-1; 
       xcoord = [Ue(1,1)+lx Ue(3,1)+lx+1 Ue(5,1)+lx+1 Ue(7,1)+lx ]';
       ycoord = [-Ue(2,1)-ly -Ue(4,1)-ly -Ue(6,1)-ly-1 -Ue(8,1)-ly-1]'; 
       patch(xcoord,ycoord,-rho(ely,elx)) 
    end 
end 
drawnow;

end

function [stresses_y, stresses_x, stresses_z]=stress(nelx,nely,U,rho,nu,E,sigma,penal);
% Finds the stress at the Gauss points

area = 1e-6; % assuming each element is 1m x 1mm (as provided in the ps)

% initializing variables
K = zeros(8); % total stiffness matrix
K_star = zeros(8); %current stiffness matrix
d = zeros(8,1); % nodal displacement vector

stresses_y = zeros(nelx,nely,2,2);
stresses_x = zeros(nelx,nely,2,2);
stresses_z = zeros(nelx,nely,2,2);

D = (E/(1-nu^2)) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu) / 2];

% Gauss points:
xi(1) = -1/sqrt(3); xi(2) = 1/sqrt(3);
eta(1) = -1/sqrt(3); eta(2) = 1/sqrt(3);
w(1) = 1; w(2) = 1; % weights
J_det = area/4; % determinant of the jacobian (Area/4)

% The following is for the four Gauss points in a *single* element

for ely = 1:nely % incrementing element rows
    for elx = 1:nelx % incrementing element columns
        % Below this is for a single element
        x_e = [elx elx elx+1 elx+1];
        y_e = [ely ely+1 ely ely+1];
 
        for i = 1:2
            for j = 1:2
                x = 1+(0.5 * xi(i) + 0.5);
                y = 3*(0.5 * eta(j) + 0.5);

                n1 = (nely+1)*(elx-1)+ely; 
                n2 = (nely+1)* elx +ely; 
                d = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1); % Displacement at Gauss points

                H = (1 / area) * [(ely - y_e(4)), 0, -(ely - y_e(4)), 0, (ely - y_e(1)), 0 , -(ely - y_e(1)), 0; 
                    0, (elx - x_e(2)), 0, -(elx - x_e(1)), 0, (elx - x_e(1)), 0, -(elx - x_e(2)); 
                    (elx - x_e(2)), (ely - y_e(4)), -(elx - x_e(1)), -(ely - y_e(4)), (elx - x_e(1)), (ely - y_e(1)), -(elx - x_e(2)), -(ely - y_e(1))];

                stresses = D * H * d; % the von mises stress for this particular element
                %disp(strain)

                stresses_gauss_y(elx,ely,i,j) = stresses(1,:);
                stresses_gauss_x(elx,ely,i,j) = stresses(2,:);
                stresses_gauss_z(elx,ely,i,j) = stresses(3,:);
            end
        end
    end
end

% Take the mean of the stresses at the gauss points to get a suitable value
% to represent the stress of the element as a whole
stresses_y = mean(mean(stresses_gauss_y,4),3);
stresses_x = mean(mean(stresses_gauss_x,4),3);
stresses_z = mean(mean(stresses_gauss_z,4),3);

% The per-element von mises stress
von_mises_stress = (((stresses_x-stresses_y).^2 + (stresses_y-stresses_z).^2 + (stresses_z-stresses_x).^2 + (stresses_x.^2+stresses_y.^2+stresses_z.^2).*6).*0.5).^(1/2);

% Prints the stress values for each element into a 50x50 CSV file
Ty = array2table(stresses_y);
writetable(Ty,'stresses_y.csv');
Tx = array2table(stresses_x);
writetable(Tx,'stresses_x.csv');
Tz = array2table(stresses_z);
writetable(Tz,'stresses_z.csv');

Tvon = array2table(von_mises_stress);
writetable(Tvon,'von_mises_stresses.csv');

% The below code plots the contour map of the von Mises stress
%  countour_von = contourf(von_mises_stress) % Used this to obtain contour plots
%  colormap(hot)
%  set(gca, 'ydir', 'reverse');
p_norm(nely,nelx,von_mises_stress,rho,sigma,penal);

end

function [strain_energy_y,strain_energy_x,strain_energy_z]=strain_energy(stresses_y,stresses_x,stresses_z,E,volume);

% The strain energy for each element is given by U = 1/2 * V/E * stress^2
strain_energy_y = 0.5*(volume/E)*(stresses_y.^2);
strain_energy_x = 0.5*(volume/E)*(stresses_x.^2);
strain_energy_z = 0.5*(volume/E)*(stresses_z.^2);

% Prints the strain energy values for each element into a 50x50 CSV file
Ty = array2table(strain_energy_y);
writetable(Ty,'strain_energy_y.csv');
Tx = array2table(strain_energy_x);
writetable(Tx,'strain_energy_x.csv');
Tz = array2table(strain_energy_z);
writetable(Tz,'strain_energy_z.csv');

% Total strain energy
total_strain_energy_y = sum(sum(strain_energy_y)); % sums twice because first only sums rows
total_strain_energy_x = sum(sum(strain_energy_x));
total_strain_energy_z = sum(sum(strain_energy_z));

disp("Total strain energy in Y-direction:");
disp(total_strain_energy_y);
disp("Total strain energy in X-direction:");
disp(total_strain_energy_z);
disp("Total strain energy in Z-direction:");
disp(total_strain_energy_z);

end

function compliance(nely,nelx,rho,U,E,nu,penal);
[KE] = lk(E,nu);
c_total = 0.;
for ely = 1:nely
   for elx = 1:nelx
     n1 = (nely+1)*(elx-1)+ely; 
     n2 = (nely+1)* elx   +ely;
     Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
     c_el = rho(ely,elx)^penal*Ue'*KE*Ue; % Compliance for this element
     c_total = c_total + c_el; % The total compliance
     dc(ely,elx) = -penal*rho(ely,elx)^(penal-1)*Ue'*KE*Ue; % This is the gradient of the compliance for the current element
   end
end

Tdc = array2table(dc);
writetable(Tdc,'grad_compliance.csv');

% Contour plot for the gradient of the compliance
contour_dc = contourf(dc)
colormap(hot)
set(gca, 'ydir', 'reverse');
end

function [KE]=lk(E,nu);
k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ... 
   -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];
KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
                  k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
                  k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
                  k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
                  k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
                  k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
                  k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
                  k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
end

function p_norm(nely,nelx,von_mises_stress,rho,sigma,penal);
i=1;
boundary = 0 % Set this to zero to observe P-norm for q4, set it to 1 to observe p-norm while disregarding elements on the boundaries i.e. q5
for P = 1:.1:5.4
    p = 0;
    for ely = 1+boundary:nely-boundary
        for elx = 1+boundary:nelx-boundary
            p_inner = von_mises_stress(ely,elx)/(rho(ely,elx).^penal * sigma);% + 1000 - 1000/rho(ely,elx);
            p = p + max(0,(p_inner))^P;
        end
    end
   p_norm(i) = p^(1/P)
   P_index(i) = P;
   i=i+1;
end
% 
% scatter(P_index,p_norm,'filled');
% xlabel('P value');
% ylabel('P-norm');

end
