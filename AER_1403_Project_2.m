clear;
clc;
close all;

% ele_size = [10,5,4,2.5,1.25,1,0.625,0.5]*10^-3;
% for i = 1:length(ele_size)
%     [F{i},d{i},elem_num(i),von_mises{i},cnvg_d(i)] = FEA(ele_size(i),0);
% end
% figure;
% plot(elem_num,cnvg_d,'-','Color','r');
% grid on;
% grid minor;
% ax = gca;
% ax.GridLineStyle = '-';
% ax.GridAlpha = 0.2;
% ax.MinorGridLineStyle ='-';
% ax.MinorGridAlpha = 0.2;
% xlabel('Number of Elements');
% ylabel('Y Displacement at Top Right Corner');

[F,d,elem_num,von_mises,cnvg_d] = FEA(2*10^-3,1);

% s = 2*10^-3;
% figureOption = 1;

function [F,d,elem_num,von_mises,cnvg_d] = FEA(s,figureOption)
tic

% Dimensions
x_1 = 80*10^-3; % m
y_1 = 40*10^-3; % m
x_2 = 20*10^-3; % m
y_2 = 100*10^-3; % m

% Number of Elements on each side
x1 = x_1/s;
y1 = y_1/s;
x2 = x_2/s;
y2 = y_2/s;
elem_num = x1*y1 + x2*y2;

% Loading Conditions
P = 5000/2; % N

% Stiffness, Force, and Displacement Matrix Preallocation
K = sparse(2*(y1+1)*(x1+1) + 2*(y2+1)*(x2+1) - 2*(y1+1),2*(y1+1)*(x1+1) + 2*(y2+1)*(x2+1) - 2*(y1+1));
F = sparse(2*(y1+1)*(x1+1) + 2*(y2+1)*(x2+1) - 2*(y1+1),1);
d = zeros(2*(y1+1)*(x1+1) + 2*(y2+1)*(x2+1) - 2*(y1+1),1);

% Strain x,y,xy Matrix Preallocation
epsilon_x = zeros((y1+1)*(x1+1) + (y2+1)*(x2+1) - (y1+1),1);
epsilon_y = zeros((y1+1)*(x1+1) + (y2+1)*(x2+1) - (y1+1),1);
epsilon_xy = zeros((y1+1)*(x1+1) + (y2+1)*(x2+1) - (y1+1),1);

% Stress x,y,xy Matrix Preallocation
sigma_x = zeros((y1+1)*(x1+1) + (y2+1)*(x2+1) - (y1+1),1);
sigma_y = zeros((y1+1)*(x1+1) + (y2+1)*(x2+1) - (y1+1),1);
sigma_xy = zeros((y1+1)*(x1+1) + (y2+1)*(x2+1) - (y1+1),1);
% von_mises = zeros((y1+1)*(x1+1) + (y2+1)*(x2+1) - (y1+1),1);

% Element Stiffness Matrix
[K_e,D_e,H_e] = stiffCalc(s);

% Coordinate of Each node
xc = zeros(1,(y1+1)*(x1+1) + (y2+1)*(x2+1) - (y1+1));
yc = zeros(1,(y1+1)*(x1+1) + (y2+1)*(x2+1) - (y1+1));

% Part 1 Stiffness Matrix, Node, DOF
counter = 0;
for i = 1:x1
    for j = 1:y1
        counter = counter + 1;
        node = [(y1+1)*(i-1)+j, (y1+1)*(i-1)+j+1, (y1+1)*i+j+1, (y1+1)*i+j];
        node_M2{counter} = node;
        dof = [2*((y1+1)*(i-1)+j)-1, 2*((y1+1)*(i-1)+j), 2*((y1+1)*(i-1)+j+1)-1, 2*((y1+1)*(i-1)+j+1), 2*((y1+1)*i+j+1)-1, 2*((y1+1)*i+j+1), 2*((y1+1)*i+j)-1, 2*((y1+1)*i+j)];
        dof_M{counter} = dof;
        K(dof,dof) = K(dof,dof) + K_e;
        xc(node) = [(i-1)*s,(i-1)*s,i*s,i*s];
        yc(node) = 0.1 - [(j-1)*s,j*s,j*s,(j-1)*s];
        if i == x1 && j == 1
            dof_num_1 = max(dof) - 4;
            node_num_1 = max(node) - 2;
        end
    end
end

% Part 2 Stiffness Matrix, Node, DOF
for i = 1:x2
    for j = 1:y2
        counter = counter + 1;
        node = node_num_1 + [(y2+1)*(i-1)+j, (y2+1)*(i-1)+j+1, (y2+1)*i+j+1, (y2+1)*i+j];       
        node_M2{counter} = node;
        dof = dof_num_1 + [2*((y2+1)*(i-1)+j)-1, 2*((y2+1)*(i-1)+j), 2*((y2+1)*(i-1)+j+1)-1, 2*((y2+1)*(i-1)+j+1), 2*((y2+1)*i+j+1)-1, 2*((y2+1)*i+j+1), 2*((y2+1)*i+j)-1, 2*((y2+1)*i+j)];
        dof_M{counter} = dof;
        K(dof,dof) = K(dof,dof) + K_e;
        xc(node) = x1*s + [(i-1)*s,(i-1)*s,i*s,i*s];
        yc(node) = 0.1 - [(j-1)*s,j*s,j*s,(j-1)*s];
        if i == x2 && j == 1
            cnvg_node = max(node) - 1;
        end
    end
end

% Force and Displacement Constraint
node_F = -P/(x2 + 1 + x2 + 1 - 2);
F(dof_num_1 + 2*(y2+1):2*(y2+1):dof_num_1 + 2*(x2+1)*(y2+1)) = 2*node_F;
F(dof_num_1 + 2*(x2+1)*(y2+1)) = node_F;
F(dof_num_1 + 2*(y2+1)) = node_F;

% Setting Constrained, All, and Normal Degrees of Freedom
dof_rstn = [1:1:2*(y1+1)];
dof_rstn_2 = [2*(y1+1)*(x1+1) + 2*(y2+1)*x2 - 2*(y1+1)+1:2:2*(y1+1)*(x1+1) + 2*(y2+1)*(x2+1) - 2*(y1+1)];
dof_all = [1:1:2*(y1+1)*(x1+1) + 2*(y2+1)*(x2+1) - 2*(y1+1)];
dof_norm = setdiff(dof_all,dof_rstn,'stable');
dof_norm = setdiff(dof_norm,dof_rstn_2,'stable');

% Solving for Displacement of Eahc DOF
d(dof_norm) = K(dof_norm,dof_norm) \ F(dof_norm);

% Rearrange the format of Node and DOF Matrices 
node_M2 = cell2mat(node_M2);
node_M2 = reshape(node_M2,4,[]);
dof_M = cell2mat(dof_M);
dof_M = reshape(dof_M,8,[]);

% Caculate Strain and Stress
counter = 0;
for i = 1:length(node_M2)
    node = node_M2(:,i);
%     xcc = xc(node);
%     ycc = yc(node);

    for j = 1:length(node)
        counter = counter + 1;
        epsilon = H_e{j}*d(dof_M(:,i));
        epsilon_M{counter} = epsilon;
        epsilon_x(node(j)) = epsilon(1);
        epsilon_y(node(j)) = epsilon(2);
        epsilon_xy(node(j)) = epsilon(3);

        sigma = D_e*epsilon;
        sigma_M{counter} = sigma;
        sigma_x(node(j)) = sigma(1);
        sigma_y(node(j)) = sigma(2);
        sigma_xy(node(j)) = sigma(3);
    end
end

epsilon_M = cell2mat(epsilon_M);
sigma_M = cell2mat(sigma_M);

epsilon_x_avg = reshape(epsilon_M(1,:),4,[]);
epsilon_x_avg = mean(epsilon_x_avg);
epsilon_y_avg = reshape(epsilon_M(2,:),4,[]);
epsilon_y_avg = mean(epsilon_y_avg);
epsilon_xy_avg = reshape(epsilon_M(3,:),4,[]);
epsilon_xy_avg = mean(epsilon_xy_avg);

sigma_x_avg = reshape(sigma_M(1,:),4,[]);
sigma_x_avg = mean(sigma_x_avg);
sigma_y_avg = reshape(sigma_M(2,:),4,[]);
sigma_y_avg = mean(sigma_y_avg);
sigma_xy_avg = reshape(sigma_M(3,:),4,[]);
sigma_xy_avg = mean(sigma_xy_avg);

for i = 1:length(sigma_x_avg)
    von_mises(i) = sqrt(((sigma_x_avg(i) - sigma_y_avg(i))^2 + sigma_y_avg(i)^2 + sigma_x_avg(i)^2 + 6*sigma_xy_avg(i)^2)/2);
end

for i = 1:length(d)
    if rem(i,2) == 1
        dx(ceil(i/2)) = d(i);
    else
        dy(i/2) = d(i);
    end
end

xc_d = xc + dx;
yc_d = yc + dy;
dd = sqrt(dx.^2 + dy.^2);
adjust_fac = 10;
xc_daf = xc + dx*adjust_fac;
yc_daf = yc + dy*adjust_fac;

% Plot the Graphs
if figureOption == 1

    % Original Shape Plot and After Shape Plot
    figure;
    plot(xc(node_M2),yc(node_M2),'-','Color','b');
    hold on;
    plot(xc(1:(y1+1):(y1+1)*(x1+1) - y1),yc(1:(y1+1):(y1+1)*(x1+1) - y1),'-','Color','b');
    hold on;
    plot(xc((y1+1)*(x1+1) - y1:(y2+1):(y1+1)*(x1+1) + (y2+1)*(x2+1) - (y1+1) - y2),yc((y1+1)*(x1+1) - y1:(y2+1):(y1+1)*(x1+1) + (y2+1)*(x2+1) - (y1+1) - y2),'-','Color','b');
    hold on;
    plot(xc_d(node_M2),yc_d(node_M2),'-','Color','r');
    hold on;
    plot(xc_d(1:(y1+1):(y1+1)*(x1+1) - y1),yc_d(1:(y1+1):(y1+1)*(x1+1) - y1),'-','Color','r');
    hold on;
    plot(xc_d((y1+1)*(x1+1) - y1:(y2+1):(y1+1)*(x1+1) + (y2+1)*(x2+1) - (y1+1) - y2),yc_d((y1+1)*(x1+1) - y1:(y2+1):(y1+1)*(x1+1) + (y2+1)*(x2+1) - (y1+1) - y2),'-','Color','r');
    axis equal;
    
    % adjusting factor added
    figure;
    plot(xc(node_M2),yc(node_M2),'-','Color','b');
    hold on;
    plot(xc(1:(y1+1):(y1+1)*(x1+1) - y1),yc(1:(y1+1):(y1+1)*(x1+1) - y1),'-','Color','b');
    hold on;
    plot(xc((y1+1)*(x1+1) - y1:(y2+1):(y1+1)*(x1+1) + (y2+1)*(x2+1) - (y1+1) - y2),yc((y1+1)*(x1+1) - y1:(y2+1):(y1+1)*(x1+1) + (y2+1)*(x2+1) - (y1+1) - y2),'-','Color','b');
    hold on;
    plot(xc_daf(node_M2),yc_daf(node_M2),'-','Color','r');
    hold on;
    plot(xc_daf(1:(y1+1):(y1+1)*(x1+1) - y1),yc_daf(1:(y1+1):(y1+1)*(x1+1) - y1),'-','Color','r');
    hold on;
    plot(xc_daf((y1+1)*(x1+1) - y1:(y2+1):(y1+1)*(x1+1) + (y2+1)*(x2+1) - (y1+1) - y2),yc_daf((y1+1)*(x1+1) - y1:(y2+1):(y1+1)*(x1+1) + (y2+1)*(x2+1) - (y1+1) - y2),'-','Color','r');
    axis equal;
    
    figure;
    patch(xc(node_M2),yc(node_M2),dx(node_M2));
    colorbar;
    title('x Displacement Distribution');
    colormap jet
    axis equal;
    
    figure;
    patch(xc(node_M2),yc(node_M2),dy(node_M2));
    colorbar;
    title('y Displacement Distribution');
    colormap jet
    axis equal;
    
    figure;
    patch(xc(node_M2),yc(node_M2),dd(node_M2));
    colorbar;
    title('Total Displacement Distribution');
    colormap jet
    axis equal;

    figure;
    patch(xc(node_M2),yc(node_M2),epsilon_x_avg);
    title('x Strain Distribution');
    colorbar;
    colormap jet
    axis equal;
    
    figure;
    patch(xc(node_M2),yc(node_M2),epsilon_y_avg);
    title('y Strain Distribution');
    colorbar;
    colormap jet
    axis equal;
    
    figure;
    patch(xc(node_M2),yc(node_M2),epsilon_xy_avg);
    title('xy Strain Distribution');
    colorbar;
    colormap jet
    axis equal;
    
    figure;
    patch(xc(node_M2),yc(node_M2),sigma_x_avg);
    title('x stress Distribution');
    colorbar;
    colormap jet
    axis equal;
    
    figure;
    patch(xc(node_M2),yc(node_M2),sigma_y_avg);
    title('y stress Distribution');
    colorbar;
    colormap jet
    axis equal;
    
    figure;
    patch(xc(node_M2),yc(node_M2),sigma_xy_avg);
    title('xy stress Distribution');
    colorbar;
    colormap jet
    axis equal;
    
    figure;
    patch(xc(node_M2),yc(node_M2),von_mises);
    title('Von Mises Stress Distribution');
    colorbar;
    colormap jet
    axis equal;
    
end

% Convergence Testing
cnvg_d = dy(cnvg_node);

toc
end


function [K,D,H] = stiffCalc(s)
E = 70*10^9; % Pa
nu = 0.3;
t = 2.5*10^-3; % m

D = (E/(1-nu^2))*[1,nu,0;nu,1,0;0,0,(1-nu)/2];

xc = [0,0,s,s].';
yc = [0,-s,-s,0].';

xyc = [xc, yc];

xiM = [-1/sqrt(3),-1/sqrt(3),1/sqrt(3),1/sqrt(3)];
zetaM = [1/sqrt(3),-1/sqrt(3),-1/sqrt(3),1/sqrt(3)];
w_xi = [1,1,1,1];
w_zeta = [1,1,1,1];

K = 0;
for i = 1:length(xiM)
%     dN1_dxi = (1/4)*(zetaM(i) - 1);
%     dN1_dzeta = (1/4)*(xiM(i) - 1);
%     dN2_dxi = (1/4)*(1 - zetaM(i));
%     dN2_dzeta = (1/4)*(-1 - xiM(i));
%     dN3_dxi = (1/4)*(1 + zetaM(i));
%     dN3_dzeta = (1/4)*(1 + xiM(i));
%     dN4_dxi = (1/4)*(-1 - zetaM(i));
%     dN4_dzeta = (1/4)*(1 - xiM(i));

    dN1_dxi = (1/4)*(-1 - zetaM(i));
    dN1_dzeta = (1/4)*(1 - xiM(i));
    dN2_dxi = (1/4)*(zetaM(i) - 1);
    dN2_dzeta = (1/4)*(xiM(i) - 1);
    dN3_dxi = (1/4)*(1 - zetaM(i));
    dN3_dzeta = (1/4)*(-1 - xiM(i));
    dN4_dxi = (1/4)*(1 + zetaM(i));
    dN4_dzeta = (1/4)*(1 + xiM(i));
    
    GN_4Q = [dN1_dxi,dN2_dxi,dN3_dxi,dN4_dxi;
        dN1_dzeta,dN2_dzeta,dN3_dzeta,dN4_dzeta];
    
    J = GN_4Q*xyc;
    detJ = det(J);
    Hstar = J^-1*GN_4Q;
    H{i} = [Hstar(1,1), 0, Hstar(1,2), 0, Hstar(1,3), 0, Hstar(1,4), 0;
        0, Hstar(2,1), 0, Hstar(2,2), 0, Hstar(2,3), 0, Hstar(2,4);
        Hstar(2,1),Hstar(1,1),Hstar(2,2),Hstar(1,2),Hstar(2,3),Hstar(1,3),Hstar(2,4),Hstar(1,4)];
    K = K + w_xi(i)*w_zeta(i)*detJ*H{i}.'*D*H{i};
end
K = double(K)*t;
end
