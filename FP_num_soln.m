%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FP numerical simulation of corrugated nanopores using FDM
% Course project : Transport Phenomena I
% Contributors   : Sthavishtha
% works in MATLAB R2017b version and above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;

% Simulation parameters
XMIN = -1.5; DX = 0.01; XMAX = +1.5;
YMIN = -2.1; DY = 0.01; YMAX = 2.1;
NTIME = 5000; DT = 0.0001;
NX = (XMAX - XMIN)/DX + 1;
NY = (YMAX - YMIN)/DY + 1;
Xedges = XMIN:DX:XMAX;
Yedges = YMIN:DY:YMAX;

% grid-size related parameters
inv_DX = 1./(DX);
inv_DY = 1./(DY);
inv_DX2 = 1./(DX*DX);
inv_DY2 = 1./(DY*DY);

% Top & Bottom wall parameters
% wall offset
% amp = (1/2)/pi;
amp = 1;
% wall amplitude
% off = (1.02/2)/pi;
off = 1.05;

syms B_u(t) B_d(t)
% top wall variation
B_u(t) = amp*cos(2*t*pi) + off;
% B_u(t) = amp + off;
% bottom wall variation
B_d(t) = -B_u(t);
% derivatives of wall variations
diff_B_u = diff(B_u, t);
diff_B_d = diff(B_d, t);

% Force along x-direction
fx_mag = 1;

lin_x = Xedges;
lin_y = Yedges;

% initial prob density distributions
p_old = zeros(NX, NY);
p_old((NX + 1)/2, (NY + 1)/2) = 1.; % delta distribution

% absorbing bcs (left & right walls)
p_old([1 NX], :) = 0;
% prob at a new time step
% p_new = p_old;

% plotting - for simple visualization
% figure(1);
% [X, Y] = meshgrid(1:NX, 1:NY);
% h = surf(X, Y, p_old');

% diffusion coeff extracted from BD simulations 
% Deffx = 0.325;
Deffy = 0.6119;
Deffx = 0.1375;
% Deffx = 0.343;
% Deffy = 0.833;
% Deffx = 0.26625;
% Deffy = 0.7494;

% error tolerance and iteration number
tol = 0.;
iter = 0;

% prob summed up along y direction
sum_px = zeros(NX, 1);

for i = 1:NX
    for j = 1:NY
        sum_px(i) = p_old(i, j) + sum_px(i);
    end
end

% initial delta prob dist plot
figure(2);
hold on;
plot(lin_x, B_u(lin_x), 'black');
plot(lin_x, B_d(lin_x), 'black');
plot(lin_x, sum_px);
box on;
hold off;

outside_old = zeros(NX, NY); % 0 : inside, -1 : on bdary, 1 : outside bdary

% flagging inside/outside the bdary domains
for i = 1:NX
    for j = 1:NY
        
        x = XMIN + (i - 1)*DX;
        y = YMIN + (j - 1)*DY;
        
        if amp*cos(2*x*pi) + off < abs(y) % true for outside the domain
            outside_old(i, j) = 1;
        end
    end
end

outside = outside_old;

% flagging the bdary/near bdary nodes
for i = 1:NX
    for j = 1:NY
        
        if outside_old(i, j) == 0
            if i < NX
                if outside_old(i + 1, j) == 1
                    outside(i, j) = -1;
                end
            end
            
            if i > 1
                if outside_old(i - 1, j) == 1
                    outside(i, j) = -1;
                end
            end
            
            if j < NY
                if outside_old(i, j + 1) == 1
                    outside(i, j) = -1;
                end
            end
            
            if j > 1
                if outside_old(i, j - 1) == 1
                    outside(i, j) = -1;
                end
            end
        end
    end
end

% clear distinctive visualization of outside/inside nodes
figure(3);
imagesc(outside');

p = zeros(NX + 2, NY);
p(NX + 2, :) = p_old(1, :); 
p(1, :) = p_old(NX, :);
p(2:NX + 1, :) = p_old;
p_new = p;

outside_new = zeros(NX + 2, NY);
outside_new(2:NX + 1, :) = outside;

tic
% looping across the domain
while iter < NTIME
    
    %     p(1, :) = p(NX, :); %periodic bcs
    p([2 NX + 1], :) = 0;
    
    for i = 2:NX + 1
        for j = 2:NY - 1
            
            % on the bdary
            if outside_new(i, j) == -1 
                if outside_new(i + 1, j) == 1 % outside the domain (right side)
%                     p(i + 1, j) = -p(i - 1, j)*(fx_mag + inv_DX)/(-inv_DX + fx_mag);
                      p(i + 1, j) = p(i - 1, j);
                end
                    
                if outside_new(i - 1, j) == 1 % outside the domain (left side)
%                     p(i - 1, j) = -p(i + 1, j)*(-inv_DX + fx_mag)/(fx_mag + inv_DX);   
                      p(i - 1, j) = p(i + 1, j);   
                end
                
                if outside_new(i, j + 1) == 1 % outside the domain (top side)
                    p(i, j + 1) = p(i, j - 1);
                end
                
                if outside_new(i, j - 1) == 1 % outside the domain (bottom side)
                    p(i, j - 1) = p(i, j + 1);
                end
            end
                            
            if outside_new(i, j) == 0 || outside_new(i, j) == -1 % in the nanopore comp domain
                p_new(i, j) = p(i, j) - fx_mag*DT*inv_DX*0.5*(p(i + 1, j) - p(i - 1, j)) + DT*inv_DX2*Deffx*( ...
                    p(i + 1, j) - 2*p(i, j) + p(i - 1, j)) + DT*inv_DY2*Deffy*( ...
                    p(i, j + 1) - 2*p(i, j) + p(i, j - 1));
                
%                 p_new(i, j) = (1. - Deffx*DT*inv_DX2 - Deffy*DT*inv_DY2 + fx_mag*inv_DX*DT)*p_old(i, j) ...
%                     + p_old(i + 1, j)*DT*(0.5*Deffx*inv_DX2 - fx_mag*inv_DX) + p_old(i - 1, j)*DT*(0.5*Deffx*inv_DX2) ...
%                     + (p_old(i, j + 1) + p_old(i, j - 1))*0.5*DT*Deffy*inv_DY2;
                
%             elseif up == 0
%                 % no-flux bc at top wall
%                 %                 p_new(i, j) = (p_old(i, j - 1)*inv_DY - amp*2.*pi*sin(2*x*pi)*0.5*inv_DX*...
%                 %                     (p_old(i + 1, j) - p_old(i - 1, j)))/(inv_DY - fx_mag*amp*2.*pi*sin(2*x*pi));
%                 p_new(i, j) = (-p_old(i, j - 1)*inv_DY - amp*2.*pi*sin(2*x*pi)*inv_DX*p_old(i + 1, j))/(inv_DY - (fx_mag + inv_DX)*amp*2.*pi*sin(2*x*pi));
                
                
%             else
%                 % no-flux bc at bottom wall
%                 %                 p_new(i, j) = (-p_old(i, j + 1)*inv_DY + amp*2.*pi*sin(2*x*pi)*0.5*inv_DX*...
%                 %                     (p_old(i + 1, j) - p_old(i - 1, j)))/(-inv_DY + fx_mag*amp*2.*pi*sin(2*x*pi));
%                 p_new(i, j) = (p_old(i, j + 1)*inv_DY + amp*2.*pi*sin(2*x*pi)*inv_DX*p_old(i + 1, j))/(-inv_DY + (fx_mag + inv_DX)*amp*2.*pi*sin(2*x*pi));
                
            end
            
            tol = tol + abs(p_new(i, j) - p(i, j));
            p(i, j) = p_new(i, j);
            
        end
    end
%     p = p_new;
    
    if mod(iter, 100) == 0
        iter
        tol        
    end
    
    iter = iter + 1;
    
    tol = 0.;
    
end
toc

p_old = p(2:NX + 1, :);

set(groot,'defaultAxesTickLabelInterpreter','latex');   

% Steady state prob surf distribution
% final prob dist surf plot
% figure(4);
% [X, Y] = meshgrid(1:NX, 1:NY);
% h = surf(X, Y, p_old');

sum_px = zeros(NX, 1);
no = zeros(NX, 1);

for i = 1:NX
    for j = 1:NY
        sum_px(i) = p_old(i, j) + sum_px(i);        
    end
    no(i) = size(find(p_old(i, :)), 2);
end

figure(5);
hold on;
Tname = title('X Probability distribution');
xname = xlabel('X direction');
yname = ylabel('Probability density');
set([Tname, xname, yname], 'interpreter', 'latex');
plot(lin_x, sum_px/(2.2*sum(sum_px)/NX), 'r-', 'LineWidth',2); % normalized x probability distriburion 
box on;
hold off;

sum_py = zeros(NY, 1);

for j = 1:NY
    for i = 1:NX
        sum_py(j) = p_old(i, j) + sum_py(j);
    end
end

fname = strcat('FP_t_', num2str(NTIME*DT), '_ntime_', num2str(NTIME), '_dt_', num2str(DT), '_amp_', num2str(amp), '_off_', num2str(off), '_f_', num2str(fx_mag));
save(strcat(fname, '.mat'));

figure(6);
Tname = title('Y Probability distribution');
xname = xlabel('Y direction');
yname = ylabel('Probability density');
set([Tname, xname, yname], 'interpreter', 'latex');
hold on;
plot(lin_y, sum_py/(4*sum(sum_py)/NY), 'r-', 'LineWidth',2); % normalized y probability distriburion 
box on;
hold off;

% figure;
% hold on;
% plot(lin_y, p((NX + 1)/2, :)/(4*sum(p((NX + 1)/2, :))/NY));
% box on;
% hold off;
