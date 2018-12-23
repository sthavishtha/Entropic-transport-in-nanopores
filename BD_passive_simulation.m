%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Brownian dynamics simulation in corrugated channels
% Course project : Transport Phenomena I
% Contributors   : Christoph, Sthavishtha, Tess
% Copyright
% All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;

% Simulation parameters
NTRA = 5e3; NTIME = 1000; NHIST = 100; DT = 0.001;
XMIN = -1.5; DX = 0.005; XMAX = 1.5;
YMIN = -2.1; DY = 0.005; YMAX = 2.1;
Xedges = XMIN:DX:XMAX;
Yedges = YMIN:DY:YMAX;
Xcenters = XMIN+DX/2:DX:XMAX-DX/2;
Ycenters = YMIN+DY/2:DY:YMAX-DY/2;

% Top & Bottom wall parameters
% wall offset
amp = 1;
% wall amplitude
off = 1.05;  

syms B_u(t) B_d(t)
% top wall variation
B_u(t) = amp*cos(2*t*pi) + off;
% B_u(t) = amp*sin(2*pi*(t - 0.25)) + off;      
% bottom wall variation
B_d(t) = -B_u(t);                                                           
                                                        
% Force along x-direction
fx_mag = [0];                                                  

lin_t = 0:DT:NTIME*DT; 
lin_x = Xedges;

% Generation of initial NTRA trajectories
x = zeros(NTRA, NTIME);
y = zeros(NTRA, NTIME);

% Effective diffusion coefficients
Deffx = zeros(NTIME + 1, NHIST);
Deffy = zeros(NTIME + 1, NHIST);
Var_x = zeros(NTIME + 1, NHIST);
Var_y = zeros(NTIME + 1, NHIST);

% for the absorbing boundary conditions into account, removing these particles away
AbsParts = false(NTRA,1); 
outside = false(NTRA, 1);

% Mobility coefficients
Mobx = zeros(NTIME + 1, NHIST);
Moby = zeros(NTIME + 1, NHIST);

for i = 1:size(fx_mag,2)
    fx = fx_mag(i)*ones(NTRA, 1);
    
    if fx_mag(i) ~= 0
          Dx = 0.1375;
          Dy = 0.6119;
        
    else
        Dx = 1;
        Dy = 1;
    end
    
% Repeat trajectories NHIST times
for k = 1:NHIST
k
    % Iterate across time until steady state (NTIME)
    tic
    for j = 1:NTIME
        
        % Difference equation
        x(:, j + 1) = x(:, j) + fx*DT + sqrt(Dx*DT)*randn(NTRA, 1);
        y(:, j + 1) = y(:, j) + sqrt(Dy*DT)*randn(NTRA, 1);
        
        % Reflective boundary conditions at top/bottom walls
        % if x(j+1) lies outside the domain --> reset to x(j) [Burada's thesis, 2008]        
        BU = double(amp.*cos(2*pi*x(:, j + 1)) + off*ones(NTRA, 1)) - abs(y(:, j + 1));
        
        % returns logical TRUE if particle is outside
        out = BU < 0;        
        % indices of brownian particles lying outside 
        out_ind = find(out);
        x(out_ind, j + 1) = x(out_ind, j);
        y(out_ind, j + 1) = y(out_ind, j);
                
        % Absorbing BCs --> letting the particles escape out of the wall
        ind_out = (x(:,j) > XMAX) | (x(:,j) < XMIN);
        AbsParts = or(AbsParts, ind_out); % finds those particles which have escaped    
         
        % Diffusion coefficient for nonadsorbed particles (to be used in absorbing bcs only)
        Deffx(j+1, k) = var(x(not(AbsParts), j + 1))/(j*DT);
        Deffy(j+1, k) = var(y(not(AbsParts), j + 1))/(j*DT);
        Var_x(j+1, k) = var(x(not(AbsParts), j + 1));
        Var_y(j+1, k) = var(y(not(AbsParts), j + 1));
        
    end
        
    px(k, :) = histc(x(not(AbsParts), j + 1), Xedges)/(DX*size(find(not(AbsParts)), 1));
    py(k, :) = histc(y(not(AbsParts), j + 1), Yedges)/(DY*size(find(not(AbsParts)), 1));
    
end
toc

    % avg diffusion coefficients and mobility currents
    Deffx_avg = mean(Deffx, 2);
    Deffy_avg = mean(Deffy, 2);
    Varx_avg = mean(Var_x, 2);
    Vary_avg = mean(Var_y, 2);
    
end

% Post-processing
% Plot of simulation results and initial condition
% FName=strcat('_',num2str(NTRA),'Particles_',num2str(NTIME),'x',num2str(DT),'Steps_Geom',num2str(amp),'x',num2str(off),'_f_',num2str(fx));
fname = strcat('BD_t_', num2str(NTIME*DT), '_ntra_', num2str(NTRA), '_ntime_', num2str(NTIME), '_dt_', num2str(DT), ....
    '_amp_', num2str(amp), '_off_', num2str(off), '_f_', num2str(fx_mag));
save(strcat(fname, '.mat'));

steps = NTIME;
Times=[ 1                      %Initial Distribution
        ceil(steps*0.01) + 1   %After 1% of the time  
        ceil(steps*0.1)        %After 10% of the time...
        ceil(steps*0.25)
        ceil(steps*0.50)
        floor(steps*1)];
    
set(groot,'defaultAxesTickLabelInterpreter','latex');      

% Plot of initial particle positions
figure(1)
Tname = title(strcat('NTRA = ', num2str(NTRA),', NTIME = ', num2str(NTIME), ...
        'NHIST = ', num2str(NHIST),', DT = ',num2str(DT), ...
        'Force = ', num2str(fx_mag)));
set(Tname, 'interpreter', 'latex');        
subplot(2,1,1)
Tname = title('t = 0 (Delta distribution)');
yname = ylabel('Y direction');
xname = xlabel('X direction');
set([Tname, xname, yname], 'interpreter', 'latex');
hold on
plot(lin_x, B_u(lin_x),'black');
plot(lin_x, B_d(lin_x),'black');
plot(x(:, 1),y(:, 1),'.','DisplayName',strcat('t = ',num2str(1*DT)));
box on
grid on
hold off

% Plot of steady state particle positions
subplot(2,1,2)
hold on
Tname = title(strcat('t = ', num2str(DT*NTIME),' (final state)'));
yname = ylabel('Y direction');
xname = xlabel('X direction');
set([Tname, xname, yname], 'interpreter', 'latex');
plot(lin_x, B_u(lin_x), 'black');
plot(lin_x, B_d(lin_x), 'black'); 
plot(x(not(AbsParts), NTIME + 1), y(not(AbsParts), NTIME + 1), '.r', 'DisplayName', strcat('t = ',num2str(NTIME*DT)));
grid on
box on
hold off

% Plot of Variance distributions
figure(3)
Tname = title('Temporal variation of variance');
yname = ylabel('Variance');
xname = xlabel('Time');
set([Tname, xname, yname], 'interpreter', 'latex');
hold on
plot(lin_t, Varx_avg, '-')
plot(lin_t, Vary_avg, '-')
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
lname = legend('Variance x', 'Variance y', 'Variance2 x', 'Variance2 y');
set(lname, 'interpreter', 'latex');
axis([0 NTIME*DT 0 inf])
grid on
box on
hold off

% Plot of steady state normalized prob distributions
figure(4)
hold on
Tname = title('Steady state probability distributions');
yname = ylabel('Probability');
xname = xlabel('X direction');
set([Tname, xname, yname], 'interpreter', 'latex');
plot([Xcenters NaN], mean(px)/(2*sum(mean(px))/size(px, 2)), 'k--', 'LineWidth',2);
errorbar([Xcenters NaN], mean(px)/(2*sum(mean(px))/size(px, 2)), std(px)/sqrt(NHIST)/(2*sum(mean(px))/size(px, 2)), 'LineStyle', 'none','Color', [  0    0.4470    0.7410]);
grid on
box on
hold off;

% xa = [Xcenters NaN];
% e = std(px)/sqrt(NHIST)/(2*sum(mean(px))/size(px, 2));
% lo = mean(px)/(2*sum(mean(px))/size(px, 2)) - e;
% hi = mean(px)/(2*sum(mean(px))/size(px, 2)) + e;
% hp = patch([xa; xa(end:-1:1); xa(1)], [lo; hi(end:-1:1); lo(1)], 'r');
% set(hp, 'facecolor', [1 0.8 0.8], 'edgecolor', 'none');

