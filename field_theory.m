% =========================================================================
% VECTOR FIELD THEORY IN MOTION: REVEALING LATENT POTENTIALS IN SOCCER
% =========================================================================
%
% DESCRIPTION:
%   This script implements the computational framework to analyze football 
%   dynamics using vector field theory. It performs the following steps:
%       1. Generates a vector field (Synthetic Dipole for demo).
%       2. Computes the Curl to analyze field properties.
%       3. Verifies Gauss's Theorem.
%       4. Reconstructs the Scalar Potential (Phi).
%
% USAGE:
%   - Run the script to visualize the theoretical framework using the 
%     included synthetic data.
%   - To use real data: Replace Section 2 with your own 'Wx' and 'Wy' 
%     matrices derived from tracking data.
%
% DEPENDENCIES:
%   - Requires 'brewermap' function for color schemes. If not available,
%     replace with standard Matlab colormaps.
%
% =========================================================================

clearvars; close all; clc;

%% 1. GRID INITIALIZATION
% Define the spatial domain representing the football pitch
M = 20;                 % Number of cells per axis (Grid resolution)
x_min = 0; x_max = 100; % Pitch length (x)
y_min = 0; y_max = 100; % Pitch width (y)

% Generate cell edges and grid (21 points -> 20 cells)
x_edges = linspace(x_min, x_max, M+1);          % Cell edges in X (21x1)
y_edges = linspace(y_min, y_max, M+1);          % Cell edges in Y (21x1)
[X_grid, Y_grid] = meshgrid(x_edges, y_edges);  % 21x21 grid

% Compute cell centers 
x_centers = (x_edges(1:end-1) + x_edges(2:end)) / 2; % X centers of cells
y_centers = (y_edges(1:end-1) + y_edges(2:end)) / 2; % Y centers of cells
[centerX, centerY] = meshgrid(x_centers, y_centers); % 20x20 grid

%% 2. DEFINE VECTOR FIELD (SYNTHETIC DEMO)
% =========================================================================
% NOTE: This section generates an ARTIFICIAL DIPOLE field.
% In a real application, 'Wx' and 'Wy' should be loaded here from 
% spatially binned average velocity matrices derived from player tracking.
% =========================================================================

% A. Dipole Setup (Source vs Sink)
offset = 10;
rs = [x_min + offset, (y_max - y_min) / 2]; % Source position (Left)
rt = [x_max - offset, (y_max - y_min) / 2]; % Sink position (Right)

% B. Calculate vector directions
Rx_s = centerX - rs(1);
Ry_s = centerY - rs(2);
Rs2  = Rx_s.^2 + Ry_s.^2 + 1e-6; % Avoid division by zero

Rx_t = centerX - rt(1);
Ry_t = centerY - rt(2);
Rt2  = Rx_t.^2 + Ry_t.^2 + 1e-6;

% Superposition of Source and Sink
Wx_base = (Rx_s ./ Rs2 - Rx_t ./ Rt2);
Wy_base = (Ry_s ./ Rs2 - Ry_t ./ Rt2);

% Normalize to get pure direction
norm_base = sqrt(Wx_base.^2 + Wy_base.^2);
Wx_dir = Wx_base ./ (norm_base + eps);
Wy_dir = Wy_base ./ (norm_base + eps);

% C. Define Velocity
% We use a displaced Gaussian bell curve to break symmetry for Gauss theorem
v_max_val = 11;
v_min_val = 2;

peak_x = 35;  % Velocity peaks at x=35 (not the center)
width  = 30;  % Width of the velocity bell curve


exponent = -((centerX - peak_x).^2) / (2 * width^2);
v_ball = v_min_val + (v_max_val - v_min_val) * exp(exponent);

% D. Apply Magnitude to Direction
Wx = Wx_dir .* v_ball;
Wy = Wy_dir .* v_ball;

% Expand arrays for visualization (padding for pcolor/surf edges)
v_exp = padarray(v_ball, [1 1], 'replicate', 'post');

%% 3. VISUALIZATION: VECTOR FIELD
figure;
% Heatmap of velocity magnitude
h = surf(X_grid, Y_grid, zeros(size(v_exp)), v_exp);
set(gca, 'YDir', 'normal');         
h.FaceColor = 'interp';             
h.EdgeColor = 'none';                   
colormap(brewermap(256, 'ylgnbu'));       

cb = colorbar;
cb.Label.FontSize = 16;
cb.Label.FontName = 'Helvetica';
cb.FontSize = 14;

axis equal;                                            
hold on;                                                

% Draw grid lines
for i = 1:M+1
    line([x_edges(i), x_edges(i)], [y_min, y_max], ...
         'Color', [0.5 0.5 0.5], 'LineStyle', '-');
    line([x_min, x_max], [y_edges(i), y_edges(i)], ...
         'Color', [0.5 0.5 0.5], 'LineStyle', '-');
end

% Overlay Quiver (Arrows)
quiver(centerX, centerY, Wx, Wy, 'k');

% Overlay football field
plot_field();


view(2);
xlim([x_min, x_max]);
ylim([y_min, y_max]);
axis off;     

%% 4. CURL COMPUTATION (ROTATIONAL ANALYSIS)
% Comparison of the field's curl against a randomized Null Model.

% A. Null Model Generation
% Randomize directions while preserving magnitude
mag = sqrt(Wx.^2 + Wy.^2);                           % Magnitude of each vector in the field
theta = 2*pi * rand(size(Wx));                       % Random angles uniformly distributed in [0, 2π]

% Reassign vector components: preserve magnitude but assign random directions
Wx_null = mag .* cos(theta);                         % x-component with randomized direction
Wy_null = mag .* sin(theta);                         % y-component with randomized direction

% B. Compute Curl
[curlz, cav] = curl(centerX, centerY, Wx, Wy);       % Curl of original field
mod_curl_W = abs(curlz);                             % Magnitude of curl
[curlz_null, cav_null] = curl(centerX, centerY, Wx_null, Wy_null); 
mod_curl_W_null = abs(curlz_null);

% Mean values for comparison
curlz_mean = mean(mod_curl_W(:));
curlz_null_mean = mean(mod_curl_W_null(:));

% C. Histogram / PDF Analysis
min_val = min([min(curlz(:)), min(curlz_null(:))]);
max_val = max([max(curlz(:)), max(curlz_null(:))]);
num_bins = 50;
edges = linspace(min_val, max_val, num_bins+1);

[counts_field_raw, ~] = histcounts(curlz,     edges, 'Normalization', 'pdf');
[counts_null_raw,  ~] = histcounts(curlz_null, edges, 'Normalization', 'pdf');

% Normalize distributions to range [0, 1] for comparison
max_val_global = max([counts_field_raw, counts_null_raw]);
counts_field   = counts_field_raw / max_val_global;
counts_null    = counts_null_raw  / max_val_global;

centers = (edges(1:end-1) + edges(2:end)) / 2;

% Plot PDF Comparison
figure;
hold on;
color_field = [0.2, 0.6, 0.8];   % Blue for field
color_null  = [0.8, 0.4, 0.2];   % Orange for null model
plot(centers, counts_field, '-o', 'DisplayName', 'Field', ...
    'Color', color_field, 'MarkerSize', 10, 'LineWidth', 1, ...
    'MarkerFaceColor', color_field, 'MarkerEdgeColor', 'k');
plot(centers, counts_null, '-o', 'DisplayName', 'Null model', ...
    'Color', color_null, 'MarkerSize', 10, 'LineWidth', 1, ...
    'MarkerFaceColor', color_null, 'MarkerEdgeColor', 'k');
xlabel('Curl', 'FontSize', 30);
ylabel('PDF', 'FontSize', 30, 'Rotation', 0);
title('Normalized curl distribution', 'FontSize', 26);
legend('Location', 'northeast', 'FontSize', 24, 'Box', 'off');
set(gca, 'FontSize', 22, 'FontName', 'Helvetica');
grid on;
box on;
xlim([min(centers), max(centers)]);
ylim([0, max([counts_field, counts_null])]);
hold off;

% Plot Spatial Curl Magnitude
figure();

% Original field curl magnitude
subplot(1,2,1)
pcolor(mod_curl_W);
hold on;
title('$\|\nabla \times \vec{\mathbf{W}}\|$', 'Interpreter', 'latex');
xlabel('x');
ylabel('y');
axis equal;
xlim([1 20]);
ylim([1 20]);

% Null model curl magnitude
subplot(1,2,2)
pcolor(mod_curl_W_null);
hold on;
title('$\|\nabla \times \vec{\mathbf{W_{null}}}\|$', 'Interpreter', 'latex');
xlabel('x');
ylabel('y');
axis equal;
xlim([1 20]);
ylim([1 20]);

min_value = min(min(mod_curl_W(:)), min(mod_curl_W_null(:)));
max_value = max(max(mod_curl_W(:)), max(mod_curl_W_null(:)));

% Shared colorbar for both subplots
c = colorbar('Position', [0.93, 0.11, 0.02, 0.815]);
c.Label.String = 'Curl magnitude';
caxis([min_value, max_value]);

subplot(1,2,1); caxis([min_value, max_value]);
subplot(1,2,2); caxis([min_value, max_value]);

%% 5. GAUSS'S THEOREM VERIFICATION

% A. Setup
h = (X_grid(1,M+1) - X_grid(1,1)) / M;                 % Grid spacing

x_a = 50; y_a = 50;                                    % Center coordinates
divW = divergence(centerX, centerY, Wx, Wy);           % Divergence of vector field
radios = linspace(0, 50, 20);                          % Radios to evaluate

integrales_div   = zeros(size(radios));                % Divergence integrals
flujos_contorno  = zeros(size(radios));                % Boundary fluxes

% B. Calculation loop (Circular Domain)
for r_idx = 1:length(radios)
    R = radios(r_idx);
    
    % 1. Volume Integral (Sum of Divergence inside)
    mask = ((centerX-x_a).^2 + (centerY-y_a).^2 <= R^2);
    integral_div = sum(divW(mask)) * h^2; 
    integrales_div(r_idx) = integral_div;

    % 2. Surface Integral (Flux across boundary)
    theta = linspace(0, 2*pi, 100); 
    x_circ = x_a + R * cos(theta);
    y_circ = y_a + R * sin(theta);

    % Interpolate field at boundary
    Wx_circ = interp2(centerX, centerY, Wx, x_circ, y_circ, 'linear'); 
    Wy_circ = interp2(centerX, centerY, Wy, x_circ, y_circ, 'linear'); 

    % Outward unit normals (radial direction)
    normales_x = x_circ - x_a; 
    normales_y = y_circ - y_a;  

    norm_base = sqrt(normales_x.^2 + normales_y.^2); 

    normales_x = normales_x ./ norm_base; 
    normales_y = normales_y ./ norm_base; 

    % Line integral
    flux_contorno = 0;
    for i = 1:length(x_circ)-1
        dx = x_circ(i+1) - x_circ(i);
        dy = y_circ(i+1) - y_circ(i);
        longitud_segmento = sqrt(dx^2 + dy^2);
        
        flujo_segmento = (Wx_circ(i) * normales_x(i) + Wy_circ(i) * normales_y(i)) * longitud_segmento;
        flux_contorno = flux_contorno + flujo_segmento;
    end

    flujos_contorno(r_idx) = flux_contorno;
end


% C. Visualization of Gauss Theorem
figure();
hold on;

color1 = [0.2, 0.6, 0.8];  % Blue
color2 = [0.8, 0.4, 0.2];  % Orange

% Surface
plot(radios, flujos_contorno, '-s', 'Color', color1, 'MarkerSize', 18, 'LineWidth', 1, ...
    'MarkerFaceColor', color1, 'MarkerEdgeColor', 'k', 'DisplayName', 'Surface');

% Volume
plot(radios, integrales_div, '-^', 'Color', color2, 'MarkerSize', 18, 'LineWidth', 1, ...
    'MarkerFaceColor', color2, 'MarkerEdgeColor', 'k', 'DisplayName', 'Volume');

% Compute coefficient R^2
valid_indices = ~isnan(flujos_contorno) & ~isnan(integrales_div);
flujos_validos = flujos_contorno(valid_indices);
integrales_validos = integrales_div(valid_indices);

p = polyfit(integrales_validos, flujos_validos, 1);
flujos_p = polyval(p, integrales_validos);

SS_res = sum((flujos_validos - flujos_p).^2);
SS_tot = sum((flujos_validos - mean(flujos_validos)).^2);
R2 = 1 - (SS_res / SS_tot);

corr_coef = corrcoef(flujos_validos, integrales_validos);

text(0.65 * max(radios), 0.15 * max(flujos_contorno), ...
     sprintf('$R^2 = %.2f$', R2), 'FontSize', 26, ...
     'Interpreter', 'latex', 'Color', [0.2 0.2 0.2]);

set(gca, 'FontSize', 22, 'FontName', 'Helvetica');
xlabel('$R$', 'Interpreter', 'latex', 'FontSize', 30);
ylabel('$\phi$', 'Interpreter', 'latex', 'FontSize', 30, 'Rotation',0);
legend('Location', 'best', 'FontSize', 24, 'Box', 'off');

ax = gca;
minVal = min(radios);
maxVal = max(radios);
numTicks = 6; % Número de ticks que deseas

ax.XTick = linspace(minVal, maxVal, numTicks);
ax.XTickLabel = arrayfun(@(x) sprintf('%.0f', x), ax.XTick, 'UniformOutput', false);

box on;
hold off;


% Square domains
for r_idx = 1:length(radios)
    R = radios(r_idx);
    
    mask = (centerX >= (x_a - R) & centerX <= (x_a + R) & ...
            centerY >= (y_a - R) & centerY <= (y_a + R));
    
    integral_div = sum(divW(mask)) * h^2; 
    integrales_div(r_idx) = integral_div;

    x_cuadrado = [x_a - R, x_a + R, x_a + R, x_a - R, x_a - R];
    y_cuadrado = [y_a - R, y_a - R, y_a + R, y_a + R, y_a - R];

    N = 25; 
    flux_contorno = 0;            
    
    for side = 1:4
        if side == 1  % Left side
            x_points = linspace(x_a - R, x_a - R, N);
            y_points = linspace(y_a - R, y_a + R, N);
            normales_x = -1; normales_y = 0;  
        elseif side == 2  % Top side
            x_points = linspace(x_a - R, x_a + R, N);
            y_points = linspace(y_a + R, y_a + R, N);
            normales_x = 0; normales_y = 1;  
        elseif side == 3  % Right side
            x_points = linspace(x_a + R, x_a + R, N);
            y_points = linspace(y_a + R, y_a - R, N);
            normales_x = 1; normales_y = 0;
        else  % Bottom side
            x_points = linspace(x_a + R, x_a - R, N);
            y_points = linspace(y_a - R, y_a - R, N);
            normales_x = 0; normales_y = -1;
        end
        % Interpolate vector field on side
        Wx_segment = interp2(centerX, centerY, Wx, x_points, y_points, 'linear'); 
        Wy_segment = interp2(centerX, centerY, Wy, x_points, y_points, 'linear'); 

        % Flux per segment (segment length = 2R/N)
        for i = 1:N
            L = 2 * R / N; 
            flujo_segmento = (Wx_segment(i) * normales_x + Wy_segment(i) * normales_y) * L; 
            flux_contorno = flux_contorno + flujo_segmento;
        end
    end

    flujos_contorno(r_idx) = flux_contorno;

end



% Visualization (square domains)
figure();
hold on;

% Surface
plot(radios, flujos_contorno, '-s', 'Color', color1, 'MarkerSize', 18, 'LineWidth', 1, ...
    'MarkerFaceColor', color1, 'MarkerEdgeColor', 'k', 'DisplayName', 'Surface');

% Volume
plot(radios, integrales_div, '-^', 'Color', color2, 'MarkerSize', 18, 'LineWidth', 1, ...
    'MarkerFaceColor', color2, 'MarkerEdgeColor', 'k', 'DisplayName', 'Volume');

% R^2 calculation
valid_indices = ~isnan(flujos_contorno) & ~isnan(integrales_div);
flujos_validos = flujos_contorno(valid_indices);
integrales_validos = integrales_div(valid_indices);

p = polyfit(integrales_validos, flujos_validos, 1);
flujos_p = polyval(p, integrales_validos);

SS_res = sum((flujos_validos - flujos_p).^2);
SS_tot = sum((flujos_validos - mean(flujos_validos)).^2);
R2 = 1 - (SS_res / SS_tot);
corr_coef = corrcoef(flujos_validos, integrales_validos);

text(0.65 * max(radios), 0.15 * max(flujos_contorno), ...
     sprintf('$R^2 = %.2f$', R2), 'FontSize', 26, ...
     'Interpreter', 'latex', 'Color', [0.2 0.2 0.2]);

set(gca, 'FontSize', 22, 'FontName', 'Helvetica');
xlabel('$R$', 'Interpreter', 'latex', 'FontSize', 30);
ylabel('$\Phi$', 'Interpreter', 'latex', 'FontSize', 30, 'Rotation',0);
legend('Location', 'best', 'FontSize', 24, 'Box', 'off');

ax = gca;
minVal = min(radios);
maxVal = max(radios);
numTicks = 6;

ax.XTick = linspace(minVal, maxVal, numTicks);
ax.XTickLabel = arrayfun(@(x) sprintf('%.0f', x), ax.XTick, 'UniformOutput', false);

box on;
hold off;

%% 6. POTENTIAL CALCULATION

% A. Pad Data for Solver
% We add padding to handle boundary conditions 
x_min_ext = x_min - h;
x_max_ext = x_max + h;
y_min_ext = y_min - h;
y_max_ext = y_max + h;
[X, Y] = meshgrid(linspace(x_min_ext, x_max_ext, M+3), linspace(y_min_ext, y_max_ext, M+3));

centerX = (X(1:end-1, 1:end-1) + X(2:end, 2:end)) / 2;
centerY = (Y(1:end-1, 1:end-1) + Y(2:end, 2:end)) / 2;  

% Extend arrays with zero padding
Z = padarray(v_ball, [1,1], 0, 'both');
Wx = padarray(Wx, [1,1], 0, 'both');
Wy = padarray(Wy, [1,1], 0, 'both');


% B. Compute Source Term (Divergence)
nx = size(Wx, 2);
ny = size(Wx, 1);
h = (X(1,M+3)-X(1,1))/(M+2);
Phi = zeros(ny, nx);         
divV = zeros(ny, nx);

% Central differences for divergence
divV(:, 2:end-1) = (Wx(:, 3:end) - Wx(:, 1:end-2)) / (2*h);
divV(2:end-1, :) = divV(2:end-1, :) + (Wy(3:end, :) - Wy(1:end-2, :)) / (2*h);

% Poisson Equation Source: RHS = -Div(V) * h^2
RHS = -divV * h^2;

% C. Jacobi Iteration
max_iter = 5000;
tol = 1e-6;

for k = 1:max_iter
    Phi_old = Phi;
    
    % Vectorized Update: Average of 4 neighbors
    Phi(2:end-1, 2:end-1) = 0.25 * (Phi_old(1:end-2, 2:end-1) + ... % Up/Down
                                    Phi_old(3:end,   2:end-1) + ... 
                                    Phi_old(2:end-1, 1:end-2) + ... % Left
                                    Phi_old(2:end-1, 3:end)   - ... % Right
                                    RHS(2:end-1, 2:end-1));
    
    % Check Convergence
    diff = abs(Phi - Phi_old);
    if max(diff(:)) < tol
        break;
    end
end


%% 7. VISUALIZATION: 3D SCALAR POTENTIAL
fig  = figure;

surf(centerX, centerY, Phi, Phi);     % phi is also used as colormap
shading interp;     
hold on;
[Xp, Yp] = meshgrid(0:100, 0:100);
Zp = zeros(size(Xp)); 
surf(Xp, Yp, Zp, 'FaceAlpha', 0.2, 'FaceColor', 'k', 'EdgeColor', 'none');  
contour3(centerX, centerY, Phi, 30, 'k', 'LineWidth', 1);  

colormap(brewermap(256, 'RdBU'));  

col_idx = round(nx * 0.2);



function [] = plot_field()
% PLOT_FIELD Plots a football pitch scaled to 0-100 units.

    % 1. AXIS SETUP AND ASPECT RATIO
    daspect([0.68 1.05 1]);
    axis([0 100 0 100]);
    
    % 2. SCALING FACTORS
    xlimit1=0;
    ylimit1=0;
    xlimit2=100;
    ylimit2=100;
    xscale=xlimit2/100/1.05;
    yscale=ylimit2/100/0.68;
    
    hold on;

    % 3. LINE STYLES
    color_linea = [.7 .7 .7];
    estil_linea = "-";
    ancho_linea = 0.1;
    
    % 4. CENTER CIRCLE
    x=50;
    y=50;
    rx=9.15/1.05;
    ry=9.15/0.68;
    
    th = 0:pi/50:2*pi;
    xunit = rx * cos(th) + x;
    yunit = ry * sin(th) + y;
    plot(xunit, yunit, 'Color', color_linea, 'LineWidth', ancho_linea, LineStyle=estil_linea);
    
    % 5. PENALTY ARCS
    x=11/1.05;
    y=50;
    rx=9.15/1.05;
    ry=9.15/0.68;
    th = -0.93:1/100:0.93;
    xunit = rx * cos(th) + x;
    yunit = ry * sin(th) + y;
    plot(xunit, yunit, 'Color', color_linea, 'LineWidth', ancho_linea, LineStyle=estil_linea);
    
    x=100-11/1.05;
    y=50;
    rx=9.15/1.05;
    ry=9.15/0.68;
    th = pi-0.93:1/100:pi+0.93;
    xunit = rx * cos(th) + x;
    yunit = ry * sin(th) + y;
    plot(xunit, yunit, 'Color', color_linea, 'LineWidth', ancho_linea, LineStyle=estil_linea);
    
    
    % 6. PENALTY AREAS
    x=[xlimit1 xlimit1+(16.5*xscale) xlimit1+(16.5*xscale) xlimit1];
    y=[ylimit2/2+(20*yscale) ylimit2/2+(20*yscale) ylimit2/2-(20*yscale) ylimit2/2-(20*yscale)];
    plot(x, y, 'Color', color_linea, 'LineWidth', ancho_linea, LineStyle=estil_linea);
    
    x=[xlimit1 xlimit1+(5.5*xscale) xlimit1+(5.5*xscale) xlimit1];
    y=[ylimit2/2+(9*yscale) ylimit2/2+(9*yscale) ylimit2/2-(9*yscale) ylimit2/2-(9*yscale)];
    plot(x, y, 'Color', color_linea, 'LineWidth', ancho_linea, LineStyle=estil_linea);
    
    x=[xlimit2 xlimit2-(16.5*xscale) xlimit2-(16.5*xscale) xlimit2];
    y=[ylimit2/2+(20*yscale) ylimit2/2+(20*yscale) ylimit2/2-(20*yscale) ylimit2/2-(20*yscale)];
    plot(x,y, 'Color', color_linea, 'LineWidth', ancho_linea, LineStyle=estil_linea);
    
    x=[xlimit2 xlimit2-(5.5*xscale) xlimit2-(5.5*xscale) xlimit2];
    y=[ylimit2/2+(9*yscale) ylimit2/2+(9*yscale) ylimit2/2-(9*yscale) ylimit2/2-(9*yscale)];
    plot(x, y, 'Color', color_linea, 'LineWidth', ancho_linea, LineStyle=estil_linea);
    
    % 7. OUTER BOUNDARIES & CENTER LINE
    color_linea = 'k';%[0.5 0.5 0.5];

    x=[xlimit1 xlimit2 xlimit2 xlimit1 xlimit1];
    y=[ylimit1 ylimit1 ylimit2 ylimit2 ylimit1];
    plot(x, y, 'Color', color_linea, 'LineWidth', ancho_linea, LineStyle=estil_linea);
    
    x=[xlimit2/2 xlimit2/2];
    y=[ylimit1 ylimit2];
    plot(x, y, 'Color', color_linea, 'LineWidth', 0.1, LineStyle=estil_linea);
    
    hold off;


end

