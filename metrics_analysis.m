%% === MONTE CARLO COMPARISON WITH D* AND BURNT AREA ANALYSIS ===
function comparison_results = monte_carlo_comparison(obstacleFile, num_trials, p, num_fires)
% Run Monte Carlo comparison between RRT, A*, and D*
% Now includes burnt area analysis

if nargin < 4, num_fires = 1; end
if nargin < 3, p = 0.45; end
if nargin < 2, num_trials = 1; end
if nargin < 1, obstacleFile = 'obstacleMap.mat'; end

% Initialise storage
results_rrt = [];
results_astar = [];
results_dstar = [];

fprintf('Progress: [');
for trial = 1:num_trials
    if mod(trial-1, max(1,floor(num_trials/50))) == 0
        fprintf('█');
    end
    
    seed = trial * 100;
    
    % --- RRT Trial ---
    try
        rng(seed);
        metrics_rrt = osm_fire_sim(obstacleFile, p, num_fires, true, 'rrt');
        close all;
        results_rrt = [results_rrt; metrics_rrt];
    catch ME
        warning('RRT trial %d failed: %s', trial, ME.message);
        nan_row = table(NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
    'VariableNames', {'Time_s','Path_m','AvgTurn_rad','MaxTurn_rad','ThermalDanger_s','BurntArea_m2','OES'});

        results_rrt = [results_rrt; nan_row];
    end
    
    % --- A* Trial (same seed) ---
    try
        rng(seed);
        metrics_astar = osm_fire_sim(obstacleFile, p, num_fires, true, 'astar');
        close all;
        results_astar = [results_astar; metrics_astar];
    catch ME
        warning('A* trial %d failed: %s', trial, ME.message);
        nan_row = table(NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
    'VariableNames', {'Time_s','Path_m','AvgTurn_rad','MaxTurn_rad','ThermalDanger_s','BurntArea_m2','OES'});

        results_astar = [results_astar; nan_row];
    end
    
    % --- D* Trial (same seed) ---
    try
        rng(seed);
        metrics_dstar = osm_fire_sim(obstacleFile, p, num_fires, true, 'dstar');
        close all;
        results_dstar = [results_dstar; metrics_dstar];
    catch ME
        warning('D* trial %d failed: %s', trial, ME.message);
        nan_row = table(NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
    'VariableNames', {'Time_s','Path_m','AvgTurn_rad','MaxTurn_rad','ThermalDanger_s','BurntArea_m2','OES'});

        results_dstar = [results_dstar; nan_row];
    end
end
fprintf('] Complete!\n\n');

% Calculate statistics
comparison_results = struct();
comparison_results.num_trials = num_trials;
comparison_results.rrt = results_rrt;
comparison_results.astar = results_astar;
comparison_results.dstar = results_dstar;

% Compute metrics
metrics = results_rrt.Properties.VariableNames;
for i = 1:length(metrics)
    metric = metrics{i};
    comparison_results.stats.rrt.(metric).mean = mean(results_rrt.(metric), 'omitnan');
    comparison_results.stats.rrt.(metric).std = std(results_rrt.(metric), 'omitnan');
    comparison_results.stats.astar.(metric).mean = mean(results_astar.(metric), 'omitnan');
    comparison_results.stats.astar.(metric).std = std(results_astar.(metric), 'omitnan');
    comparison_results.stats.dstar.(metric).mean = mean(results_dstar.(metric), 'omitnan');
    comparison_results.stats.dstar.(metric).std = std(results_dstar.(metric), 'omitnan');
    
    % ANOVA test 
    valid_idx = ~isnan(results_rrt.(metric)) & ~isnan(results_astar.(metric)) & ~isnan(results_dstar.(metric));
    if sum(valid_idx) > 3
        groups = [ones(sum(valid_idx),1); 2*ones(sum(valid_idx),1); 3*ones(sum(valid_idx),1)];
        values = [results_rrt.(metric)(valid_idx); results_astar.(metric)(valid_idx); results_dstar.(metric)(valid_idx)];
        [p_anova, ~, stats_anova] = anova1(values, groups, 'off');
        comparison_results.stats.(metric).anova_p = p_anova;
        
        % Post-hoc pairwise t-tests
        [h_rrt_astar, p_rrt_astar] = ttest(results_rrt.(metric)(valid_idx), results_astar.(metric)(valid_idx));
        [h_rrt_dstar, p_rrt_dstar] = ttest(results_rrt.(metric)(valid_idx), results_dstar.(metric)(valid_idx));
        [h_astar_dstar, p_astar_dstar] = ttest(results_astar.(metric)(valid_idx), results_dstar.(metric)(valid_idx));
        
        comparison_results.stats.(metric).pairwise.rrt_astar_p = p_rrt_astar;
        comparison_results.stats.(metric).pairwise.rrt_dstar_p = p_rrt_dstar;
        comparison_results.stats.(metric).pairwise.astar_dstar_p = p_astar_dstar;
    else
        comparison_results.stats.(metric).anova_p = NaN;
    end
end

% Display results
display_monte_carlo_results(comparison_results);

% Plot
plot_monte_carlo_comparison(comparison_results);

% Save
timestamp = datestr(now, 'yyyymmdd_HHMMSS');


end

%% Display Function
function display_monte_carlo_results(results)
    fprintf('\n╔═══════════════════════════════════════════════════════════╗\n');
    fprintf('║                    RESULTS SUMMARY                        ║\n');
    fprintf('╚═══════════════════════════════════════════════════════════╝\n\n');
    
    metrics = {'Time_s', 'Path_m', 'AvgTurn_rad', 'MaxTurn_rad', ...
           'ThermalDanger_s', 'BurntArea_m2', 'OES'};

    names = {'Computation Time (s)', 'Path Length (m)', 'Avg Turn (rad)', ...
         'Max Turn (rad)', 'Thermal Danger (s)', 'Burnt Area (m²)', 'Opportunity Exposure Score'};
    
    for i = 1:length(metrics)
        metric = metrics{i};
        fprintf('─── %s ───\n', names{i});
        fprintf('  RRT:  %.3f ± %.3f\n', ...
            results.stats.rrt.(metric).mean, results.stats.rrt.(metric).std);
        fprintf('  A*:   %.3f ± %.3f\n', ...
            results.stats.astar.(metric).mean, results.stats.astar.(metric).std);
        fprintf('  D*:   %.3f ± %.3f\n', ...
            results.stats.dstar.(metric).mean, results.stats.dstar.(metric).std);
        
        p_val = results.stats.(metric).anova_p;
        if ~isnan(p_val)
            if p_val < 0.001
                fprintf('  *** HIGHLY SIGNIFICANT DIFFERENCE (ANOVA p < 0.001)\n');
            elseif p_val < 0.05
                fprintf('  ** SIGNIFICANT DIFFERENCE (ANOVA p = %.4f)\n', p_val);
            else
                fprintf('  No significant difference (ANOVA p = %.4f)\n', p_val);
            end
            
            % Pairwise comparisons
            if isfield(results.stats.(metric), 'pairwise')
                fprintf('  Pairwise comparisons:\n');
                if results.stats.(metric).pairwise.rrt_astar_p < 0.05
                    fprintf('    • RRT vs A*: p = %.4f **\n', results.stats.(metric).pairwise.rrt_astar_p);
                end
                if results.stats.(metric).pairwise.rrt_dstar_p < 0.05
                    fprintf('    • RRT vs D*: p = %.4f **\n', results.stats.(metric).pairwise.rrt_dstar_p);
                end
                if results.stats.(metric).pairwise.astar_dstar_p < 0.05
                    fprintf('    • A* vs D*:  p = %.4f **\n', results.stats.(metric).pairwise.astar_dstar_p);
                end
            end
        end
        
        % Which algorithm is best
        planners = {'RRT', 'A*', 'D*'};
        vals = [results.stats.rrt.(metric).mean, ...
                results.stats.astar.(metric).mean, ...
                results.stats.dstar.(metric).mean];
        
        if strcmp(metric, 'OES')
            % For OES, HIGHER is better
            [best_val, best_idx] = max(vals);
        else
            % For all other metrics, LOWER is better
            [best_val, best_idx] = min(vals);
        end
        
        fprintf('  → BEST: %s (%.3f)\n', planners{best_idx}, best_val);

        fprintf('\n');
    end
end

%% Plotting Function
function plot_monte_carlo_comparison(results)

    figure('Position', [100 100 1600 900], ...
           'Name', 'Monte Carlo Comparison - RRT vs A* vs D*');

    % Create clean layout
    tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Metrics
    metrics = {'Time_s', 'Path_m', 'AvgTurn_rad', 'MaxTurn_rad', ...
               'ThermalDanger_s', 'OES'};
    titles = {'Computation Time (s)', 'Path Length (m)', 'Avg Turn (rad)', ...
              'Max Turn (rad)', 'Thermal Danger (s)', ...
              'Opportunity Exposure Score'};

    num_metrics = length(metrics);

    for i = 1:num_metrics
        nexttile;

        metric = metrics{i};

        %get data
        rrt_data   = results.rrt.(metric);
        astar_data = results.astar.(metric);
        dstar_data = results.dstar.(metric);

        % Filter NaNs
        rrt_valid   = rrt_data(~isnan(rrt_data));
        astar_valid = astar_data(~isnan(astar_data));
        dstar_valid = dstar_data(~isnan(dstar_data));

        % Combine for boxplot
        data = [rrt_valid; astar_valid; dstar_valid];
        group = [ones(length(rrt_valid),1);
                 2*ones(length(astar_valid),1);
                 3*ones(length(dstar_valid),1)];

        % Plot boxplot
        boxplot(data, group, ...
                'Labels', {'RRT', 'A*', 'D*'}, ...
                'Colors', 'rbg');

        ylabel(titles{i});
        grid on;
        hold on;

        % Add mean markers
        plot(1, mean(rrt_valid),   'r*', 'MarkerSize', 15, 'LineWidth', 2);
        plot(2, mean(astar_valid), 'b*', 'MarkerSize', 15, 'LineWidth', 2);
        plot(3, mean(dstar_valid), 'g*', 'MarkerSize', 15, 'LineWidth', 2);

        % Add p-value to title
        p_val = results.stats.(metric).anova_p;
        if ~isnan(p_val)
            if p_val < 0.001
                title(sprintf('%s\n(ANOVA p < 0.001 ***)', titles{i}));
            elseif p_val < 0.05
                title(sprintf('%s\n(ANOVA p = %.4f **)', titles{i}, p_val));
            else
                title(sprintf('%s\n(ANOVA p = %.4f)', titles{i}, p_val));
            end
        else
            title(titles{i});
        end
    end

    % Main title 
    sgtitle(sprintf('Monte Carlo Comparison (%d Trials)', results.num_trials), ...
            'FontSize', 16, 'FontWeight', 'bold');

end

function metrics = osm_fire_sim(obstacleFile, p, num_fires, use_uav, planner_type)
% OSM Fire Simulation with Thermal Camera and UAV Path Planning
% obstacleFile: path to MAT file with obstacles (e.g., 'obstacleMap.mat')
% p: fire spread probability (default 0.4)
% num_fires: number of initial fire points (default 1)
% use_uav: enable UAV path planning (default true)
% planner_type: 'rrt' or 'astar' (default 'rrt')

if nargin < 5, planner_type = 'astar'; end
if nargin < 4, use_uav = true; end
if nargin < 3, num_fires = 1; end
if nargin < 2, p = 0.4; end
if nargin < 1, obstacleFile = 'obstacleMap.mat'; end

sim_start_time = tic;

% Validate planner type
planner_type = lower(planner_type);
if ~ismember(planner_type, {'rrt','astar','dstar'})
    error('planner_type must be ''rrt'', ''astar'', or ''dstar''');
end
disp(['Using path planner: ', upper(planner_type)]);

%% === 1. Load Pre‑computed Obstacle Map ===
disp(['Loading obstacle map from: ', obstacleFile]);
try
    loaded = load(obstacleFile);
    obstacles = loaded.obstacles;
    refLat    = loaded.refLat;
    refLon    = loaded.refLon;
catch ME
    error(['Failed to load obstacle map: ', ME.message]);
end
disp(['Loaded ', num2str(numel(obstacles)), ' building obstacles']);

xLimits = [-500 500];   yLimits = [-500 500];

%% === 2. Create Grid Representation ===
gridRes = 10;
xMin = xLimits(1); xMax = xLimits(2);
yMin = yLimits(1); yMax = yLimits(2);
N_x = ceil((xMax-xMin)/gridRes);
N_y = ceil((yMax-yMin)/gridRes);

C = ones(N_x,N_y);                 % 0=empty,1=veg,2=building,3=fire,4=burned
T = ones(N_x,N_y)*20;              % ambient temperature
H = zeros(N_x,N_y);                % building height
buildingIndex = zeros(N_x,N_y);
uav_coverage = false(N_x,N_y);
buildingStates = zeros(numel(obstacles),1);

%% === 3. Rasterize Buildings ===
disp('Rasterizing buildings...');
for i = 1:numel(obstacles)
    verts = obstacles(i).Vertices;
    if ~isempty(verts)
        x_coords = verts(:,1); y_coords = verts(:,2); z_coords = verts(:,3);
        h = max(z_coords);
        ix1 = max(1, floor((min(x_coords)-xMin)/gridRes)+1);
        ix2 = min(N_x, ceil((max(x_coords)-xMin)/gridRes)+1);
        iy1 = max(1, floor((min(y_coords)-yMin)/gridRes)+1);
        iy2 = min(N_y, ceil((max(y_coords)-yMin)/gridRes)+1);
        for ix = ix1:ix2
            for iy = iy1:iy2
                wx = xMin + (ix-0.5)*gridRes;
                wy = yMin + (iy-0.5)*gridRes;
                if inpolygon(wx,wy,x_coords,y_coords)
                    C(ix,iy) = 2; H(ix,iy) = h; buildingIndex(ix,iy) = i;
                end
            end
        end
    end
end
disp(['Grid size: ', num2str(N_x), 'x', num2str(N_y)]);
disp(['Buildings: ', num2str(sum(C(:)==2))]);

%% === 4. Initialise Fires ===
disp('Initializing fires...');
targetLat = 63.4215217; targetLon = 10.3997097;
mPerDegLat = 111320;
mPerDegLon = 111320 * cosd(refLat);
dx = (targetLon-refLon)*mPerDegLon;
dy = (targetLat-refLat)*mPerDegLat;
fx = round((dx-xMin)/gridRes)+1;
fy = round((dy-yMin)/gridRes)+1;

fires = []; fire_locations = [];

% ---- forced fire (if inside grid) ----
if fx>=1 && fx<=N_x && fy>=1 && fy<=N_y
    if C(fx,fy)~=3 && C(fx,fy)~=4
        C(fx,fy) = 3;
        T(fx,fy) = 600 + 300*(H(fx,fy)>0);
        if H(fx,fy)>0, buildingStates(buildingIndex(fx,fy))=1; end
        fires = [fires; fx fy];
        fire_z = max(H(fx,fy),5);
        fire_locations = [fire_locations; dx dy fire_z];
        disp(['Placed forced fire at lat/lon ',num2str(targetLat),', ',num2str(targetLon)]);
    else
        warning('Forced cell already burning – using random placement.');
    end
else
    warning('Forced lat/lon outside grid – using random placement.');
end

% ---- fill remaining fires ----
if size(fires,1) < num_fires
    remaining = num_fires - size(fires,1);
    attempts = 0;
    while remaining>0 && attempts<1000
        attempts = attempts+1;
        rx = randi(N_x); ry = randi(N_y);
        if C(rx,ry)~=3 && C(rx,ry)~=4
            C(rx,ry) = 3;
            T(rx,ry) = 600 + 300*(H(rx,ry)>0);
            if H(rx,ry)>0, buildingStates(buildingIndex(rx,ry))=1; end
            fires = [fires; rx ry];
            wx = xMin + (rx-0.5)*gridRes;
            wy = yMin + (ry-0.5)*gridRes;
            fire_z = max(H(rx,ry),5);
            fire_locations = [fire_locations; wx wy fire_z];
            remaining = remaining-1;
        end
    end
end
[numfires,~] = size(fires);
disp(['Initial fires: ', num2str(numfires)]);

%% === 5. UAV Path Planning Setup ===
uav_path = []; uav_position_idx = 0;
uav_trajectory = [];  %  Store complete UAV trajectory
replanning_count = 0;  %  Track replanning events
fires_extinguished = 0;  %  Track extinguished fires
fire_count_history = [];  %  Track fire count over time
coverage_history = [];  %  Track coverage percentage
thermal_exposure_history = [];  %  Track UAV thermal exposure
path_length_history = [];  %  Track path lengths
in_fire_zone = [];  %  Track if UAV is in fire zone

if use_uav && numfires>0
    uav_start = [260 -16 0];
    gx = round((uav_start(1)-xMin)/gridRes)+1;
    gy = round((uav_start(2)-yMin)/gridRes)+1;
    if gx<1 || gx>N_x || gy<1 || gy>N_y
        gx = max(1,min(gx,N_x)); gy = max(1,min(gy,N_y));
        uav_start(1) = xMin + (gx-0.5)*gridRes;
        uav_start(2) = yMin + (gy-0.5)*gridRes;
    end
    if H(gx,gy)>0 || T(gx,gy)>200
        uav_start(3) = max(uav_start(3), max(H(gx,gy)+20,30));
    end

    disp('Planning initial UAV path...');
    [uav_path,ok] = planToAllFires(uav_start, fire_locations, planner_type, ...
                                   T, H, C, xMin, yMin, gridRes, N_x, N_y);
    
    if ok
        disp(['Initial path: ',num2str(size(uav_path,1)),' waypoints']);
        uav_position_idx = 1;
        uav_trajectory = [uav_trajectory; uav_start];  % Start trajectory
    else
        disp('Path planning failed – UAV disabled');
        use_uav = false;
    end
end

%% === 6. Figures ===
fig1 = figure(1); set(fig1,'Position',[100 100 900 700],...
    'Name',['OSM Fire Sim + UAV (',upper(planner_type),')']);
fig2 = figure(2); set(fig2,'Position',[1050 100 900 700],...
    'Name','Thermal Camera + UAV Coverage');
fig3 = figure(3); set(fig3,'Position',[100 850 900 300],...
    'Name','Real-Time Metrics Dashboard');

%% === 7. Simulation Loop ===
maxsteps = 500; k = 0;

% ---- real‑time fire‑avoidance tracking ----
close_passes_during_flight = 0;
total_waypoints_checked    = 0;

while k<maxsteps && numfires>0
    lastC = C; numfires_next = 0; fires_next = [];

    % ---- temperature update ----
    T_new = T;
    for i=1:N_x, for j=1:N_y
        if C(i,j)==3
            T_new(i,j) = 550 + rand*100 + 300*(H(i,j) > 0);
        elseif C(i,j)==4
            T_new(i,j) = max(20,T(i,j)*0.85);
        else
            T_new(i,j) = 20;
        end
    end, end
    for i=2:N_x-1, for j=2:N_y-1
        if C(i,j)==3
            T_new(i-1,j) = max(T_new(i-1,j),T_new(i,j)*0.3);
            T_new(i+1,j) = max(T_new(i+1,j),T_new(i,j)*0.3);
            T_new(i,j-1) = max(T_new(i,j-1),T_new(i,j)*0.3);
            T_new(i,j+1) = max(T_new(i,j+1),T_new(i,j)*0.3);
        end
    end, end
    T = T_new;
    
    % ---- fire spread ----
    for i=1:numfires
        fx = fires(i,1); fy = fires(i,2);
        nbr = [fx-1 fy; fx+1 fy; fx fy-1; fx fy+1];
        for n=1:4
            nx = nbr(n,1); ny = nbr(n,2);
            if nx<1||nx>N_x||ny<1||ny>N_y, continue; end
            if uav_coverage(nx,ny)
                prob = p * 0.2;
            else
                prob = p;
            end
            if (C(nx,ny)==1||C(nx,ny)==2) && rand<prob
                C(nx,ny)=3; numfires_next=numfires_next+1;
                fires_next(numfires_next,:)=[nx ny];
                if buildingIndex(nx,ny)>0, buildingStates(buildingIndex(nx,ny))=1; end
            end
        end
        C(fx,fy)=4;
        if buildingIndex(fx,fy)>0, buildingStates(buildingIndex(fx,fy))=2; end
    end
    fires = fires_next; numfires = numfires_next;

    % ---- update fire_locations ----
    fire_locations = zeros(numfires,3);
    for i=1:numfires
        fx = fires(i,1); fy = fires(i,2);
        fire_locations(i,:) = [xMin+(fx-0.5)*gridRes, yMin+(fy-0.5)*gridRes, max(H(fx,fy),5)];
    end

    % ---- UAV move + real‑time avoidance check ----
    current_thermal_exposure = 0;
    if use_uav && ~isempty(uav_path) && uav_position_idx>0
        uav_position_idx = min(uav_position_idx+1, size(uav_path,1));
        uav_pos = uav_path(uav_position_idx,:);
        uav_trajectory = [uav_trajectory; uav_pos];  %  Add to trajectory
        
        ux = round((uav_pos(1)-xMin)/gridRes)+1;
        uy = round((uav_pos(2)-yMin)/gridRes)+1;

        % Get thermal exposure
        if ux>=1 && ux<=N_x && uy>=1 && uy<=N_y
            current_thermal_exposure = T(ux,uy);
        end

        % ---- is UAV too close to an active fire? (3×3 cell radius) ----
        too_close = false;
        for dx=-1:1, for dy=-1:1
            nx = ux+dx; ny = uy+dy;
            if nx>=1&&nx<=N_x&&ny>=1&&ny<=N_y && C(nx,ny)==3
                too_close = true;
            end
        end, end
        if too_close, close_passes_during_flight = close_passes_during_flight+1; end
        total_waypoints_checked = total_waypoints_checked+1;

        % ---- suppression area ----
        supR = 1;
        for dx=-supR:supR, for dy=-supR:supR
            sx = ux+dx; sy = uy+dy;
            if sx>=1&&sx<=N_x&&sy>=1&&sy<=N_y
                uav_coverage(sx,sy)=true;
                if C(sx,sy)==3
                    T(sx,sy)=max(50,T(sx,sy)*0.3);
                    if rand<0.4 && T(sx,sy)<100
                        C(sx,sy)=1; T(sx,sy)=30;
                        fires_extinguished = fires_extinguished + 1;  
                        disp(['UAV extinguished fire at [',num2str(sx),',',num2str(sy),']']);
                    end
                end
            end
        end, end
    end
    is_in_fire_zone = false;
    if use_uav && uav_position_idx > 0
        search_radius = 3; % cells
        for dx = -search_radius:search_radius
            for dy = -search_radius:search_radius
                nx = ux + dx;
                ny = uy + dy;
                if nx >= 1 && nx <= N_x && ny >= 1 && ny <= N_y
                    if C(nx, ny) == 3  % Active fire
                        is_in_fire_zone = true;
                        break;
                    end
                end
            end
            if is_in_fire_zone
                break;
            end
        end
    end
    in_fire_zone(end+1) = is_in_fire_zone;
    % ---- replan every 3 steps ----
    if use_uav && ~isempty(uav_path) && mod(k,3)==0 && uav_position_idx>0 && numfires>0
        cur = uav_path(uav_position_idx,:);
        [newp,ok] = planToAllFires(cur, fire_locations, planner_type, T,H,C,xMin,yMin,gridRes,N_x,N_y);
        if ok && ~isempty(newp)
            uav_path = newp; uav_position_idx = 1;
            replanning_count = replanning_count + 1;  
            disp(['Replanned at step ',num2str(k), ' (Total replans: ', num2str(replanning_count), ')']);
        end
    end

    % NEW: Store metrics history
    fire_count_history(end+1) = sum(C(:)==3);
    coverage_history(end+1) = 100*sum(uav_coverage(:))/(N_x*N_y);
    thermal_exposure_history(end+1) = current_thermal_exposure;
    % Compute turn angles for entire trajectory so far
    if size(uav_trajectory,1) >= 3
        diffs = diff(uav_trajectory(:,1:3),1,1);
        dirs = diffs ./ vecnorm(diffs,2,2);
        turn_angles = acos(max(-1,min(1,dot(dirs(1:end-1,:), dirs(2:end,:), 2))));
        % Prepend 0 for initial steps 
        turn_angle_time_series = [0; turn_angles];
    else
        turn_angle_time_series = NaN(size(uav_trajectory,1),1);
    end


    if use_uav && ~isempty(uav_path)
        path_length_history(end+1) = size(uav_path,1);
    else
        path_length_history(end+1) = 0;
    end

    %% === 8. 3‑D visualisation ===
    figure(fig1); clf; hold on;
    [Xg,Yg] = meshgrid(xMin:50:xMax, yMin:50:yMax);
    surf(Xg,Yg,zeros(size(Xg)),'FaceColor',[0.3 0.6 0.3],'EdgeColor','none','FaceAlpha',0.4);

    % UAV coverage
    for i=1:N_x, for j=1:N_y, if uav_coverage(i,j)
        wx = xMin+(i-0.5)*gridRes; wy = yMin+(j-0.5)*gridRes;
        patch([wx-gridRes/2 wx+gridRes/2 wx+gridRes/2 wx-gridRes/2],...
              [wy-gridRes/2 wy-gridRes/2 wy+gridRes/2 wy+gridRes/2],...
              [0.1 0.1 0.1 0.1],[0.2 0.6 1.0],'EdgeColor','none','FaceAlpha',0.3);
    end, end, end

    % Buildings
    for i=1:numel(obstacles)
        if ~isempty(obstacles(i).Vertices)
            if buildingStates(i)==0
                patch('Vertices',obstacles(i).Vertices,'Faces',obstacles(i).Faces,...
                      'FaceColor',[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.8);
            elseif buildingStates(i)==1
                patch('Vertices',obstacles(i).Vertices,'Faces',obstacles(i).Faces,...
                      'FaceColor',[0.3 0.15 0.05],'EdgeColor','none','FaceAlpha',0.9);
            else
                v=obstacles(i).Vertices; v(:,3)=v(:,3)*0.7;
                patch('Vertices',v,'Faces',obstacles(i).Faces,...
                      'FaceColor',[0.1 0.1 0.1],'EdgeColor','none','FaceAlpha',0.7);
            end
        end
    end

    % Fires / burned ground
    for i=1:N_x, for j=1:N_y
        wx = xMin+(i-0.5)*gridRes; wy = yMin+(j-0.5)*gridRes;
        if C(i,j)==3 && H(i,j)>0, drawFire(wx,wy,H(i,j),gridRes*1.5);
        elseif C(i,j)==3, drawVegetationFire(wx,wy,gridRes);
        elseif C(i,j)==4 && H(i,j)==0, drawBurnedGround(wx,wy,gridRes);
        end
    end, end

   % === UAV Trajectory with Colour Gradient + Colourbar ===
    if use_uav && size(uav_trajectory,1) > 1
        num_points = size(uav_trajectory,1);
    
        % Build colormap 
        cmap = zeros(num_points, 3);
        for i = 1:num_points
            t = (i-1) / (num_points-1); % 0 → 1 along path
            
            if t < 0.5
                % Red → Orange → Yellow
                red   = 0.6 + 0.4 * (t/0.5);
                green = 0.2 + 0.5 * (t/0.5);
                blue  = 0.0;
            else
                % Yellow → Light Green
                red   = 1.0 - 0.8 * ((t-0.5)/0.5);
                green = 0.7 + 0.3 * ((t-0.5)/0.5);
                blue  = 0.0 + 0.2 * ((t-0.5)/0.5);
            end
            
            cmap(i,:) = [red, green, blue];
        end
    
        % Plot trajectory segments using the same colours
        for i = 1:num_points-1
            plot3(uav_trajectory(i:i+1,1), ...
                  uav_trajectory(i:i+1,2), ...
                  uav_trajectory(i:i+1,3), ...
                  'Color', cmap(i,:), 'LineWidth', 2.5);
            hold on;
        end
    
        % --- Add colourbar showing Start -> End ---
        colormap(cmap);
    
                % Invisible dummy surface to activate colour scaling
        cdata = linspace(0,1,num_points);
        surf([0 0; 0 0], [0 0; 0 0], [0 0; 0 0], [0 1; 0 1], ...
             'EdgeColor', 'none', 'FaceColor', 'none');
        
        c = colorbar;
        
        % Label for the colourbar
        c.Label.String = 'Path Progress';
        
        % Remove default tick labels so we can add manual text
        c.TickLabels = [];
        
        % Get colourbar position
        cbPos = c.Position;  % [x y width height]
        
        % Add manual 'Start' label at the bottom right, just right of the colourbar
        text(cbPos(1) + cbPos(3) + 0.25, cbPos(2)-0.3, 'Start', ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
             'Units', 'normalized');
        
        % Add manual 'End' label at the top right, just right of the colourbar
        text(cbPos(1) + cbPos(3) + 0.25, cbPos(2) + cbPos(4)+0.3, 'End', ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
             'Units', 'normalized');


    end


    % UAV path (planned path in BLUE DASHED)
    if use_uav && ~isempty(uav_path)
        plot3(uav_path(:,1),uav_path(:,2),uav_path(:,3),'b--','LineWidth',2);
        if uav_position_idx>0
            drawUAV(uav_path(uav_position_idx,1),uav_path(uav_position_idx,2),uav_path(uav_position_idx,3));
        end
        plot3(uav_path(1,1),uav_path(1,2),uav_path(1,3),'go','MarkerSize',10,'MarkerFaceColor','g');
        plot3(uav_path(end,1),uav_path(end,2),uav_path(end,3),'ro','MarkerSize',10,'MarkerFaceColor','r');
    end

    view(3); axis equal; camlight; lighting gouraud; grid on;
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    if use_uav
        if total_waypoints_checked > 0
            avoidPct = 100 * (1 - close_passes_during_flight / total_waypoints_checked);
        else
            avoidPct = 0;
        end
        title(sprintf('Step %d | Fires:%d | Burned:%d | UAV cov:%d%% | Avoid:%.1f%% | Replans:%d',...
            k,sum(C(:)==3),sum(C(:)==4),round(100*sum(uav_coverage(:))/(N_x*N_y)),avoidPct,replanning_count));
    else
        title(sprintf('Step %d | Fires:%d | Burned:%d',k,sum(C(:)==3),sum(C(:)==4)));
    end
    xlim([xMin xMax]); ylim([yMin yMax]); zlim([0 max(max(H(:))+50,80)]);

    %% === 9. Thermal view ===
    figure(fig2); clf;
    comb = T';
    for i=1:N_x, for j=1:N_y, if uav_coverage(i,j), comb(j,i)=-50; end, end, end
    imagesc(xMin:gridRes:xMax, yMin:gridRes:yMax, comb);
    colormap([0.2 0.6 1.0; hot(64)]); caxis([-100 1100]); hold on;
    if use_uav && uav_position_idx>0
        plot(uav_pos(1),uav_pos(2),'wo','MarkerSize',15,'LineWidth',3);
        plot(uav_pos(1),uav_pos(2),'go','MarkerSize',10,'MarkerFaceColor','g');
        if size(uav_trajectory,1) > 1
            plot(uav_trajectory(:,1),uav_trajectory(:,2),'g-','LineWidth',2);
        end
        plot(uav_path(:,1),uav_path(:,2),'w--','LineWidth',1);
    end
    patch([xMin+50 xMin+70 xMin+70 xMin+50],[yMin+50 yMin+50 yMin+70 yMin+70],[0.2 0.6 1.0],...
          'EdgeColor','none','FaceAlpha',0.5);
    text(xMin+75,yMin+60,'UAV Coverage','Color','white','FontWeight','bold');
    colorbar; axis equal tight; set(gca,'YDir','normal');
    xlabel('X [m]'); ylabel('Y [m]');
    title(sprintf('Thermal + Coverage – Step %d | Fires:%d | Suppressed:%d | Extinguished:%d',...
        k,sum(C(:)==3),sum(uav_coverage(:)),fires_extinguished));

    %% === 10.  Metrics Dashboard ===
  
    figure(fig3); clf;
    
    % Plot 1: Thermal Damage (Exposure)
    subplot(1,2,1);
    plot(thermal_exposure_history,'m-','LineWidth',2);
    grid on; xlabel('Step'); ylabel('Temperature (°C)');
    title('UAV Thermal Exposure');
    hold on;
    plot([0 length(thermal_exposure_history)], [150 150], 'r--', 'LineWidth', 1);
    legend('Exposure', 'Danger Threshold', 'Location', 'best');
    
    % Plot 2: Turn Angles
    subplot(1,2,2);
    plot(rad2deg(turn_angle_time_series), 'b-', 'LineWidth', 2);
    grid on;
    xlabel('Step');
    ylabel('Turn Angle (deg)');
    title('Turn Angle Over Time');
    
    
    drawnow; pause(0.1);
    k = k+1;
    if isequal(C,lastC), disp('Fire stopped spreading.'); break; end
end

%% === Final Metrics (Thermal + Smoothness + OES) ===
%% Ensure sim_start_time exists
if ~exist('sim_start_time','var')
    sim_time_sec = tic;   % start timing now
end

sim_time_sec = toc(sim_start_time);

% Calculate burnt area (cells with C==4, multiplied by cell area)
burnt_cells = sum(C(:)==4);
burnt_area_m2 = burnt_cells * (gridRes^2);


if use_uav && exist('uav_trajectory','var') && ~isempty(uav_trajectory)
    % Path Length
    diffs = diff(uav_trajectory(:,1:3),1,1);
    total_path_length = sum(vecnorm(diffs,2,2));

    % Smoothness
    if size(diffs,1) >= 2
        dirs = diffs ./ vecnorm(diffs,2,2);
        turn_angles = acos(max(-1, min(1, dot(dirs(1:end-1,:), dirs(2:end,:), 2))));
        avg_smoothness = mean(abs(turn_angles));
        max_turn_angle = max(abs(turn_angles));
    else
        avg_smoothness = NaN;
        max_turn_angle = NaN;
    end

    % === Thermal Danger Time ===
    if exist('thermal_exposure_history','var') && ~isempty(thermal_exposure_history)
        danger_threshold = 150; % °C
        time_per_step = sim_time_sec / length(thermal_exposure_history);

        dangerous_steps = sum(thermal_exposure_history > danger_threshold);
        thermal_danger_time = dangerous_steps * time_per_step;

        fprintf('Thermal danger analysis: %d steps above %d°C = %.2f seconds\n', ...
                dangerous_steps, danger_threshold, thermal_danger_time);
    else
        thermal_danger_time = 0;
        warning('Thermal exposure history not available for danger calculation');
    end

    % === Opportunity Exposure Score (OES) ===
    % Arrival index = first time step inside the fire zone
    % note you must have a boolean vector "in_fire_zone" same length as thermal_exposure_history
    if exist('in_fire_zone','var') && ~isempty(in_fire_zone)
        arrival_index = find(in_fire_zone == true, 1, 'first');
    else
        arrival_index = [];
    end

    if isempty(arrival_index)
        OES = 0;   % never entered fire zone
    else
        % Split thermal exposure
        potential_exposure = sum(thermal_exposure_history(1:arrival_index));      % before arrival
        actual_exposure     = sum(thermal_exposure_history(arrival_index:end));   % after arrival

        % Avoid division by zero
        OES = potential_exposure / (potential_exposure + actual_exposure + 1e-6);
    end

else
    total_path_length = NaN;
    avg_smoothness = NaN;
    max_turn_angle = NaN;
    thermal_danger_time = NaN;
    OES = NaN;
end

%% === Add to Metrics Table ===
metrics = table(sim_time_sec, total_path_length, avg_smoothness, max_turn_angle, ...
    thermal_danger_time, burnt_area_m2, OES, ...
    'VariableNames', {'Time_s','Path_m','AvgTurn_rad','MaxTurn_rad','ThermalDanger_s','BurntArea_m2','OES'});

disp('====== UAV Simulation Metrics ======');
disp(metrics);

%% === Save results ===
results_struct = struct(...
    'timestamp', datetime('now'), ...
    'planner_type', planner_type, ...
    'metrics', metrics, ...
    'trajectory', uav_trajectory, ...
    'thermal_history', thermal_exposure_history, ...
    'OES', OES);


end



% -------------------------------------------------------------------------

%% === Multi-Goal Planning Helper ===
function [uav_path, success] = planToAllFires(current_pos, fire_locations, planner_type, T, H, C, xMin, yMin, gridRes, N_x, N_y)    
    if isempty(fire_locations)
        uav_path = [];
        success = false;
        return;
    end
    
    % Sort fires by distance from current position (in XY plane)
    dists = vecnorm(fire_locations(:,1:2) - current_pos(1:2), 2, 2);
    [~, sort_idx] = sort(dists);
    sorted_fires = fire_locations(sort_idx, :);
    
    % Adjust target heights
    sorted_fires(:,3) = sorted_fires(:,3) + 30;
    
    % Chain paths to each goal sequentially
    path_segments = {};
    next_start = current_pos;
    success = true;
    
    for g = 1:size(sorted_fires,1)
        goal = sorted_fires(g,:);
        if strcmp(planner_type, 'rrt')
            [seg_path, planning_success] = planPathRRT3D(next_start, goal, T, H, C, xMin, yMin, gridRes, N_x, N_y);
        elseif strcmp(planner_type, 'dstar')
            [seg_path, planning_success] = planPathDStar3D(next_start, goal, T, H, C, xMin, yMin, gridRes, N_x, N_y);
        else
            [seg_path, planning_success] = planPathAStar3D(next_start, goal, T, H, C, xMin, yMin, gridRes, N_x, N_y);
        end
        
        if ~planning_success
            success = false;
            break;
        end
        
        path_segments{end+1} = seg_path;
        next_start = seg_path(end,:);
    end
    
    if ~success
        uav_path = [];
        return;
    end
    
    % Concatenate segments
    uav_path = path_segments{1};
    for s = 2:length(path_segments)
        uav_path = [uav_path; path_segments{s}(2:end,:)]; % Skip duplicate start point
    end
end

%% ------------------ D* 3D Planner  A*-fallback ------------------
function [path, success] = planPathDStar3D(start, goal, T, H, C, xMin, yMin, gridRes, N_x, N_y)
fprintf('\n=== D* 3D – Robust variant with A* fallback ===\n');

% Use same occupancy logic / z-resolution as A*
temp_threshold = 150;
zMax = 120;
zRes = 10;
N_z = ceil(zMax / zRes);

% Build occupancy grid (same logic as A*/RRT to keep behaviours consistent)
occGrid = false(N_x, N_y, N_z);
for i = 1:N_x
    for j = 1:N_y
        if H(i,j) > 0
            occGrid(i,j,1:N_z) = true;
        end
        if T(i,j) > temp_threshold
            z_idx_max = min(N_z, ceil(80 / zRes));
            occGrid(i, j, 1:z_idx_max) = true;
        end
        if C(i,j) == 3
            occGrid(i,j,1:min(N_z,ceil(40/zRes))) = true;
        end
        if C(i,j) == 4
            occGrid(i,j,1:min(N_z,ceil(15/zRes))) = true;
        end
    end
end

% Convert start/goal to grid indices consistent with A*/RRT
startIdx = worldToGrid([start(1), start(2), start(3)], xMin, yMin, 0, gridRes, zRes);
goalIdx  = worldToGrid([goal(1),  goal(2),  goal(3)],  xMin, yMin, 0, gridRes, zRes);
startIdx = max([1 1 1], min([N_x N_y N_z], round(startIdx)));
goalIdx  = max([1 1 1], min([N_x N_y N_z], round(goalIdx)));

% If start or goal inside obstacle, try small local adjustments (same approach as RRT/A*)
if occGrid(startIdx(1), startIdx(2), startIdx(3))
    found_safe = false;
    for radius = 1:3
        for dx = -radius:radius
            for dy = -radius:radius
                if dx == 0 && dy == 0, continue; end
                new_x = startIdx(1) + dx;
                new_y = startIdx(2) + dy;
                if new_x >= 1 && new_x <= N_x && new_y >= 1 && new_y <= N_y
                    for test_z = [startIdx(3), startIdx(3)+1:min(N_z, startIdx(3)+3)]
                        if ~occGrid(new_x, new_y, test_z)
                            startIdx = [new_x, new_y, test_z];
                            fprintf('D*: Adjusted start to [%d,%d,%d]\n', new_x, new_y, test_z);
                            found_safe = true; break;
                        end
                    end
                end
                if found_safe, break; end
            end
            if found_safe, break; end
        end
        if found_safe, break; end
    end
    if ~found_safe
        warning('D*: start contained in obstacles and no nearby free cell found — falling back to A*');
        [path, success] = planPathAStar3D(start, goal, T, H, C, xMin, yMin, gridRes, N_x, N_y);
        return;
    end
end

if occGrid(goalIdx(1), goalIdx(2), goalIdx(3))
    % If goal blocked, find nearest free
    newGoal = findNearestFreeCell3D(occGrid, goalIdx, 10);
    if isempty(newGoal)
        warning('D*: goal blocked and no nearby free cell — falling back to A*');
        [path, success] = planPathAStar3D(start, goal, T, H, C, xMin, yMin, gridRes, N_x, N_y);
        return;
    else
        fprintf('D*: Fire goal blocked, using nearest free cell [%d, %d, %d]\n', newGoal);
        goalIdx = newGoal;
    end
end

% --- Initialize rhs/gScore arrays  ---
rhs = inf(N_x, N_y, N_z);
gScore = inf(N_x, N_y, N_z);
rhs(goalIdx(1), goalIdx(2), goalIdx(3)) = 0;

% openSet columns: i j k key
openSet = [double(goalIdx(1)), double(goalIdx(2)), double(goalIdx(3)), 0];
success = false;

% heuristic for keys (Euclidean in grid coords)
h = @(a,b) norm(double(a)-double(b));

% Helper: update or insert into openSet with given key
    function upsertOpen(n, keyVal)
        if isempty(openSet)
            openSet = double([n, keyVal]);
            return;
        end
        % find existing entry
        idxs = find(all(openSet(:,1:3) == n, 2));
        if ~isempty(idxs)
            % keep smallest key
            if keyVal < openSet(idxs(1),4)
                openSet(idxs(1),4) = keyVal;
            end
        else
            openSet = [openSet; double([n, keyVal])];
        end
    end

% --- Propagation loop (limited, robust) ---
maxIter = 200000;
iter = 0;
while ~isempty(openSet) && iter < maxIter
    iter = iter + 1;
    % Pop node with smallest key
    [~, idmin] = min(openSet(:,4));
    u = round(openSet(idmin,1:3));
    openSet(idmin,:) = [];
    
    % Stop condition: when we have propagated to start 
    if isfinite(rhs(startIdx(1), startIdx(2), startIdx(3)))
        success = true;
        break;
    end
    
    % For each neighbor of u, attempt to relax rhs
    for di = -1:1
        for dj = -1:1
            for dk = -1:1
                if di==0 && dj==0 && dk==0, continue; end
                n = u + [di dj dk];
                if any(n < 1) || any(n > [N_x N_y N_z]), continue; end
                if occGrid(n(1), n(2), n(3)), continue; end
                cost = norm([di dj dk]);
                newrhs = cost + rhs(u(1), u(2), u(3));
                if newrhs < rhs(n(1), n(2), n(3))
                    rhs(n(1), n(2), n(3)) = newrhs;
                    keyVal = rhs(n(1), n(2), n(3)) + h(n, startIdx);
                    upsertOpen(n, keyVal);
                end
            end
        end
    end
end

% If propagation failed to reach start, fallback to A*
if ~success
    warning('D*: propagation failed (iter=%d) — falling back to A*', iter);
    [path, success] = planPathAStar3D(start, goal, T, H, C, xMin, yMin, gridRes, N_x, N_y);
    return;
end

% --- Path reconstruction from startIdx following rhs gradient to goalIdx ---
pathGrid = startIdx;
current = startIdx;
maxSteps = 20000;
steps = 0;
while ~isequal(current, goalIdx) && steps < maxSteps
    steps = steps + 1;
    bestNext = [];
    bestVal = inf;
    for di = -1:1
        for dj = -1:1
            for dk = -1:1
                if di==0 && dj==0 && dk==0, continue; end
                n = current + [di dj dk];
                if any(n < 1) || any(n > [N_x N_y N_z]), continue; end
                if occGrid(n(1), n(2), n(3)), continue; end
                v = rhs(n(1), n(2), n(3));
                if v < bestVal
                    bestVal = v;
                    bestNext = n;
                end
            end
        end
    end
    if isempty(bestNext) || ~isfinite(bestVal)
        % failed to find next step — fallback to A*
        warning('D*: reconstruction failed at step %d — falling back to A*', steps);
        [path, success] = planPathAStar3D(start, goal, T, H, C, xMin, yMin, gridRes, N_x, N_y);
        return;
    end
    pathGrid = [pathGrid; bestNext];
    current = bestNext;
end

% Convert to world coordinates and return
path = gridToWorld(pathGrid, xMin, yMin, 0, gridRes, zRes);

% If final waypoint is not close enough to true goal, append the exact goal as A*/RRT do
if ~isempty(path) && norm(path(end,:) - goal) < 30
    path = [path; goal];
end

success = true;
end

%% === RRT Path Planning Function ===
function [path, success] = planPathRRT3D(start, goal, T, H, C, xMin, yMin, gridRes, N_x, N_y)
    % RRT path planning using KD-tree nearest neighbor search

    temp_threshold = 150;
    zMax = 120;
    zRes = 10;
    N_z = ceil(zMax / zRes);




    % --- Build occupancy grid ---
    occGrid = false(N_x, N_y, N_z);
    for i = 1:N_x
        for j = 1:N_y
            if H(i,j) > 0
                occGrid(i,j,1:N_z) = true;
            end
            if T(i,j) > temp_threshold
                z_idx_max = min(N_z, ceil(80 / zRes));
                occGrid(i, j, 1:z_idx_max) = true;
            end
            if C(i,j) == 3
                occGrid(i,j,1:min(N_z,ceil(40/zRes))) = true;
            end
            if C(i,j) == 4
                occGrid(i,j,1:min(N_z,ceil(15/zRes))) = true;
            end
        end
    end

    % --- Convert start/goal to grid indices ---
    startIdx = worldToGrid(start, xMin, yMin, 0, gridRes, zRes);
    goalIdx  = worldToGrid(goal,  xMin, yMin, 0, gridRes, zRes);
    startIdx = max([1 1 1], min([N_x N_y N_z], startIdx));
    goalIdx  = max([1 1 1], min([N_x N_y N_z], goalIdx));

    if occGrid(startIdx(1), startIdx(2), startIdx(3))
            found_safe = false;
            for radius = 1:3  % Search in increasing radii
                for dx = -radius:radius
                    for dy = -radius:radius
                        if dx == 0 && dy == 0, continue; end
                        
                        new_x = startIdx(1) + dx;
                        new_y = startIdx(2) + dy;
                        
                        if new_x >= 1 && new_x <= N_x && new_y >= 1 && new_y <= N_y
                            % Try current height and slightly above
                            for test_z = [startIdx(3), startIdx(3)+1:min(N_z, startIdx(3)+3)]
                                if ~occGrid(new_x, new_y, test_z)
                                    startIdx = [new_x, new_y, test_z];
                                    fprintf('Adjusted start to [%d,%d,%d]\n', new_x, new_y, test_z);
                                    found_safe = true;
                                    break;
                                end
                            end
                            if found_safe, break; end
                        end
                    end
                    if found_safe, break; end
                end
                if found_safe, break; end
            end
        end

    % --- Parameters ---
    maxIter = 5000;
    stepSize = 5;
    goalBias = 0.2;

    nodes = startIdx;
    parents = zeros(maxIter,1);
    success = false;

    closestNode = startIdx;
    closestDistance = norm(startIdx - goalIdx);
    closestNodeIdx = 1;

    % --- Initialise KD-tree (rebuild periodically) ---
    Mdl = KDTreeSearcher(nodes);

    for iter = 1:maxIter
        % Random sample (goal bias)
        if rand < goalBias
            randPt = goalIdx;
        else
            randPt = [randi(N_x), randi(N_y), randi([ceil(N_z*0.6), N_z])];
        end

        % Rebuild KD-tree every 25 iterations for performance
        if mod(iter,25)==0
            Mdl = KDTreeSearcher(nodes);
        end

        % --- Nearest neighbor via KD-tree ---
        nearestIdx = knnsearch(Mdl, randPt);
        nearest = nodes(nearestIdx,:);

        dir = randPt - nearest;
        dnorm = norm(dir);
        if dnorm == 0, continue; end
        newNode = round(nearest + min(stepSize, dnorm)*dir/dnorm);
        newNode = max([1 1 1], min([N_x N_y N_z], newNode));

        % --- Collision check ---
        if ~checkCollisionSimple(nearest, newNode, occGrid, N_x, N_y, N_z)
            nodes = [nodes; newNode];
            parents(size(nodes,1)) = nearestIdx;

            currentDistance = norm(newNode - goalIdx);
            if currentDistance < closestDistance
                closestDistance = currentDistance;
                closestNode = newNode;
                closestNodeIdx = size(nodes,1);
            end

            if norm(newNode - goalIdx) <= stepSize
                success = true;
                break;
            end
        end
    end

    % --- Reconstruct path ---
    if ~success && closestDistance < 15
        success = true;
    end
    if ~success
        path = []; return;
    end

    pathIdx = closestNodeIdx;
    pathGrid = [];
    while pathIdx > 0
        pathGrid = [nodes(pathIdx,:); pathGrid];
        pathIdx = parents(pathIdx);
    end

    path = gridToWorld(pathGrid, xMin, yMin, 0, gridRes, zRes);
    if closestDistance <= 10
        path = [path; goal];
    end
end

%% === Simple Collision Check ===
function collision = checkCollisionSimple(pt1, pt2, occGrid, N_x, N_y, N_z)
    % Simple linear interpolation collision check
    steps = max(abs(pt2 - pt1));
    if steps == 0
        collision = occGrid(pt1(1), pt1(2), pt1(3));
        return;
    end
    
    for s = 0:steps
        t = s/steps;
        checkPt = round(pt1 * (1-t) + pt2 * t);
        checkPt = max([1 1 1], min([N_x N_y N_z], checkPt));
        if occGrid(checkPt(1), checkPt(2), checkPt(3))
            collision = true;
            return;
        end
    end
    collision = false;
end

%% === A* Path Planning Function ===
function [path, success] = planPathAStar3D(start, goal, T, H, C, xMin, yMin, gridRes, N_x, N_y)
    % Parameters
    temp_threshold = 150;
    zMax = 120;
    zRes = 10;
    N_z = ceil(zMax / zRes);
    
    % Build 3D occupancy grid
    occGrid = false(N_x, N_y, N_z);
    for i = 1:N_x
        for j = 1:N_y
            if H(i,j) > 0
                occGrid(i,j,1:N_z) = true;
            end
            if T(i,j) > temp_threshold
                z_idx_max = min(N_z, ceil(80 / zRes));
                occGrid(i, j, 1:z_idx_max) = true;
            end
            if C(i,j) == 3
                occGrid(i,j,1:min(N_z,ceil(40/zRes))) = true;
            end
            if C(i,j) == 4
                occGrid(i,j,1:min(N_z,ceil(15/zRes))) = true;
            end
        end
    end
    
    % Convert start and goal positions to grid indices
    startIdx = worldToGrid([start(1), start(2), start(3)], xMin, yMin, 0, gridRes, zRes);
    startIdx = max([1 1 1], min([N_x N_y N_z], startIdx));
    
    goalIdx = worldToGrid([goal(1), goal(2), goal(3)], xMin, yMin, 0, gridRes, zRes);
    goalIdx = max([1 1 1], min([N_x N_y N_z], goalIdx));
    
    % If goal cell is blocked, find nearest free neighbor
    if occGrid(goalIdx(1), goalIdx(2), goalIdx(3))
        newGoal = findNearestFreeCell3D(occGrid, goalIdx, 10);
        if isempty(newGoal)
            warning('No reachable free cell near fire goal found');
            path = [];
            success = false;
            return;
        else
            fprintf('Fire goal blocked, using nearest free cell [%d, %d, %d]\n', newGoal);
            goalIdx = newGoal;
        end
    end
    
    % Initialise A* data structures
    openSet = [double(startIdx(1)), double(startIdx(2)), double(startIdx(3)), 0, heuristic(startIdx, goalIdx)];
    closedSet = false(N_x, N_y, N_z);
    cameFrom = zeros(N_x, N_y, N_z, 3);
    gScore = inf(N_x, N_y, N_z);
    fScore = inf(N_x, N_y, N_z);
    
    gScore(startIdx(1), startIdx(2), startIdx(3)) = 0;
    fScore(startIdx(1), startIdx(2), startIdx(3)) = heuristic(startIdx, goalIdx);
    
    success = false;
    maxIter = 50000;
    iter = 0;
    
    while ~isempty(openSet) && iter < maxIter
        iter = iter + 1;
        [~, idx] = min(openSet(:,5));
        current = round(openSet(idx, 1:3));
        openSet(idx,:) = [];
        
        if isequal(current, goalIdx)
            success = true;
            break;
        end
        
        closedSet(current(1), current(2), current(3)) = true;
        
        for di = -1:1
            for dj = -1:1
                for dk = -1:1
                    if di == 0 && dj == 0 && dk == 0
                        continue;
                    end
                    
                    ni = current(1) + di;
                    nj = current(2) + dj;
                    nk = current(3) + dk;
                    
                    if ni < 1 || ni > N_x || nj < 1 || nj > N_y || nk < 1 || nk > N_z
                        continue;
                    end
                    
                    if occGrid(ni, nj, nk) || closedSet(ni, nj, nk)
                        continue;
                    end
                    
                    dist = norm([di, dj, dk]);
                    tentative_gScore = gScore(current(1), current(2), current(3)) + dist;
                    
                    if tentative_gScore < gScore(ni, nj, nk)
                        cameFrom(ni, nj, nk, :) = current;
                        gScore(ni, nj, nk) = tentative_gScore;
                        fScore(ni, nj, nk) = tentative_gScore + heuristic([ni, nj, nk], goalIdx);
                        
                        if ~any(all(openSet(:,1:3) == [ni, nj, nk], 2))
                            openSet = [openSet; double(ni), double(nj), double(nk), tentative_gScore, fScore(ni, nj, nk)];
                        end
                    end
                end
            end
        end
    end
    
    if ~success
        path = [];
        return;
    end
    
    % Reconstruct path
    pathGrid = goalIdx;
    current = goalIdx;
    maxPathSteps = 10000;
    pathSteps = 0;
    
    while ~isequal(current, startIdx) && pathSteps < maxPathSteps
        pathSteps = pathSteps + 1;
        parent_vec = squeeze(cameFrom(current(1), current(2), current(3), :));
        parent = round(parent_vec)';
        if all(parent == 0)
            warning('A* path reconstruction failed - broken parent chain');
            path = [];
            success = false;
            return;
        end
        pathGrid = [parent; pathGrid];
        current = parent;
    end
    
    path = gridToWorld(pathGrid, xMin, yMin, 0, gridRes, zRes);
end


function nearestFree = findNearestFreeCell3D(occGrid, goalIdx, maxRadius)
    if nargin < 3
        maxRadius = 10; % default max search radius
    end
    N_x = size(occGrid,1);
    N_y = size(occGrid,2);
    N_z = size(occGrid,3);
    for r = 1:maxRadius
        for dx = -r:r
            for dy = -r:r
                for dz = -r:r
                    nx = goalIdx(1) + dx;
                    ny = goalIdx(2) + dy;
                    nz = goalIdx(3) + dz;
                    if nx < 1 || nx > N_x || ny < 1 || ny > N_y || nz < 1 || nz > N_z
                        continue;
                    end
                    if ~occGrid(nx, ny, nz)
                        nearestFree = [nx, ny, nz];
                        return;
                    end
                end
            end
        end
    end
    nearestFree = [];
end



%% === Helper Functions ===

function h = heuristic(pt, goal)
    h = norm(pt - goal);
end

function gridIdx = worldToGrid(worldPt, xMin, yMin, zMin, gridRes, zRes)
    gridIdx = round([(worldPt(1)-xMin)/gridRes, (worldPt(2)-yMin)/gridRes, (worldPt(3)-zMin)/zRes]) + 1;
end

function worldPt = gridToWorld(gridIdx, xMin, yMin, zMin, gridRes, zRes)
    worldPt = [(gridIdx(:,1)-1)*gridRes + xMin, (gridIdx(:,2)-1)*gridRes + yMin, (gridIdx(:,3)-1)*zRes + zMin];
end

function collision = checkCollision(pt1, pt2, occGrid, N_x, N_y, N_z)
    steps = max(abs(pt2 - pt1));
    if steps == 0
        collision = occGrid(pt1(1), pt1(2), pt1(3));
        return;
    end
    
    for s = 0:steps
        checkPt = round(pt1 + s/steps * (pt2 - pt1));
        checkPt = max([1 1 1], min([N_x N_y N_z], checkPt));
        if occGrid(checkPt(1), checkPt(2), checkPt(3))
            collision = true;
            return;
        end
    end
    collision = false;
end

function drawBuilding(x, y, height, size, color)
    half = size/2;
    verts = [x-half y-half 0; x+half y-half 0; x+half y+half 0; x-half y+half 0; ...
             x-half y-half height; x+half y-half height; x+half y+half height; x-half y+half height];
    faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
    patch('Vertices', verts, 'Faces', faces, 'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
end

function drawFire(x, y, base_height, flame_height)
    [X, Y, Z] = cylinder([1 0.2], 20);
    Z = Z * flame_height + base_height;
    surf(X*5+x, Y*5+y, Z, 'FaceColor', [1 0.3 0], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
end

function drawVegetationFire(x, y, size)
    half = size/2;
    patch([x-half x+half x+half x-half], [y-half y-half y+half y+half], ...
          [0.1 0.1 0.1 0.1], [1 0.5 0], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    [X, Y, Z] = cylinder([0.5 0.1], 20);
    Z = Z * 3 + 0.1;
    surf(X*size/2+x, Y*size/2+y, Z, 'FaceColor', [1 0.3 0], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
end

function drawBurnedGround(x, y, size)
    half = size/2;
    patch([x-half x+half x+half x-half], [y-half y-half y+half y+half], ...
          [0.05 0.05 0.05 0.05], [0.15 0.15 0.15], 'EdgeColor', 'none');
end

function out = clamp(val, minVal, maxVal)
    out = max(minVal, min(maxVal, val));
end

function drawUAV(x, y, z)
    x = double(x);
    y = double(y);
    z = double(z);
    body_size = 3;
    [X, Y, Z_cyl] = cylinder(1, 20);
    Z_cyl = Z_cyl * 0.5 + z;
    surf(X*body_size+x, Y*body_size+y, Z_cyl, 'FaceColor', [0.2 0.2 0.8], 'EdgeColor', 'none');
    rotor_positions = [body_size*1.2 0; -body_size*1.2 0; 0 body_size*1.2; 0 -body_size*1.2];
    for i = 1:4
        [X_r, Y_r, Z_r] = cylinder(0.8, 20);
        Z_r = Z_r * 0.1 + z + 0.5;
        surf(X_r*2+x+rotor_positions(i,1), Y_r*2+y+rotor_positions(i,2), Z_r, ...
             'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    end
end


