function [lat_interp, lon_interp, depth_path] = interpolate_along_bathymetry(segment_data, lat1, lon1, lat2, lon2, juld_query, juld_start, juld_end,PARAM, varargin)

%------------------------------------------------------------------------
%INTERPOLATE_ALONG_BATHYMETRY Interpolation following bathymetry using
% etopo2
%
% This function interpolates between two positions following
% bathymetric contours and produce a smooth path
%
% SYNTAX:
%   [lat_interp, lon_interp] = interpolate_along_bathymetry(lat1, lon1, lat2, lon2, juld_query, juld_start, juld_end,PARAM)
%   [lat_interp, lon_interp, depth_path] = interpolate_along_bathymetry(..., 'param', value, ...)
%
% INPUTS:
%   lat1, lon1   : Starting position (degrees)
%   lat2, lon2   : Ending position (degrees)
%   juld_query   : Dates where to calculate positions
%   juld_start   : Start date
%   juld_end     : End date
%   PARAM        : global parameters
% OPTIONAL PARAMETERS:
%   'method'           : 'isobath', 'adaptive_isobath' (default: 'isobath')
%   'target_depth'     : Target depth in meters (default: mean depth)
%   'depth_tolerance'  : Depth tolerance in meters (default: 200m)
%   'n_waypoint'       : Number of waypoints (default: 20)
%   'usegrounded'      : 0 if grounded are not used to determine initial path, 1 otherwise
% OUTPUTS:
%   lat_interp   : Interpolated latitudes
%   lon_interp   : Interpolated longitudes
%   depth_path   : Depths along the path (optional)
%
% AUTHORS  : ccabanes 01/2026
% ------------------------------------------------------------------------

% Input validation and default parameters

 n=length(varargin);

if n/2~=floor(n/2)
    error('Check the imput arguments')
end


f=varargin(1:2:end);
c=varargin(2:2:end);
s = cell2struct(c,f,2);

% default CONFIG
%params.method          = 'adaptive_isobath';
params.target_depth    = [];
params.depth_tolerance = 50;
params.n_waypoints     = 15;
params.adaptation = 0;
params.usegrounded = 0;
% Input CONFIG
if isfield(s,'method')==1;params.method=s.method;end;
if isfield(s,'target_depth')==1;params.target_depth=s.target_depth;end;
if isfield(s,'depth_tolerance')==1;params.depth_tolerance=s.depth_tolerance;end;
if isfield(s,'n_waypoints')==1;params.n_waypoints=s.n_waypoints;end;
if isfield(s,'usegrounded')==1;params.usegrounded=s.usegrounded;end;


% Compute bathymetric grid
 [lat_grid, lon_grid, depth_grid] = get_bathymetry_grid(lat1, lon1, lat2, lon2,PARAM);


% GEBCO_FILE = '/home/ccabanes/provisoire/_ressources/GEBCO_2020/GEBCO_2020.nc';
% 
%  [depth_grid, lon_grid, lat_grid] = get_gebco_elev_zone(lon1, lon2, lat1, lat2, GEBCO_FILE);
%   
% Determine target depth
if isempty(params.target_depth)
    % Target depth = mean of start and end point depths
    depth1 = interp2(lon_grid, lat_grid, depth_grid, lon1, lat1, 'linear', NaN);
    depth2 = interp2(lon_grid, lat_grid, depth_grid, lon2, lat2, 'linear', NaN);

    if isnan(depth1) || isnan(depth2)
        warning('Unable to interpolate depth at start/end points');
%         [lat_interp, lon_interp] = interpolate_geopos_sphere(lat1, lon1, lat2, lon2, juld_query, juld_start, juld_end);
        [lon_interp, lat_interp] = interpolate_between_2_locations_geo(juld_start, lon1, lat1, juld_end, lon2, lat2, juld_query);

        depth_path = [];
        return;
    end

    params.target_depth = (depth1 + depth2) / 2;
    %params.target_depth = -2500;

end

if params.n_waypoints>4

    fprintf('Bathymetric interpolation: method %s, target depth %.0f m\n', params.method, params.target_depth);


    switch params.method
        case 'isobath'
            [waypoints_lat, waypoints_lon] = find_isobath_path(lat1, lon1, lat2, lon2, lat_grid, lon_grid, depth_grid, params);

        case 'adaptive_isobath'
            [waypoints_lat, waypoints_lon] = find_adaptive_isobath_path(segment_data, lat1, lon1, lat2, lon2, lat_grid, lon_grid, depth_grid, params);

        otherwise
            error('Unknown method: %s', params.method);
    end
else
    fprintf('n_waypoints<=4 (short path), linear interpolation used');
%     [lat_interp, lon_interp] = interpolate_geopos_sphere(lat1, lon1, lat2, lon2, juld_query, juld_start, juld_end);
    [lon_interp, lat_interp] = interpolate_between_2_locations_geo(juld_start, lon1, lat1, juld_end, lon2, lat2, juld_query);

    depth_path = [];
    return;
end

% Check that a valid path was found
if isempty(waypoints_lat) || length(waypoints_lat) < 2
    fprintf('Unable to find a bathymetric path. Linear interpolation used.');
%     [lat_interp, lon_interp] = interpolate_geopos_sphere(lat1, lon1, lat2, lon2, juld_query, juld_start, juld_end);
    [lon_interp, lat_interp] = interpolate_between_2_locations_geo(juld_start, lon1, lat1, juld_end, lon2, lat2, juld_query);
    depth_path = [];
    return;
end

% Interpolation along the found path
[lat_interp, lon_interp] = interpolate_along_waypoints(waypoints_lat, waypoints_lon, juld_query, juld_start, juld_end);

% Compute depths along the path 
if nargout > 2
    depth_path = zeros(size(lat_interp));
    for i = 1:length(lat_interp)
        depth_path(i) = interp2(lon_grid, lat_grid, depth_grid, lon_interp(i), lat_interp(i), 'linear', NaN);
    end
end

end


% METHOD 1: PATH FOLLOWING AN ISOBATH at target depth
% ========================================================================
%------------------------------------------------------------------------
function [waypoints_lat, waypoints_lon] = find_isobath_path(lat1, lon1, lat2, lon2, lat_grid, lon_grid, depth_grid, params)
% FIND_ISOBATH_PATH Finds a path following an isobath
%------------------------------------------------------------------------

fprintf('  Searching for a path along the %.0f m isobath...\n', params.target_depth);

% Create a finer grid (*2) for contour search
[lon_fine, lat_fine] = meshgrid(linspace(min(lon_grid(:)), max(lon_grid(:)), size(lon_grid, 2)*1), ...
    linspace(min(lat_grid(:)), max(lat_grid(:)), size(lat_grid, 1)*1));
depth_fine = interp2(lon_grid, lat_grid, depth_grid, lon_fine, lat_fine, 'spline');


% Extract contours at target depth
target_depths = [params.target_depth - params.depth_tolerance, ...
    params.target_depth, ...
    params.target_depth + params.depth_tolerance];

best_path_lat = [];
best_path_lon = [];
min_distance = inf;

for target_dep = target_depths
    try
        % Extract contour
        %figure(2)
        %C = contour(lon_fine, lat_fine, depth_fine, [target_dep target_dep]);
       C = contourc(lon_grid(1,:), lat_grid(:,1), depth_grid,  [target_dep target_dep]);

        if size(C, 2) < 3
            continue;
        end

        % Analyze extracted contours
        [contour_segments] = parse_contour_matrix(C);

        % Find the contour closest to start/end points
        for i = 1:length(contour_segments)
            segment_lat = contour_segments{i}(:, 2);
            segment_lon = contour_segments{i}(:, 1);

            if length(segment_lat) < 3
                continue;
            end

            % Compute distance to start and end points
            dist_start = min(sqrt((segment_lat - lat1).^2 + (segment_lon - lon1).^2));
            dist_end = min(sqrt((segment_lat - lat2).^2 + (segment_lon - lon2).^2));
            total_dist = dist_start + dist_end;

            if total_dist < min_distance
                min_distance = total_dist;
                best_path_lat = segment_lat;
                best_path_lon = segment_lon;
            end
        end

    catch ME
        fprintf('    Error while extracting contour %.0f m: %s\n', target_dep, ME.message);
        continue;
    end
end

% If a contour was found, create the waypoints
if ~isempty(best_path_lat)
    % Connect start point to nearest contour segment
    [~, idx_start] = min(sqrt((best_path_lat - lat1).^2 + (best_path_lon - lon1).^2));
    [~, idx_end] = min(sqrt((best_path_lat - lat2).^2 + (best_path_lon - lon2).^2));

    if idx_start >= idx_end
        best_path_lat=flipud(best_path_lat);
        best_path_lon=flipud(best_path_lon);
    end

    [~, idx_start] = min(sqrt((best_path_lat - lat1).^2 + (best_path_lon - lon1).^2));
    [~, idx_end] = min(sqrt((best_path_lat - lat2).^2 + (best_path_lon - lon2).^2));

    % Create the path following the contour
    contour_lat = best_path_lat(idx_start:idx_end);
    contour_lon = best_path_lon(idx_start:idx_end);

    % Complete path
    waypoints_lat = [lat1; contour_lat; lat2];
    waypoints_lon = [lon1; contour_lon; lon2];


    % --- Compute cumulative distance along the path  ---
    distances = distance_lpo(waypoints_lat, waypoints_lon)';
    cumdist = [0; cumsum(distances(:))];

    % Remove duplicates or zero-length segments
    valid_idx = [true; diff(cumdist) > 1e-6];
    cumdist = cumdist(valid_idx);
    waypoints_lat = waypoints_lat(valid_idx);
    waypoints_lon = waypoints_lon(valid_idx);

    % --- Define target distances for uniform spacing ---
    target_dist = linspace(0, cumdist(end), params.n_waypoints+2);

    % --- Interpolate lat/lon positions at uniform distances ---
    waypoints_lat = interp1(cumdist, waypoints_lat, target_dist, 'linear');
    waypoints_lon = interp1(cumdist, waypoints_lon, target_dist, 'linear');

    % Optional: smooth the path to avoid abrupt changes
    if length(waypoints_lat) > 3
        try
            % Simple moving average smoothing
            smooth_window = min(10, floor(length(waypoints_lat)/3));
            if smooth_window >= 3
                waypoints_lat(2:end-1) = smooth(waypoints_lat(2:end-1), smooth_window);
                waypoints_lon(2:end-1) = smooth(waypoints_lon(2:end-1), smooth_window);
            end
            fprintf('    Path smoothing - Simple moving average \n');
    catch
        % If smoothing fails, keep original points
        fprintf('    Path smoothing failed, keeping original points\n');
    end
end

    fprintf('    Path found with %d waypoints\n', length(waypoints_lat));
else
    fprintf('    No path found along the isobath\n');
    waypoints_lat = [];
    waypoints_lon = [];
end

end


% ========================================================================
% METHOD 2: PATH FOLLOWING AN ISOBATH THAT ADAPTS BETWEEN START AND END POINTS
% ========================================================================
%------------------------------------------------------------------------
function [waypoints_lat, waypoints_lon] = find_adaptive_isobath_path(segment_data, lat1, lon1, lat2, lon2, lat_grid, lon_grid, depth_grid, params)
% Creates a trajectory that progressively follows changing
% isobaths between the start and end points
%------------------------------------------------------------------------

fprintf('  Searching for an adaptive path between isobaths...\n');

% Compute depths at start and end points
depth_start = interp2(lon_grid, lat_grid, depth_grid, lon1, lat1, 'linear', NaN);
depth_end = interp2(lon_grid, lat_grid, depth_grid, lon2, lat2, 'linear', NaN);

if isnan(depth_start) || isnan(depth_end)
    fprintf('    Unable to retrieve depths at start/end points\n');
    waypoints_lat = [];
    waypoints_lon = [];
    return;
end

fprintf('    Start depth: %.0f m, End depth: %.0f m\n', depth_start, depth_end);

% find subsegment for adaptation of the isobath from depth_start eo depth_end (take also into account 
% grounded depths)
[grp_start, grp_end] = findContiguousOnes(segment_data.grounded);
[ind_max, cycle_max,profPresMax] = maxDepthPerGroundedGroup(segment_data.grounded, segment_data.cycleNumber, segment_data.profPresMax);
[subseg_start, subseg_end] = buildSubSegmentsFromIndMax(ind_max, length(segment_data.cycleNumber), segment_data.profPresMax,10);

% Number of segments for adaptation of the isobath
% n_segments = round(max(1, params.n_waypoints/10));
% segment_ratio = linspace(0, 1, n_segments + 1);
n_segments =length(subseg_start);
% Initialize waypoints
all_waypoints_lat = lat1;
all_waypoints_lon = lon1;

% Current position
current_lat = lat1;
current_lon = lon1;

% first guess path (test)
[first_waypoints_lat, first_waypoints_lon] = find_isobath_path(lat1, lon1, lat2, lon2, lat_grid, lon_grid, depth_grid, params);
% Re-sample first guess path to match segment ratio points
% ref_indices = round(linspace(1, length(first_waypoints_lat), n_segments + 1));
% first_waypoints_lon = first_waypoints_lon(ref_indices);
% first_waypoints_lat = first_waypoints_lat(ref_indices);
first_waypoints_lon = first_waypoints_lon([subseg_start(1) subseg_end]);
first_waypoints_lat = first_waypoints_lat([subseg_start(1) subseg_end]);


target_depth_segment=-segment_data.profPresMax(subseg_end);
target_depth_segment(segment_data.grounded(subseg_end)==0) = NaN;
if params.usegrounded==0  % do not use grounded depths (=> will interpolate instead)
    target_depth_segment(segment_data.grounded(subseg_end)==1) = NaN;
end
if isnan(target_depth_segment(end))
    target_depth_segment(end) = depth_end;
end



x = [segment_data.time(1) segment_data.time(subseg_end)];
y = [depth_start target_depth_segment];

% linear interp
target_depth_interp = interp1( ...
    x(~isnan(y)), y(~isnan(y)), ...
    x, 'linear');

%close all
%figure
lat_max=max(lat1,lat2); lat_min=min(lat1,lat2);
lon_max=max(lon1,lon2); lon_min=min(lon1,lon2);
m_proj('stereographic','lat',min((lat_max+lat_min)/2+2,90),'lon',max((lon_min+lon_max)/2-2,-180),'radius',min(lat_max-lat_min+5,30));
[depth_map,lon_map,lat_map] = m_tbase([-180 180 70 90]); 
%m_contour(lon_map,lat_map,depth_map,[-2000 -1000 0]);
%hold on
%m_grid 

if ~isfield(params,'adaptation')
        params.adaptation = 0;   % default = standard behavior
end
if n_segments<3
   params.adaptation = 0;
end

%m_plot(first_waypoints_lon,first_waypoints_lat,'r')
% 
% Create each segment with a progressively changing target depth
for i = 1:n_segments

        if i==1
            target_depth_segment = double((target_depth_interp(i+1)+target_depth_interp(i)))/2; % smooth start
        else
            target_depth_segment = double((target_depth_interp(i+1)));
        end   
        % Compute target depth for this segment (progressive interpolation)
        % target_depth_segment = depth_start + (depth_end - depth_start) * segment_ratio(i+1);
%         m_contour(lon_map,lat_map,depth_map,[depth_start depth_start],'c')
% 
%         m_contour(lon_map,lat_map,depth_map,[target_depth_segment target_depth_segment]);
%         m_etopo2('contour',[target_depth_segment target_depth_segment])
    % Compute target geographic position for this segment from the first
    % guess path (test)
    target_lat_segment = first_waypoints_lat(i+1);
    target_lon_segment = first_waypoints_lon(i+1);

%     % Compute target geographic position for this segment
%         target_lat_segment = lat1 + (lat2 - lat1) * segment_ratio(i+1);
%         target_lon_segment = lon1 + (lon2 - lon1) * segment_ratio(i+1);
% m_plot(current_lon,current_lat,'r^') 
% m_plot(lon2,lat2,'vr') 
% m_plot(target_lon_segment,target_lat_segment,'+r') 

    %     % Handle -180/180 crossing  %  cc pas ici
%     dlon = target_lon_segment - current_lon;
%     if abs(dlon) > 180
%         if dlon > 0
%             current_lon = current_lon + 360;
%         else
%             target_lon_segment = target_lon_segment + 360;
%         end
%     end

    fprintf('    Segment %d/%d: target depth %.0f m\n', i, n_segments, target_depth_segment);
    % Search for contours near the target depth
    search_depths = [target_depth_segment - params.depth_tolerance/2, ...
        target_depth_segment, ...
        target_depth_segment + params.depth_tolerance/2];

    segment_waypoints_lat = [];
    segment_waypoints_lon = [];
    min_path_length = inf;

    % Test different depths for this segment
    for test_depth = search_depths
        try
            % Extract contour
            C = contourc(lon_grid(1,:), lat_grid(:,1), depth_grid, [test_depth test_depth]);

            if size(C, 2) < 3
                continue;
            end

            % Analyze extracted contours
            [contour_segments] = parse_contour_matrix(C);

            % Find the best contour segment for this path section
            for j = 1:length(contour_segments)
                segment_lat = contour_segments{j}(:, 2);
                segment_lon = contour_segments{j}(:, 1);

                if length(segment_lat) < 3
                    continue;
                end

                % Compute distances to segment start and end point

                 dist_from_current_all = arrayfun(@(la,lo) ...
                     distance_lpo([current_lat la], [current_lon lo]), ...
                     segment_lat, segment_lon);
                [dist_from_current,idx_start] = min(dist_from_current_all);

                dist_to_target_all = arrayfun(@(la,lo) ...
                    distance_lpo([target_lat_segment la], [target_lon_segment lo]), ...
                    segment_lat, segment_lon);


                [dist_to_target,idx_end] = min(dist_to_target_all);


                % Selection best contour
                path_length = dist_from_current + dist_to_target;
               % if path_length < min_path_length && dist_from_current < 100000 && dist_to_target < 100000
                if path_length < min_path_length 
                    min_path_length = path_length;

                    % Find connection points on this contour
                    % idx_start and idx_end, already computed
%                    m_plot(segment_lon,segment_lat,'+y')

                    % Extract the appropriate sub-segment
                    if idx_start <= idx_end
                        selected_lat = segment_lat(idx_start:idx_end);
                        selected_lon = segment_lon(idx_start:idx_end);
                    else
                        selected_lat = segment_lat(idx_start:-1:idx_end);
                        selected_lon = segment_lon(idx_start:-1:idx_end);
                    end

                    segment_waypoints_lat = selected_lat;
                    segment_waypoints_lon = selected_lon;
                end
            end

        catch ME
            fprintf('      Contour error %.0f m: %s\n', test_depth, ME.message);
            continue;
        end
    end

    % If no contour found, use direct interpolation
    if isempty(segment_waypoints_lat)
        fprintf('      No contour found, using direct interpolation\n');
        % Direct interpolation toward target position
        n_direct_points = 3;
        direct_ratio = linspace(0, 1, n_direct_points + 1)';
        direct_ratio = direct_ratio(2:end); % Exclude starting point (already included)

        segment_waypoints_lat = current_lat + (target_lat_segment - current_lat) * direct_ratio;
        segment_waypoints_lon = current_lon + (target_lon_segment - current_lon) * direct_ratio;
    end

    % Add waypoints from this segment (without duplicating the first point)
    if length(segment_waypoints_lat) > 1
        all_waypoints_lat = [all_waypoints_lat; segment_waypoints_lat(2:end)];
        all_waypoints_lon = [all_waypoints_lon; segment_waypoints_lon(2:end)];

        % Update current position
        current_lat = segment_waypoints_lat(end);
        current_lon = segment_waypoints_lon(end);
    else
        % If only one point, use it as the new current position
        current_lat = target_lat_segment;
        current_lon = target_lon_segment;
        all_waypoints_lat = [all_waypoints_lat; current_lat];
        all_waypoints_lon = [all_waypoints_lon; current_lon];
    end
% m_plot(segment_waypoints_lon,segment_waypoints_lat,'m+-')
    
end

% Ensure the final point matches the destination
if abs(all_waypoints_lat(end) - lat2) > 1e-6 || abs(all_waypoints_lon(end) - lon2) > 1e-6
    all_waypoints_lat = [all_waypoints_lat; lat2];
    all_waypoints_lon = [all_waypoints_lon; lon2];
end

% Normalize longitudes %cc pas ici
all_waypoints_lon = mod(all_waypoints_lon + 180, 360) - 180;

% --- Compute cumulative distance along the path  ---
distances = distance_lpo(all_waypoints_lat', all_waypoints_lon');  

cumdist = [0; cumsum(distances(:))];

% Remove duplicates or zero-length segments
valid_idx = [true; diff(cumdist) > 1e-6];
cumdist = cumdist(valid_idx);
all_waypoints_lat = all_waypoints_lat(valid_idx);
all_waypoints_lon = all_waypoints_lon(valid_idx);

% --- Define target distances for uniform spacing ---
target_dist = linspace(0, cumdist(end), params.n_waypoints+2);

% --- Interpolate lat/lon positions at uniform distances ---
waypoints_lat = interp1(cumdist, all_waypoints_lat, target_dist, 'linear');
waypoints_lon = interp1(cumdist, all_waypoints_lon, target_dist, 'linear');


fprintf('    Adaptive path created with %d waypoints\n', length(waypoints_lat));

% Check path consistency
if length(waypoints_lat) < 2
    fprintf('    Insufficient path, linear interpolation required\n');
    waypoints_lat = [];
    waypoints_lon = [];
    return;
end

% Optional: smooth the path to avoid abrupt changes
if length(waypoints_lat) > 3
    try
        % Simple moving average smoothing
        smooth_window = min(10, floor(length(waypoints_lat)/3));
        if smooth_window >= 3
            waypoints_lat(2:end-1) = smooth(waypoints_lat(2:end-1), smooth_window);
            waypoints_lon(2:end-1) = smooth(waypoints_lon(2:end-1), smooth_window);
        end
        fprintf('    Path smoothing - Simple moving average \n');
    catch
        % If smoothing fails, keep original points
        fprintf('    Path smoothing failed, keeping original points\n');
    end
end

% m_plot(waypoints_lon,waypoints_lat,'k+-')
%close(100)

end

% ========================================================================
% OTHER FUNCTIONS
% ========================================================================

% ===================== PARSE CONTOUR MATRIX =========================
%------------------------------------------------------------------------
function segments = parse_contour_matrix(C)
%------------------------------------------------------------------------

segments = {};
k = 1;
while k < size(C,2)
    npts = C(2,k);
    pts = C(:, k+1:k+npts);
    segments{end+1} = pts';
    k = k + npts + 1;
end

end

% ===================== GET BATHYMETRY GRID MT_BASE or ETOPO2 =========================
%------------------------------------------------------------------------
function [lat_grid, lon_grid, depth_grid] = get_bathymetry_grid(lat1, lon1, lat2, lon2, PARAM)
%   Retrieves a continuous bathymetry grid
%   detects the longitude convention: [-360,0], [-180,180], [0,360]
%   and returns the grid in the same convention as the input.
%------------------------------------------------------------------------

marge = 15;

% === 1. Automatic detection of input longitude convention ===
if all([lon1, lon2] <= 0) && all([lon1, lon2] >= -360)
    normType = -360;   % => Convention [-360, 0]
elseif all([lon1, lon2] >= 0) && all([lon1, lon2] <= 360)
    normType = 360;    % => Convention [0, 360]
else
    normType = 180;    % => Convention [-180, 180]
end

% === 2. Internal normalization to [-180,180] for processing ===
lon1n = mod(lon1 + 180, 360) - 180;
lon2n = mod(lon2 + 180, 360) - 180;

% === 3. Latitude handling ===
lat_min = max(-90, min(lat1, lat2) - marge);
lat_max = min(90, max(lat1, lat2) + marge);

% === 4. Handling crossing of the dateline ===
cross_dateline = abs(lon2n - lon1n) > 180;
if cross_dateline
    if lon1n < 0
        lon1n = lon1n + 360;
    else
        lon2n = lon2n + 360;
    end
end

lon_min = min(lon1n, lon2n) - marge;
lon_max = max(lon1n, lon2n) + marge;

fprintf('  Retrieving bathymetry (%.1f°x%.1f°)...\n', lon_max - lon_min, lat_max - lat_min);

try
    % === 5. Retrieval via m_tbase ===

    if PARAM.isarctic==0; % nb: initialization of m_proj has not impact on the bathy grid: useless
        m_proj( 'mercator','longitude',[lon_min lon_max],'latitude',[lat_min lat_max]);
        %[depth_grid,lon_grid,lat_grid] = m_tbase([lon_min lon_max lat_min lat_max]);
        [depth_grid,lon_grid,lat_grid] = m_etopo2([lon_min lon_max lat_min lat_max]);
    else
        m_proj('stereographic','lat',min((lat_max+lat_min)/2+2,90),'lon',max((lon_min+lon_max)/2-2,-180),'radius',min(lat_max-lat_min+2,30));
       %[depth_grid,lon_grid,lat_grid] = m_tbase([lon_min lon_max lat_min lat_max]);
       [depth_grid,lon_grid,lat_grid] = m_etopo2([lon_min lon_max lat_min lat_max]);
    end


    % === 6. Sorting and continuity ===
    lon_row = lon_grid(1, :);
    [lon_sorted, sidx] = sort(lon_row);
    depth_sorted = depth_grid(:, sidx);

    % Detect wrapping
    d = diff(lon_sorted);
    wrap_gap = (lon_sorted(1) + 360) - lon_sorted(end);
    gaps = [d, wrap_gap];
    [~, idx_gap] = max(gaps);
    start_idx = mod(idx_gap, numel(lon_sorted)) + 1;

    % Reorder
    N = numel(lon_sorted);
    new_order = [start_idx:N, 1:start_idx-1];
    lon_cont = lon_sorted(new_order);
    depth_cont = depth_sorted(:, new_order);

    % Ensure continuous increase
    ref0 = lon_cont(1);
    for k = 1:N
        if lon_cont(k) < ref0
            lon_cont(k) = lon_cont(k) + 360;
        end
    end

    lon_grid = repmat(lon_cont, size(depth_cont, 1), 1);
    depth_grid = depth_cont;

%     % === 7. Final conversion according to input convention ===
%     switch normType
%         case -360
%             % Convert to [-360,0]
%             lon_grid = mod(lon_grid, 360);
%             lon_grid(lon_grid > 0) = lon_grid(lon_grid > 0) - 360;
% 
%         case 360
%             % Convert to [0,360]
%             lon_grid = mod(lon_grid, 360);
% 
%         otherwise  % normType == 180
%             % Convert to [-180,180]
%             lon_grid = mod(lon_grid + 180, 360) - 180;
%     end
% 
%     % Reorder columns according to increasing longitudes
%     [~, idx_final] = sort(lon_grid(1, :));
%     lon_grid = lon_grid(:, idx_final);
%     depth_grid = depth_grid(:, idx_final);

    

catch ME
    fprintf('  Error while retrieving bathymetry: %s\n', ME.message);
    lat_grid = [];
    lon_grid = [];
    depth_grid = [];
end

end


%------------------------------------------------------------------------
function [lat_interp, lon_interp] = interpolate_along_waypoints(waypoints_lat, waypoints_lon, juld_query, juld_start, juld_end)
%------------------------------------------------------------------------
% Interpolates along waypoints

dist_segments = distance_lpo(waypoints_lat, waypoints_lon);  % returns [n-1 x 1]

% --- Compute cumulative distance along path ---
distances = [0, cumsum(dist_segments)];

% --- Normalize distances to [0, 1] ---
if distances(end) == 0
    distances_norm = distances;
else
    distances_norm = distances / distances(end);
end


% Temporal interpolation factor
t = (juld_query - juld_start) / (juld_end - juld_start);
t = max(0, min(1, t)); 

% Interpolation along the path
lat_interp = interp1(distances_norm, waypoints_lat, t);
lon_interp = interp1(distances_norm, waypoints_lon, t);
end



%--------------------------------------------------------------

function [o_interpLon, o_interpLat] = interpolate_between_2_locations_geo( ...
    t1, lon1, lat1, t2, lon2, lat2, tInterp)
% ------------------------------------------------------------
% Interpolate between two geographic locations using a geodesic
% (WGS84 ellipsoid) instead of linear lat/lon interpolation.
% ------------------------------------------------------------
if exist('distance','file')&exist('reckon','file')&exist('azimuth','file')
     [o_interpLon, o_interpLat] = interpolate_between_2_locations_geo1( ...
    t1, lon1, lat1, t2, lon2, lat2, tInterp);
else
    [o_interpLon, o_interpLat] = interpolate_between_2_locations_geo2( ...
    t1, lon1, lat1, t2, lon2, lat2, tInterp);
end

end

function [o_interpLon, o_interpLat] = interpolate_between_2_locations_geo1( ...
    t1, lon1, lat1, t2, lon2, lat2, tInterp)
% ------------------------------------------------------------
% Interpolate between two geographic locations using a geodesic
% (WGS84 ellipsoid) instead of linear lat/lon interpolation.
% use mapping toolbox
% ------------------------------------------------------------

% Output initialization
o_interpLon = [];
o_interpLat = [];


% Fraction of time between the two fixes
f = (tInterp - t1) / (t2 - t1);

% If time exactly matches one endpoint
if f <= 0
    o_interpLon = lon1;
    o_interpLat = lat1;
    return
elseif f >= 1
    o_interpLon = lon2;
    o_interpLat = lat2;
    return
end

% ------------------------------------------------------------
% Geodesic calculation using WGS84 ellipsoid
% ------------------------------------------------------------

% Compute geodesic distance and azimuth from point 1 to point 2
ellipsoid = wgs84Ellipsoid('kilometer');

[dist12, az12] = distance(lat1, lon1, lat2, lon2, ellipsoid); % distance in km

% Distance to travel along the geodesic
dInterp = f * dist12;

% Move from point 1 along the geodesic
[o_interpLat, o_interpLon] = reckon(lat1, lon1, dInterp, az12, ellipsoid);


end


function [o_interpLon, o_interpLat] = interpolate_between_2_locations_geo2( ...
    t1, lon1, lat1, t2, lon2, lat2, tInterp)
% ------------------------------------------------------------
% Interpolate between two geographic locations using a geodesic
% (WGS84 ellipsoid) instead of linear lat/lon interpolation.
% does not require mapping toolbox
% use distance_lpo
% ------------------------------------------------------------

% Ensure column vector
tInterp = tInterp(:);

% Fraction of time
f = (tInterp - t1) / (t2 - t1);

% Clamp
f = max(0, min(1, f));

% ------------------------------------------------------------
% Choose N based on temporal sampling
% ------------------------------------------------------------
Nt = numel(tInterp);
oversample = 3;                 % user-tunable
N = max(2, oversample * Nt);

% ------------------------------------------------------------
% Compute geodesic
% ------------------------------------------------------------
[~, lat_gd, lon_gd] = distance_lpo( ...
    [lat1 lat2], [lon1 lon2], N, 'wgs84');

% ------------------------------------------------------------
% Map time fraction to geodesic index
%------------------------------------------------------------
idx = round(1 + f * (N-1));
idx = max(1, min(N, idx));

o_interpLat = lat_gd(idx);
o_interpLon = lon_gd(idx);

end

function [lat_interp, lon_interp] = interpolate_geopos_sphere(lat1, lon1, lat2, lon2, juld_query, juld_start, juld_end)

%INTERPOLATE_GEOPOS Linear interpolation between two (lon, lat) points
%  Haversine formulae
%
% Inputs:
%   lat1, lon1 : Starting point (in degrees)
%   lat2, lon2 : Ending point (in degrees)
%   juld_query : Dates (Julian days) at which positions are required
%   juld_start : Start date of the trajectory
%   juld_end   : End date of the trajectory
%
% Outputs:
%   lat_interp, lon_interp : Interpolated positions (in degrees)
% ------------------------------------------------------------------------

% Approximate Earth radius in km
R = 6371;

% Conversion to radians
lat1_rad = deg2rad(lat1); lon1_rad = deg2rad(lon1);
lat2_rad = deg2rad(lat2); lon2_rad = deg2rad(lon2);

% Handle longitudes to avoid ±180° discontinuity
% If the longitude difference > 180°, adjust longitudes
dlon = lon2_rad - lon1_rad;
if abs(dlon) > pi
    if dlon > 0
        lon1_rad = lon1_rad + 2*pi;
    else
        lon2_rad = lon2_rad + 2*pi;
    end
end

% 3D vectors on the sphere
p1 = R * [cos(lat1_rad)*cos(lon1_rad), cos(lat1_rad)*sin(lon1_rad), sin(lat1_rad)];
p2 = R * [cos(lat2_rad)*cos(lon2_rad), cos(lat2_rad)*sin(lon2_rad), sin(lat2_rad)];

% Angle between the two vectors (in radians)
omega = acos(dot(p1, p2) / (R^2));

% Fraction of the trajectory for each date
f = (juld_query - juld_start) / (juld_end - juld_start);

lat_interp = zeros(size(f));
lon_interp = zeros(size(f));

for i = 1:length(f)
    if omega == 0
        p = p1;
    else
        % Spherical linear interpolation
        A = sin((1 - f(i)) * omega) / sin(omega);
        B = sin(f(i) * omega) / sin(omega);
        p = A * p1 + B * p2;
    end

    % Convert back to lat/lon
    x = p(1); y = p(2); z = p(3);
    lat_interp(i) = rad2deg(asin(z / R));
    lon_interp(i) = rad2deg(atan2(y, x));
end

% Final conversion of longitudes into [-180, 180]
%lon_interp = mod(lon_interp + 180, 360) - 180;
end

function [grp_start, grp_end] = findContiguousOnes(isgrounded)
% FINDCONTIGUOUSONES
% Finds contiguous groups of 1s in a binary vector
%
% INPUT:
%   isgrounded : vector containing 0 or 1
%
% OUTPUT:
%   grp_start  : indices where a group of 1s starts
%   grp_end    : indices where a group of 1s ends

isgrounded = isgrounded(:)';   % ensure row vector

% Find transitions
d = diff([0 isgrounded 0]);

% Start of group: 0 -> 1 transition
grp_start = find(d == 1);

% End of group: 1 -> 0 transition
grp_end = find(d == -1) - 1;

end




function [ind_max, cycle_max, profPresMax_max] = ...
    maxDepthPerGroundedGroup(isgrounded, cycleNumber, profPresMax)
%
% For each contiguous group of grounded==1, find the point where
% profPresMax is minimun.
%
% INPUTS:
%   isgrounded   - vector (0/1)
%   cycleNumber  - vector of cycle numbers
%   profPresMax  - vector of maximum depths
%
% OUTPUTS:
%   ind_max            - index of maximum depth in each group
%   cycle_max          - cycle number corresponding to ind_max
%   profPresMax_max    - maximum profPresMax value per group

% Find contiguous grounded groups
[grp_start, grp_end] = findContiguousOnes(isgrounded);

nGroups = numel(grp_start);

% Preallocate
ind_max = nan(1, nGroups);
cycle_max = nan(1, nGroups);
profPresMax_max = nan(1, nGroups);

for k = 1:nGroups
    % Indices of the current group
    idx = grp_start(k):grp_end(k);
    
    % Minimum pressure inside the group
    [profPresMax_max(k), imax] = min(profPresMax(idx));
    
    % Absolute index
    ind_max(k) = idx(imax);
    
    % Corresponding cycle number
    cycle_max(k) = cycleNumber(ind_max(k));
end

end


function [subseg_start, subseg_end] = buildSubSegmentsFromIndMax(ind_max, nPoints, ProfPresMax, min_point)
% Build sub-segments using ind_max as  segment ends
%
% INPUTS:
%   ind_max    - indices defining the end of each segment
%   nPoints    - total number of points in the vector
%   min_point - minimum number of points per sub-segment
%
% OUTPUTS:
%   subseg_start - indices of sub-segment starts
%   subseg_end   - indices of sub-segment ends

% Ensure row vector and sorted
ind_max = sort(ind_max(:)');

% Remove invalid or duplicate indices
ind_max = ind_max(ind_max >= 1 & ind_max <= nPoints);
ind_max = unique(ind_max);

if isempty(ind_max)
    ind_max=nPoints;
end
if ind_max(1)~=1  
ind_max=[1 ind_max];
end
if ind_max(end)~=nPoints  
ind_max=[ind_max nPoints];
end

% Initialisation avec le premier point
ind_filt = ind_max(1);
% Parcours des points intermédiaires (sauf le dernier)
for k = 2:length(ind_max)-1
    if ind_max(k) - ind_filt(end) >= min_point
        ind_filt(end+1) = ind_max(k);
    end
end
% Gestion du dernier point
if ind_max(end) - ind_filt(end) < min_point
    % Si trop proche, on supprime l'avant-dernier
    ind_filt(end) = [];
end
%ind_filt(end+1)=ind_max(end);

% enleve premeier point
if isempty(ind_filt)==0
ind_filt(1) = [];
end

subseg_start = [];
subseg_end = [];

% Define segment starts and ends
seg_start = [1, ind_filt(1:end) + 1];
seg_end   = [ind_filt nPoints] ;



% Loop on each segment
for k = 1:length(seg_start)
    
    s0 = seg_start(k);
    s1 = seg_end(k);
    seg_len = s1 - s0 + 1;
    
    if seg_len < 2*min_point
        % Segment too short → keep as is
         subseg_start(end+1) = s0;
         subseg_end(end+1)   = s1;
    else
        % Split into sub-segments
        nFull = floor(seg_len / min_point);
        
        for i = 1:nFull
            ss = s0 + (i-1)*min_point;
            se = ss + min_point - 1;
            
            % Ensure last sub-segment ends exactly at ind_max
            if i == nFull
                se = s1;
            end
            
            subseg_start(end+1) = ss;
            subseg_end(end+1)   = se;
        end
    end
end

end

function [range,A12,A21]=distance_lpo(lat,long,argu1,argu2);
% DIST    Computes distance and bearing between points on the earth
%         using various reference spheroids.
%
%         [RANGE,AF,AR]=DIST(LAT,LONG) computes the ranges RANGE between
%         points specified in the LAT and LONG vectors (decimal degrees with
%         positive indicating north/east). Forward and reverse bearings
%         (degrees) are returned in AF, AR.
%
%         [RANGE,GLAT,GLONG]=DIST(LAT,LONG,N) computes N-point geodesics
%         between successive points. Each successive geodesic occupies
%         it's own row (N>=2)
%
%         [..]=DIST(...,'ellipsoid') uses the specified ellipsoid
%         to get distances and bearing. Available ellipsoids are:
%
%         'clarke66'  Clarke 1866
%         'iau73'     IAU 1973
%         'wgs84'     WGS 1984
%         'sphere'    Sphere of radius 6371.0 km
%
%          The default is 'wgs84'.
%
%          Ellipsoid formulas are recommended for distance d<2000 km,
%          but can be used for longer distances.

%Notes: RP (WHOI) 3/Dec/91
%         Mostly copied from BDC "dist.f" routine (copied from ....?), but
%         then wildly modified to bring it in line with Matlab vectorization.
%
%       RP (WHOI) 6/Dec/91
%         Feeping Creaturism! - added geodesic computations. This turned
%         out to be pretty hairy since there were a lot of branch problems
%         with asin, atan when computing geodesics subtending > 90 degrees
%         that were ignored in the original code!
%       RP (WHOI) 15/Jan/91
%         Fixed some bothersome special cases, like when computing geodesics
%         and N=2, or LAT=0...
%	A Newhall (WHOI) Sep 1997
%	   modified and fixed a bug found in Matlab version 5
%
%		NOTE: This routine may interfere with dist that
%			is supplied with matlab's neural net toolbox.

%C GIVEN THE LATITUDES AND LONGITUDES (IN DEG.) IT ASSUMES THE IAU SPHERO
%C DEFINED IN THE NOTES ON PAGE 523 OF THE EXPLANATORY SUPPLEMENT TO THE
%C AMERICAN EPHEMERIS.
%C
%C THIS PROGRAM COMPUTES THE DISTANCE ALONG THE NORMAL
%C SECTION (IN M.) OF A SPECIFIED REFERENCE SPHEROID GIVEN
%C THE GEODETIC LATITUDES AND LONGITUDES OF THE END POINTS
%C  *** IN DECIMAL DEGREES ***
%C
%C  IT USES ROBBIN'S FORMULA, AS GIVEN BY BOMFORD, GEODESY,
%C FOURTH EDITION, P. 122.  CORRECT TO ONE PART IN 10**8
%C AT 1600 KM.  ERRORS OF 20 M AT 5000 KM.
%C
%C   CHECK:  SMITHSONIAN METEOROLOGICAL TABLES, PP. 483 AND 484,
%C GIVES LENGTHS OF ONE DEGREE OF LATITUDE AND LONGITUDE
%C AS A FUNCTION OF LATITUDE. (SO DOES THE EPHEMERIS ABOVE)
%C
%C PETER WORCESTER, AS TOLD TO BRUCE CORNUELLE...1983 MAY 27
%C

spheroid='wgs84';
geodes=0;
if (nargin >= 3),
    if (isstr(argu1)),
        spheroid=argu1;
    else
        geodes=1;
        Ngeodes=argu1;
        if (Ngeodes <2), error('Must have at least 2 points in a goedesic!');end;
        if (nargin==4), spheroid=argu2; end;
    end;
end;

if (spheroid(1:3)=='sph'),
    A = 6371000.0;
    B = A;
    E = sqrt(A*A-B*B)/A;
    EPS= E*E/(1-E*E);
elseif (spheroid(1:3)=='cla'),
    A = 6378206.4E0;
    B = 6356583.8E0;
    E= sqrt(A*A-B*B)/A;
    EPS = E*E/(1.-E*E);
elseif(spheroid(1:3)=='iau'),
    A = 6378160.e0;
    B = 6356774.516E0;
    E = sqrt(A*A-B*B)/A;
    EPS = E*E/(1.-E*E);
elseif(spheroid(1:3)=='wgs'),

    %c on 9/11/88, Peter Worcester gave me the constants for the
    %c WGS84 spheroid, and he gave A (semi-major axis), F = (A-B)/A
    %c (flattening) (where B is the semi-minor axis), and E is the
    %c eccentricity, E = ( (A**2 - B**2)**.5 )/ A
    %c the numbers from peter are: A=6378137.; 1/F = 298.257223563
    %c E = 0.081819191
    A = 6378137.;
    E = 0.081819191;
    B = sqrt(A.^2 - (A*E).^2);
    EPS= E*E/(1.-E*E);

else
    error('dist: Unknown spheroid specified!');
end;


NN=max(size(lat));
if (NN ~= max(size(long))),
    error('dist: Lat, Long vectors of different sizes!');
end

if (NN==size(lat)), rowvec=0;  % It is easier if things are column vectors,
else                rowvec=1; end; % but we have to fix things before returning!

lat=lat(:)*pi/180;     % convert to radians
long=long(:)*pi/180;

lat(lat==0)=eps*ones(sum(lat==0),1);  % Fixes some nasty 0/0 cases in the
% geodesics stuff

PHI1=lat(1:NN-1);    % endpoints of each segment
XLAM1=long(1:NN-1);
PHI2=lat(2:NN);
XLAM2=long(2:NN);

% wiggle lines of constant lat to prevent numerical probs.
if (any(PHI1==PHI2)),
    for ii=1:NN-1,
        if (PHI1(ii)==PHI2(ii)), PHI2(ii)=PHI2(ii)+ 1e-14; end;
    end;
end;
% wiggle lines of constant long to prevent numerical probs.
if (any(XLAM1==XLAM2)),
    for ii=1:NN-1,
        if (XLAM1(ii)==XLAM2(ii)), XLAM2(ii)=XLAM2(ii)+ 1e-14; end;
    end;
end;



%C  COMPUTE THE RADIUS OF CURVATURE IN THE PRIME VERTICAL FOR
%C EACH POINT

xnu=A./sqrt(1.0-(E*sin(lat)).^2);
xnu1=xnu(1:NN-1);
xnu2=xnu(2:NN);

%C*** COMPUTE THE AZIMUTHS.  A12 (A21) IS THE AZIMUTH AT POINT 1 (2)
%C OF THE NORMAL SECTION CONTAINING THE POINT 2 (1)

TPSI2=(1.-E*E)*tan(PHI2) + E*E*xnu1.*sin(PHI1)./(xnu2.*cos(PHI2));
PSI2=atan(TPSI2);

%C*** SOME FORM OF ANGLE DIFFERENCE COMPUTED HERE??

DPHI2=PHI2-PSI2;
DLAM=XLAM2-XLAM1;
CTA12=(cos(PHI1).*TPSI2 - sin(PHI1).*cos(DLAM))./sin(DLAM);
A12=atan((1.)./CTA12);
CTA21P=(sin(PSI2).*cos(DLAM) - cos(PSI2).*tan(PHI1))./sin(DLAM);
A21P=atan((1.)./CTA21P);

%C    GET THE QUADRANT RIGHT
DLAM2=(abs(DLAM)<pi).*DLAM + (DLAM>=pi).*(-2*pi+DLAM) + ...
    (DLAM<=-pi).*(2*pi+DLAM);
A12=A12+(A12<-pi)*2*pi-(A12>=pi)*2*pi;
A12=A12+pi*sign(-A12).*( sign(A12) ~= sign(DLAM2) );
A21P=A21P+(A21P<-pi)*2*pi-(A21P>=pi)*2*pi;
A21P=A21P+pi*sign(-A21P).*( sign(A21P) ~= sign(-DLAM2) );
%%A12*180/pi
%%A21P*180/pi


SSIG=sin(DLAM).*cos(PSI2)./sin(A12);
% At this point we are OK if the angle < 90...but otherwise
% we get the wrong branch of asin!
% This fudge will correct every case on a sphere, and *almost*
% every case on an ellipsoid (wrong hnadling will be when
% angle is almost exactly 90 degrees)
dd2=[cos(long).*cos(lat) sin(long).*cos(lat) sin(lat)];
dd2=sum((diff(dd2).*diff(dd2))')';
if ( any(abs(dd2-2) < 2*((B-A)/A))^2 ),
    disp('dist: Warning...point(s) too close to 90 degrees apart');
end;
bigbrnch=dd2>2;

SIG=asin(SSIG).*(bigbrnch==0) + (pi-asin(SSIG)).*bigbrnch;

SSIGC=-sin(DLAM).*cos(PHI1)./sin(A21P);
SIGC=asin(SSIGC);
A21 = A21P - DPHI2.*sin(A21P).*tan(SIG/2.0);

%C   COMPUTE RANGE

G2=EPS*(sin(PHI1)).^2;
G=sqrt(G2);
H2=EPS*(cos(PHI1).*cos(A12)).^2;
H=sqrt(H2);
TERM1=-SIG.*SIG.*H2.*(1.0-H2)/6.0;
TERM2=(SIG.^3).*G.*H.*(1.0-2.0*H2)/8.0;
TERM3=(SIG.^4).*(H2.*(4.0-7.0*H2)-3.0*G2.*(1.0-7.0*H2))/120.0;
TERM4=-(SIG.^5).*G.*H/48.0;

range=xnu1.*SIG.*(1.0+TERM1+TERM2+TERM3+TERM4);


if (geodes),

    %c now calculate the locations along the ray path. (for extra accuracy, could
    %c do it from start to halfway, then from end for the rest, switching from A12
    %c to A21...
    %c started to use Rudoe's formula, page 117 in Bomford...(1980, fourth edition)
    %c but then went to Clarke's best formula (pg 118)

    %RP I am doing this twice because this formula doesn't work when we go
    %past 90 degrees!
    Ngd1=round(Ngeodes/2);

    % First time...away from point 1
    if (Ngd1>1),
        wns=ones(1,Ngd1);
        CP1CA12 = (cos(PHI1).*cos(A12)).^2;
        R2PRM = -EPS.*CP1CA12;
        R3PRM = 3.0*EPS.*(1.0-R2PRM).*cos(PHI1).*sin(PHI1).*cos(A12);
        C1 = R2PRM.*(1.0+R2PRM)/6.0*wns;
        C2 = R3PRM.*(1.0+3.0*R2PRM)/24.0*wns;
        R2PRM=R2PRM*wns;
        R3PRM=R3PRM*wns;

        %c  now have to loop over positions
        RLRAT = (range./xnu1)*([0:Ngd1-1]/(Ngeodes-1));

        THETA = RLRAT.*(1 - (RLRAT.^2).*(C1 - C2.*RLRAT));
        C3 = 1.0 - (R2PRM.*(THETA.^2))/2.0 - (R3PRM.*(THETA.^3))/6.0;
        DSINPSI =(sin(PHI1)*wns).*cos(THETA) + ...
            ((cos(PHI1).*cos(A12))*wns).*sin(THETA);
        %try to identify the branch...got to other branch if range> 1/4 circle
        PSI = asin(DSINPSI);

        DCOSPSI = cos(PSI);
        DSINDLA = (sin(A12)*wns).*sin(THETA)./DCOSPSI;
        DTANPHI=(1.0+EPS)*(1.0 - (E^2)*C3.*(sin(PHI1)*wns)./DSINPSI).*tan(PSI);
        %C compute output latitude (phi) and long (xla) in radians
        %c I believe these are absolute, and don't need source coords added
        PHI = atan(DTANPHI);
        %  fix branch cut stuff -
        otherbrcnh= sign(DLAM2*wns) ~= sign([sign(DLAM2) diff(DSINDLA')'] );
        XLA = XLAM1*wns + asin(DSINDLA).*(otherbrcnh==0) + ...
            (pi-asin(DSINDLA)).*(otherbrcnh);
    else
        PHI=PHI1;
        XLA=XLAM1;
    end;

    % Now we do the same thing, but in the reverse direction from the receiver!
    if (Ngeodes-Ngd1>1),
        wns=ones(1,Ngeodes-Ngd1);
        CP2CA21 = (cos(PHI2).*cos(A21)).^2;
        R2PRM = -EPS.*CP2CA21;
        R3PRM = 3.0*EPS.*(1.0-R2PRM).*cos(PHI2).*sin(PHI2).*cos(A21);
        C1 = R2PRM.*(1.0+R2PRM)/6.0*wns;
        C2 = R3PRM.*(1.0+3.0*R2PRM)/24.0*wns;
        R2PRM=R2PRM*wns;
        R3PRM=R3PRM*wns;

        %c  now have to loop over positions
        RLRAT = (range./xnu2)*([0:Ngeodes-Ngd1-1]/(Ngeodes-1));

        THETA = RLRAT.*(1 - (RLRAT.^2).*(C1 - C2.*RLRAT));
        C3 = 1.0 - (R2PRM.*(THETA.^2))/2.0 - (R3PRM.*(THETA.^3))/6.0;
        DSINPSI =(sin(PHI2)*wns).*cos(THETA) + ...
            ((cos(PHI2).*cos(A21))*wns).*sin(THETA);
        %try to identify the branch...got to other branch if range> 1/4 circle
        PSI = asin(DSINPSI);

        DCOSPSI = cos(PSI);
        DSINDLA = (sin(A21)*wns).*sin(THETA)./DCOSPSI;
        DTANPHI=(1.0+EPS)*(1.0 - (E^2)*C3.*(sin(PHI2)*wns)./DSINPSI).*tan(PSI);
        %C compute output latitude (phi) and long (xla) in radians
        %c I believe these are absolute, and don't need source coords added
        PHI = [PHI fliplr(atan(DTANPHI))];
        % fix branch cut stuff
        otherbrcnh= sign(-DLAM2*wns) ~= sign( [sign(-DLAM2) diff(DSINDLA')'] );
        XLA = [XLA fliplr(XLAM2*wns + asin(DSINDLA).*(otherbrcnh==0) + ...
            (pi-asin(DSINDLA)).*(otherbrcnh))];
    else
        PHI = [PHI PHI2];
        XLA = [XLA XLAM2];
    end;

    %c convert to degrees
    A12 = PHI*180/pi;
    A21 = XLA*180/pi;
    range=range*([0:Ngeodes-1]/(Ngeodes-1));


else

    %C*** CONVERT TO DECIMAL DEGREES
    A12=A12*180/pi;
    A21=A21*180/pi;
    if (rowvec),
        range=range';
        A12=A12';
        A21=A21';
    end;
end;

end