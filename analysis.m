%% Initialization
clear;
basePath = 'raw\';
idx_list = setdiff(1:12, 3); nsample = numel(idx_list);

cnt_x1_nm = 7020; cnt_x2_nm = 48290;    
cnt_y1_nm = 7020; cnt_y2_nm = 48290;
pxsz = 21.6; % pixel size in nm for dSTORM
npxr = floor((cnt_y2_nm-cnt_y1_nm)/pxsz);
npxc = floor((cnt_x2_nm-cnt_x1_nm)/pxsz);

radius = 500;  % in nm
threshold = 10.^[2,.5,-inf];
ngroup = numel(threshold);
colors = [1, 0, 0;  0, 1, 0;  0, 0, 1];


%% Loop through all samples
[xy_tub, xy_tau, I_tub, I_tau] = deal(cell(nsample,1));
[dt, V, C, Ecc, Area, Area_1st, Area_cluster, Dist_NN1, Dist_NN10, Dist_1st, Density, group, group_cluster, coloc_0th, coloc_1st, coloc_local, coloc_cluster] = deal(cell(nsample,1));
[nsites,ncluster] = deal(zeros(nsample,1));
for n = 1:nsample
% for n = 1
    disp(n);

    % TAU DATA PROCESSING
    fname_tau = sprintf('%ssample%d_tau_registered.csv', basePath, idx_list(n));
    Pos_tau = readtable(fname_tau);
    % Imin_tau = quantile(Pos_tau{:, 4}, 0.9)
    Imin_tau = 500;
    Pos_tau = Pos_tau(Pos_tau{:, 4} >= Imin_tau, :);  % Intensity thresholding
    Pos_tau = Pos_tau(Pos_tau{:, 2} >= cnt_x1_nm & Pos_tau{:, 2} <= cnt_x2_nm & Pos_tau{:, 3} >= cnt_y1_nm & Pos_tau{:, 3} <= cnt_y2_nm, :);
    xy_tau{n}(:,1) = (Pos_tau{:, 2}-cnt_x1_nm)/pxsz;
    xy_tau{n}(:,2) = (Pos_tau{:, 3}-cnt_y1_nm)/pxsz;
    I_tau{n} = Pos_tau{:, 4};
    nsites(n) = length(xy_tau{n});

    % TUBULIN DATA PROCESSING
    fname_tub = sprintf('%ssample%d_tubulin.csv', basePath, idx_list(n)); 
    Pos_tub = readtable(fname_tub);
    Imin_tub = 1000;
    Pos_tub = Pos_tub(Pos_tub{:, 4} > Imin_tub, :);  % Intensity thresholding
    Pos_tub = Pos_tub(Pos_tub{:, 2} >= cnt_x1_nm & Pos_tub{:, 2} <= cnt_x2_nm & Pos_tub{:, 3} >= cnt_y1_nm & Pos_tub{:, 3} <= cnt_y2_nm, :);
    xy_tub{n}(:,1) = (Pos_tub{:, 2}-cnt_x1_nm)/pxsz;
    xy_tub{n}(:,2) = (Pos_tub{:, 3}-cnt_y1_nm)/pxsz;
    I_tub{n} = Pos_tub{:, 4};

    % Delaunay triangulation and Voronoi diagram calculation
    dt{n} = delaunayTriangulation(xy_tau{n}(:,1), xy_tau{n}(:,2));
    [V{n}, C{n}] = voronoiDiagram(dt{n});

    % CELL AREA, PROXIMITY, AND COLOCALIZATION RATE
    [Ecc{n}, Area{n}, Area_1st{n}, Dist_NN1{n}, Dist_NN10{n}, Dist_1st{n}, Density{n}, tubCount_0th, tubCount_1st, tubCount_local] = deal(nan(nsites(n), 1));
    for i = 1:nsites(n)
        if all(C{n}{i} ~= 1)  % Check if the cell is finite
            vertices = V{n}(C{n}{i}, :);
            centeredVertices = vertices - mean(vertices, 1);
            [~, eigenValues] = eig(cov(centeredVertices));
            lambda1 = max(diag(eigenValues));
            lambda2 = min(diag(eigenValues));
            Ecc{n}(i) = sqrt(1 - (lambda2 / lambda1));

            Area{n}(i) = polyarea(vertices(:,1), vertices(:,2)) * (pxsz/1e3)^2;  % in μm^2

            % PROXIMITY CALCULATION
            distances = sort(pdist2(xy_tau{n}(i, :), xy_tau{n}(1:end~=i, :)));
            Dist_NN1{n}(i) = distances(1) * (pxsz/1e3);  % in μm
            Dist_NN10{n}(i) = median(distances(1:10)) * (pxsz/1e3);  % in μm

            % Colocalization in the Voronoi cells
            IN = inpolygon(xy_tub{n}(:,1), xy_tub{n}(:,2), vertices(:,1), vertices(:,2));
            tubCount_0th(i) = sum(IN);

            % Colocalization around the seeds
            IN = hypot(xy_tub{n}(:,1) - xy_tau{n}(i, 1), xy_tub{n}(:,2) - xy_tau{n}(i, 2)) <= radius/pxsz;
            tubCount_local(i) = sum(IN);
        end
    end
    coloc_0th{n} = tubCount_0th ./ Area{n};
    coloc_local{n} = tubCount_local ./ (pi * (radius/1e3)^2);

    % First rank distance and density calculation (neighboring cells)        
    for i = 1:nsites(n)
        if all(C{n}{i} ~= 1)  % Check if the cell is finite
            j_list = setdiff(unique(dt{n}.ConnectivityList(any(dt{n}.ConnectivityList==i, 2),:)),i);
            Dist_1st{n}(i) = mean(hypot(xy_tau{n}(j_list,1) - xy_tau{n}(i, 1), xy_tau{n}(j_list,2) - xy_tau{n}(i, 2)));

            Area_1st{n}(i) = Area{n}(i) + sum(Area{n}(j_list));
            Density{n}(i) = (1 + numel(j_list))/Area_1st{n}(i);

            C_all = unique([C{n}{j_list}]);
            if all(C_all ~= 1)
                vertices_1st = V{n}(C_all, :);
                vertices_1st = vertices_1st(convhull(vertices_1st),:);
                IN = inpolygon(xy_tub{n}(:,1), xy_tub{n}(:,2), vertices_1st(:,1), vertices_1st(:,2));            
                tubCount_1st(i) = sum(IN);
            end
        end
    end
    coloc_1st{n} = tubCount_1st ./ Area_1st{n};

    % GROUPING OF AREA OR PROXIMITY
    group{n} = zeros(nsites(n), 1);
    sel1 = true(nsites(n), 1);
    for g = 1:ngroup
        % sel2 = Dist_NN1{n} <= quantile(Dist_NN1{n}, g/ngroup);
        sel2 = Density{n} >= threshold(g);
        % sel2 = Area{n} >= threshold(g);
        group{n}(sel1 & sel2) = g;
        sel1(sel2) = false;
    end
    for i = 1:nsites(n)
        if ~all(C{n}{i} ~= 1) || isnan(Density{n}(i))
            group{n}(i) = nan;
        end
    end

    % clustering
    group_cluster{n} = nan(nsites(n), 2);  % To store [group, cluster ID]
    clusterID = 0;  % Cluster counter
    for i = 1:nsites(n)
        if ~isnan(group{n}(i)) && isnan(group_cluster{n}(i, 2))  % Not yet assigned to a cluster (check for NaN)
            clusterID = clusterID + 1;  % Start a new cluster
            group_cluster{n}(i, :) = [group{n}(i), clusterID];
            
            % Initialize a stack to perform depth-first search (DFS)
            stack = i;
            
            while ~isempty(stack)
                idx = stack(end);  % Take the current index from the stack
                stack(end) = [];   % Pop from the stack
                
                % Find neighbors of the current site (idx)
                neighbors = setdiff(unique(dt{n}.ConnectivityList(any(dt{n}.ConnectivityList == idx, 2), :)), idx);
                
                for j = neighbors'
                    % Check if the neighbor belongs to the same group and is unassigned
                    if isnan(group_cluster{n}(j, 2)) && group{n}(j) == group{n}(i)
                        group_cluster{n}(j, :) = [group{n}(j), clusterID];  % Assign the same group and cluster ID
                        stack = [stack; j];  % Add the neighbor to the stack to continue DFS
                    end
                end
            end
        end
    end

    ncluster(n) = clusterID;
    [Area_cluster{n},coloc_cluster{n}] = deal(zeros(ncluster(n),1));
    for clusterID = 1:ncluster(n)
        sel = group_cluster{n}(:,2) == clusterID;
        Area_cluster{n}(clusterID) = sum(Area{n}(sel));
        coloc_cluster{n}(clusterID) = sum(coloc_0th{n}(sel).*Area{n}(sel))/Area_cluster{n}(clusterID);
    end
end
save('analysis');

%% check
maxfig(100); clf;

% % colocal
% n = 1;
% plot(xy_tub{n}(:,1), xy_tub{n}(:,2),'g.','markersize',5); hold all;
% plot(xy_tau{n}(:,1), xy_tau{n}(:,2),'r.','markersize',10);

% props
sel_n = setdiff(1:11,[1]);
for n = sel_n
    for pid = 1:4
        subplot(2,2,pid);
        if pid == 1
            edges = linspace(Imin_tau,5e4,100);
            % N = histcounts(I_tau{n}(:),edges);
            N = histcounts(I_tau{n}(:),edges,'norm','prob');
        elseif pid == 2
            edges = linspace(0,5e4,100);
            % N = histcounts(I_tub{n}(:),edges);
            N = histcounts(I_tub{n}(:),edges,'norm','prob');
        elseif pid == 3
            edges = linspace(-6,2,20);
            N = histcounts(log10(Area{n}(:)),edges);
        elseif pid == 4
            edges = linspace(-2,6,20);
            N = histcounts(log10(Density{n}(:)),edges);
        end
        xdat = edges(1:end-1) + diff(edges(1:2))/2;
        plot(xdat,N); hold all;
        % title(n);
        % pause;
    end
end

%% plot results: correlations between cell properties (area, distance, density)
sel_n = setdiff(1:11,[1]);
xdat_all_list = {log10(cat(1,Dist_NN1{sel_n})), log10(cat(1,Dist_NN10{sel_n})), log10(cat(1,Dist_1st{sel_n}))};
xlim_list = {[-3,1]-.5, [-3,1], [-3,1]+1.5};
xlabel_list = {'log(dist(NN1))', 'log(dist(NN10))', 'log(dist(1st))'};
ydat_all_list = {log10(cat(1,Area{sel_n})), log10(cat(1,Density{sel_n}))};
ylim_list = {[-6,2], [-2,6]};
ylabel_list = {'log(area)', 'log(1st-rank density)'};

maxfig(101); clf;
for pid = 1:3
    subplot(3,4,1+pid); 
    xdat_all = xdat_all_list{pid};
    edges = linspace(xlim_list{pid}(1),xlim_list{pid}(2),100);
    histogram(xdat_all,edges);
    xlim(xlim_list{pid});
    % ylim([0,1000]);
    xlabel(xlabel_list{pid});
end

% area or 1st-rank density
ncomp = 3;
options = statset('MaxIter',1000);
for pid1 = 1:2
    subplot(3,4,4*pid1+1);
    ydat_all = ydat_all_list{pid1};
    edges = linspace(ylim_list{pid1}(1),ylim_list{pid1}(2),100)';
    binsz = diff(edges(1:2));
    histogram(ydat_all,edges); hold all;
    % ylim([0,1000]);

    sel = ydat_all >= edges(1) & ydat_all <= edges(end);
    gmres = fitgmdist(ydat_all(sel),ncomp,'Options',options);
    xdat_fit = edges(1:end-1) + binsz/2;
    ydat_fit = diff(cdf(gmres,edges))*numel(ydat_all);
    plot(xdat_fit,ydat_fit,'r');
    for i = 1:ncomp
        ydat_fit = diff(normcdf(edges,gmres.mu(i),sqrt(gmres.Sigma(i))))*(numel(ydat_all)*gmres.ComponentProportion(i));
        plot(xdat_fit,ydat_fit);
    end
    xlim(ylim_list{pid1});
    xlabel(ylabel_list{pid1});

    for pid2 = 1:3
        subplot(3,4,4*pid1 + 1 + pid2);
        xdat_all = xdat_all_list{pid2};
        plot(xdat_all,ydat_all,'.','markersize',2);
        xlim(xlim_list{pid2});
        ylim(ylim_list{pid1});
        yxlabel(ylabel_list{pid1}, xlabel_list{pid2});
    end
end
subplot(341);
plot(ydat_all_list{1}, ydat_all_list{2},'.','markersize',2);
xlim(ylim_list{1});
ylim(ylim_list{2});
yxlabel(ylabel_list{2}, ylabel_list{1});
% saveas(gcf,'stats_90%-1000.jpg');
% saveas(gcf,'stats_1e4-all.jpg');
% saveas(gcf,'stats_1e3-all.jpg');
% saveas(gcf,'stats_1e3-all_del1.jpg');
% saveas(gcf,'stats_500-all_del1.jpg');
saveas(gcf,'stats_500-1000_del1.jpg');

%% plot results: correlations between cell properties and colocalization
xdat_all_list = {log10(cat(1,Area{sel_n})), log10(cat(1,Density{sel_n}))};
xlim_list = {[-6,2], [-2,6]};
xlabel_list = {'log(area)', 'log(1st-rank density)'};
ydat_all_list = {log10(cat(1,coloc_0th{sel_n})), log10(cat(1,coloc_1st{sel_n})), log10(cat(1,coloc_local{sel_n})), log10(cat(1,coloc_cluster{sel_n}))};
ylabel_list = {'log(coloc(0th))', 'log(coloc(1st))', 'log(coloc(local))', 'log(coloc(cluster))'};

maxfig(102); clf;
for pid = 1:4
    subplot(3,4,pid); 
    ydat_all = ydat_all_list{pid};
    edges = linspace(-2,6,100);
    histogram(ydat_all,edges);
    xlim([-2,6]);
    % ylim([0,1000]);
    xlabel(ylabel_list{pid});
    missed = sum(isinf(ydat_all)/numel(ydat_all));
    title(num2str(missed,'%.2f'));
end

% area or 1st-rank density
for pid1 = 1:2
    for pid2 = 1:3
        subplot(3,4,4*pid1 + pid2);
        xdat_all = xdat_all_list{pid1};
        ydat_all = ydat_all_list{pid2};
        plot(xdat_all,ydat_all,'.','markersize',2);
        xlim(xlim_list{pid1});
        ylim([-2,6]);
        yxlabel(ylabel_list{pid2}, xlabel_list{pid1});
    end
end

pid1 = 1; pid2 = 4;
subplot(3,4,4*pid1 + pid2);
xdat_all = log10(cat(1,Area_cluster{sel_n}));
ydat_all = ydat_all_list{pid2};
gdat_all = [];
for n = sel_n
    [~,ia] = unique(group_cluster{n}(:,2));
    gdat_all = [gdat_all; group_cluster{n}(ia(1:ncluster(n)),1)];
end
sel = gdat_all == 1;
plot(xdat_all(sel),ydat_all(sel),'.','markersize',10);
% plot(xdat_all,ydat_all,'.','markersize',10);
xlim(xlim_list{pid1});
ylim([-2,6]);
yxlabel(ylabel_list{pid2}, 'log(area(cluster))');
title('group#1');

subplot(3,4,12);
edges = linspace(xlim_list{pid1}(1),xlim_list{pid1}(2),100);
histogram(xdat_all(sel),edges);
xlim(xlim_list{pid1});
xlabel('log(area(cluster))');
title('group#1');

saveas(gcf,'colocalization_trend.jpg');

%%
maxfig(103); clf;
pid2_list = [1:6, 8+(1:6)];
xlim_list = {[-6,2], [-2,6]};
xlabel_list = {'log(area)', 'log(1st-rank density)'};
for pid1 = 1:2
    for n = 1:nsample
        for g = 1:ngroup
            sel = group{n} == g;
            if pid1 == 1
                xdat = log10(Area{n}(sel));
            else
                xdat = log10(Density{n}(sel));
            end
            ydat = log10(coloc_local{n}(sel));
            subplot(4,8,16*(pid1-1)+pid2_list(n));
            plot(xdat,ydat,'.','markersize',5,'color',colors(g,:)); hold all;
            subplot(2,4,4*(pid1-1)+4);
            plot(xdat,ydat,'.','markersize',5,'color',colors(g,:)); hold all;
        end
        subplot(4,8,16*(pid1-1)+pid2_list(n));
        xlim(xlim_list{pid1}); ylim([-2,6]);
        yxlabel('log(coloc(local))',xlabel_list{pid1});
    end
    subplot(2,4,4*(pid1-1)+4);
    xlim(xlim_list{pid1}); ylim([-2,6]);
    yxlabel('log(coloc(local))',xlabel_list{pid1});
end
saveas(gcf,'colocalization_by image and group.jpg');


maxfig(104); clf;
gdat_all = cat(1,group{sel_n});
for pid = 1:3
    subplot(2,2,pid);
    ydat_all = ydat_all_list{pid};
    sel = ~isinf(ydat_all);
    h = violinplot(ydat_all(sel),gdat_all(sel));
    ylim([-2,6]);
    
    clear missed;
    for i = 1:ngroup
        h(i).ViolinColor = {colors(i, :)};  % Wrap color in a cell array
        missed(i) = sum(gdat_all==i & isinf(ydat_all))/numel(ydat_all);
    end
    ylabel(ylabel_list{pid});
    title(num2str(missed,'%.2f '));
end

pid = 4;
subplot(2,2,pid);
gdat_all = [];
for n = sel_n
    [~,ia] = unique(group_cluster{n}(:,2));
    gdat_all = [gdat_all; group_cluster{n}(ia(1:ncluster(n)),1)];
end
ydat_all = log10(cat(1,coloc_cluster{sel_n}));
sel = ~isinf(ydat_all);
h = violinplot(ydat_all(sel),gdat_all(sel));
ylim([-2,6]);

clear missed;
for i = 1:ngroup
    h(i).ViolinColor = {colors(i, :)};  % Wrap color in a cell array
    missed(i) = sum(gdat_all==i & isinf(ydat_all))/numel(ydat_all);
end
ylabel(ylabel_list{pid});
title(num2str(missed,'%.2f '));
saveas(gcf,'colocalization_violin.jpg');


%% voronoi cell shape
xdat_all_list = {log10(cat(1,Area{sel_n})), log10(cat(1,Density{sel_n}))};
xlim_list = {[-6,2], [-2,6]};
xlabel_list = {'log(area)', 'log(1st-rank density)'};

maxfig(105); clf;
for pid = 1:2
    subplot(2,2,pid);
    xdat = xdat_all_list{pid};
    ydat = cat(1,Ecc{sel_n});
    plot(xdat,ydat,'.','markersize',5);
    xlim(xlim_list{pid});
    yxlabel('Eccentricity',xlabel_list{pid});
end

%% Voronoi results
Imin = [1e4,5e3];
Imax = [2e4,5e4];
Imax2 = Imin + (Imax-Imin)*5;
img_rgb_all = zeros(npxr,npxc,3,nsample, 'uint8');

for n = 1:nsample
% for n = 4
    disp(n);
    maxfig(n); clf;
    
    ax(1) = subplot2(1,2,1);
    % Load and preprocess the raw images
    fname_img = sprintf('%ssample%d_tau.tif', basePath, idx_list(n));
    img_tau = imread(fname_img);
    fname_img = sprintf('%ssample%d_tubulin.tif', basePath, idx_list(n));
    img_tub = imread(fname_img);

    x0 = (size(img_tau, 2) - npxc) / 2;
    y0 = (size(img_tau, 1) - npxr) / 2;

    img_r = imcrop(img_tau, [x0 y0 npxc-1 npxr-1]);
    img_g = imcrop(img_tub, [x0 y0 npxc-1 npxr-1]);
    img_r = rescale(img_r,0,255,"InputMin",Imin(1),"InputMax",Imax(1));
    img_g = rescale(img_g,0,255,"InputMin",Imin(2),"InputMax",Imax(2));     
    img_rgb = cat(3,img_r,img_g,zeros(npxr,npxc,'uint8'));
    img_rgb_all(:,:,:,n) = img_rgb;
    if n == 4
        img_rgb_selected = img_rgb;
        img_r = imcrop(img_tau, [x0 y0 npxc-1 npxr-1]);
        img_g = imcrop(img_tub, [x0 y0 npxc-1 npxr-1]);
        img_r = rescale(img_r,0,255,"InputMin",Imin(1),"InputMax",Imax2(1));
        img_g = rescale(img_g,0,255,"InputMin",Imin(2),"InputMax",Imax2(2));     
        img_rgb = cat(3,img_r,img_g,zeros(npxr,npxc,'uint8'));
        img_rgb_selected_lowcontrast = img_rgb;
    end

    imshow2(img_rgb); hold on;
    xlim([.5,npxc+.5]); ylim([.5,npxr+.5]);
    axis square on;
    title('tau (r) + tubulin (g)');
    

    ax(2) = subplot2(1,2,2); axis square ij; hold on;
    % img_rgb = cat(3,img_r,zeros(npxr,npxc,2,'uint8'));
    % imshow2(img_rgb);
    % imshow2(img_r,[0,255]);

    % % Overlay the Voronoi diagram
    % [vx, vy] = voronoi(dt{n});
    % plot(vx, vy, 'w', 'LineWidth', .5); % Plotting Voronoi lines
    
    % draw borders only
    idx = find(isnan(group{n}));
    vertices = [];
    for i = 1:numel(idx)
        tmp = V{n}(C{n}{idx(i)},:);
        vertices = [vertices; tmp; tmp(1,:); nan, nan];
    end
    plot(vertices(:,1), vertices(:,2), 'color',rgb('gray'), 'linewidth',.5);

    for g = ngroup:-1:1
        idx = find(group{n} == g);
        vertices = [];
        for i = 1:numel(idx)
            tmp = V{n}(C{n}{idx(i)},:);
            vertices = [vertices; tmp; tmp(1,:); nan, nan];
        end
        plot(vertices(:,1), vertices(:,2), 'color',colors(g,:), 'linewidth',.5);
    end

    % Scatter plot all dots colored by group
    for g = 1:ngroup
        scatter(dt{n}.Points(group{n} == g, 1), dt{n}.Points(group{n} == g, 2), 10, colors(g,:), 'filled', 'MarkerEdgeColor', 'none');
    end
    
    % % Fill each Voronoi cell with corresponding color and set transparency
    % for i = 1:nsites(n)
    %     if ~isnan(group{n}(i))
    %         fill(V{n}(C{n}{i}, 1), V{n}(C{n}{i}, 2), colors(group{n}(i), :), 'FaceAlpha', 0.1, 'edgecolor',rgb('gray'));
    %     end
    % end

    % % Fill each Voronoi cell with corresponding group color and set transparency
    % for i = 1:ncluster(n)
    %     sel = group_cluster{n}(:,2) == i;
    %     vertices_cluster = V{n}(unique(cat(2,C{n}{sel})), :);
    %     vertices_cluster = vertices_cluster(convhull(vertices_cluster),:);
    %     group_tmp = group_cluster{n}(sel,1);
    %     fill(vertices_cluster(:,1), vertices_cluster(:,2), colors(group_tmp(1),:), 'FaceAlpha', 0.5);
    % 
    %     % % Calculate the centroid of the Voronoi region
    %     % x_centroid = mean(V{n}(C{n}{i}, 1));
    %     % y_centroid = mean(V{n}(C{n}{i}, 2));
    %     % 
    %     % % Label the centroid with the cluster ID
    %     % text(x_centroid, y_centroid, num2str(group_cluster{n}(i, 2)), 'Color', 'w', 'FontSize', 10, 'HorizontalAlignment', 'center');
    % end
    
    xlim([.5,npxc+.5]); ylim([.5,npxr+.5]);
    title('Voronoi diagram with tau overlay');
    linkaxes(ax(1:2),'xy');
    set(ax(:), 'DataAspectRatio', [1 1 1]); % Fix aspect ratio for ax1
    
    % Save the figure
    finalFileName = sprintf('voronoi results\\sample%d.png', idx_list(n));
    exportgraphics(gcf, finalFileName, 'Resolution', 300);
    saveas(gcf,regexprep(finalFileName,'png','fig'));
    % pause;
end
save('analysis');