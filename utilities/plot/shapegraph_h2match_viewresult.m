% shapegraph_h2match_viewresult.m: function that generates various plots after performing shape graph
% registration using second order elastic Sobolev metrics (H2 metrics).
%
% Input:
%   optPath        : structure containing optimal deformation path of transformed source
%   transfSource   : structure containing transformed source
%   transfTarget   : structure containing transformed target
%   updatedSource  : structure containing source with reordered vertices and component curves after applying Tarjan's algorithm
%   summary        : structure containing costs, energies and information from the optimization process
%   num_time_pts   : number of (evenly spaced) time points at which to visualize geodesic path (optional, default is 5)
%   save_fig       : boolean parameter indicating whether to save figures of the geodesic path (optional, default is false)
%   file_name      : desired file name for saving figures [string]

% Ouputs: 
%   Figure(1)      : geometric plot of overlayed source (blue), transformed source (black) and target (red)
%   Figure(2)      : plot of edge weights for transformed source (black) and target (red)
%   Figure(3)      : geometric plot of overlayed source (blue), transformed source (black) and target (red), colored 
%                    with edge weights by transparency (edges with higher weights are opaque, edges with smaller
%                    weights are transparent), with geodesic path of transformed source overlayed on the plot
%   Figure(4)      : Time series plot of geodesic path of transformed source


function shapegraph_h2match_viewresult(optPath, transfSource, transfTarget, updatedSource, summary, varargin) 


% Dimension
[~,d] = size(transfTarget.x);

% Set default parameters
num_time_pts = 5;
save_fig = false;
file_name = 'matching';

% Use user-specified parameters
ii = 1;
while ii <= length(varargin)
    if (isa(varargin{ii}, 'char'))
        switch (lower(varargin{ii}))
            case 'num_time_pts'
                ii = ii + 1;
                num_time_pts = varargin{ii};
            case 'save_fig'
                ii = ii + 1;
                save_fig = varargin{ii};
            case 'file_name'
                ii = ii + 1;
                file_name = varargin{ii};
        end
    end
    ii = ii + 1;
end

% Extract objective function parameters and spline data
objfun = summary.parameters.objfun;
splineData = objfun.splineData;


%% Display H2 energy of path
energy_path = summary.costs.energy_path;
disp('Geodesic distance between source and target is:');
disp((energy_path)^(1/2));


%% Plot (a): Geometric plots of source, transformed source and target
figure(1)
hold on

% target
plot_shape(transfTarget,'colorStyle','r.-','lineWidth',2)

for k = 1:length(updatedSource.connComp)
    
    % source
    for l = 1:size(updatedSource.connComp{k}.G,1)   
        
        if d == 2
            plot([updatedSource.connComp{k}.x(updatedSource.connComp{k}.G(l,1),1) updatedSource.connComp{k}.x(updatedSource.connComp{k}.G(l,2),1)],...
                 [updatedSource.connComp{k}.x(updatedSource.connComp{k}.G(l,1),2) updatedSource.connComp{k}.x(updatedSource.connComp{k}.G(l,2),2)],...
                 'b.-','LineWidth', 2)
             
        elseif d == 3
            plot3([updatedSource.connComp{k}.x(updatedSource.connComp{k}.G(l,1),1) updatedSource.connComp{k}.x(updatedSource.connComp{k}.G(l,2),1)],...
                  [updatedSource.connComp{k}.x(updatedSource.connComp{k}.G(l,1),2) updatedSource.connComp{k}.x(updatedSource.connComp{k}.G(l,2),2)],...
                  [updatedSource.connComp{k}.x(updatedSource.connComp{k}.G(l,1),3) updatedSource.connComp{k}.x(updatedSource.connComp{k}.G(l,2),3)],...
                  'b.-','LineWidth', 2)
            grid on
            view(d)
            
        end
        
    end
    
    % transformed source
    for l = 1:size(transfSource.connComp{k}.G,1)
        
        if d == 2
            plot([transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,1),1) transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,2),1)],...
                 [transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,1),2) transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,2),2)],...
                 'k.-','LineWidth', 2)    
             
        elseif d == 3
            plot3([transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,1),1) transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,2),1)],...
                  [transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,1),2) transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,2),2)],...
                  [transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,1),3) transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,2),3)],...
                  'k.-','LineWidth', 2) 
            grid on
            view(d)
             
        end
        
    end    
    
end
axis equal off
hold off
set(gcf, 'Color', 'w')


%% Geometric plots with weights
switch lower(objfun.edge_weight)
    
    case {'both','source','target','fixed_weights'}

        % Plot (b): Plot of edge weights for transformed source and target
        figure(2)
        hold on
        num_edges = 0;
        
        % transformed source edge weights
        for k = 1:length(updatedSource.connComp)  
            plot(num_edges+1:num_edges+length(transfSource.connComp{k}.rho), 1+transfSource.connComp{k}.rho,'k.-','LineWidth', 1.2)  
            num_edges = num_edges + length(transfSource.connComp{k}.rho);
        end
        
        % target edge weights
        plot(1:length(transfTarget.rho), 1+transfTarget.rho,'r.-','LineWidth', 1.2) 
        
        % Adding legend to plot
        h = zeros(2, 1);
        h(1) = plot(NaN,NaN,'k.-'); 
        h(2) = plot(NaN,NaN,'r.-'); 
        legend(h,'transformed source','target')
        xlabel('edge')
        ylabel('edge weight') 
        hold off
        set(gcf, 'Color', 'w')
        

        % Plot (c): Geometric plots colored with edge weights, with correspondence lines and geodesic
        col_num = 40;
        minrho_source = min([0 min(1+updatedSource.rho)-0.05]);
        maxrho_source = max([1 max(1+updatedSource.rho)+0.05]);
        minrho_target = min([0 min(1+transfTarget.rho)-0.05]);
        maxrho_target = max([1 max(1+transfTarget.rho)+0.05]);
        minrho_endcurve = min([0 min(1+transfSource.rho)-0.05]);
        maxrho_endcurve = max([1 max(1+transfSource.rho)+0.05]);
        
        transparency = linspace(0.1,1,col_num)'; 
        
        % Deletion of creation of mass parameter (for plotting of geodesics)
        delete_mass = zeros(length(updatedSource.connComp),1);
        
        for k = 1:length(updatedSource.connComp)
            if 1+min(updatedSource.connComp{k}.rho) > 1+min(transfSource.connComp{k}.rho)
                delete_mass(k) = true;
            
            else
                delete_mass(k) = false;
            end
        end
        
            
        figure(3)
        hold on
            
        % target
        for l = 1:size(transfTarget.G,1)
            col_ind = round(col_num*(1+transfTarget.rho(l)-minrho_target)/(maxrho_target-minrho_target)); % target edge weight color
                    
            if d == 2
                line([transfTarget.x(transfTarget.G(l,1),1) transfTarget.x(transfTarget.G(l,2),1)],...
                     [transfTarget.x(transfTarget.G(l,1),2) transfTarget.x(transfTarget.G(l,2),2)],...
                     'color', [1, 0, 0, transparency(max(min(col_num, col_ind+1),1))], 'LineWidth', 2)
                     
            elseif d == 3
                line([transfTarget.x(transfTarget.G(l,1),1) transfTarget.x(transfTarget.G(l,2),1)],...
                     [transfTarget.x(transfTarget.G(l,1),2) transfTarget.x(transfTarget.G(l,2),2)],...
                     [transfTarget.x(transfTarget.G(l,1),3) transfTarget.x(transfTarget.G(l,2),3)],...
                     'color', [1, 0, 0, transparency(max(min(col_num, col_ind+1),1))], 'LineWidth', 2)
                 
                 grid on
                 view(d)
                     
            end
                    
        end
                
        % source
        for l = 1:size(updatedSource.G,1)
            col_ind = round(col_num*(1+updatedSource.rho(l)-minrho_source)/(maxrho_source-minrho_source)); % source edge weight color
                    
            if d == 2
                plot([updatedSource.x(updatedSource.G(l,1),1) updatedSource.x(updatedSource.G(l,2),1)],...
                     [updatedSource.x(updatedSource.G(l,1),2) updatedSource.x(updatedSource.G(l,2),2)],...
                     'color', [0, 0, 1, transparency(max(min(col_num, col_ind+1),1))],'LineWidth', 2)    
                    
            elseif d == 3
                plot3([updatedSource.x(updatedSource.G(l,1),1) updatedSource.x(updatedSource.G(l,2),1)],...
                      [updatedSource.x(updatedSource.G(l,1),2) updatedSource.x(updatedSource.G(l,2),2)],...
                      [updatedSource.x(updatedSource.G(l,1),3) updatedSource.x(updatedSource.G(l,2),3)],...
                      'color', [0, 0, 1, transparency(max(min(col_num, col_ind+1),1))],'LineWidth', 2)    
                 
            grid on
            view(d)
                    
            end
                    
        end
                
        % transformed source
        for k = 1:length(transfSource.connComp)
                    
            for l = 1:size(transfSource.connComp{k}.G,1)
                col_ind = round(col_num*(1+transfSource.connComp{k}.rho(l)-minrho_endcurve)/(maxrho_endcurve-minrho_endcurve)); % transformed source edge weight color        
                        
                if d == 2
                    plot([transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,1),1) transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,2),1)],...
                         [transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,1),2) transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,2),2)],...
                         'color', [0, 0, 0, transparency(max(min(col_num,col_ind+1),1))], 'LineWidth', 1.5)
                             
                elseif d == 3
                	plot3([transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,1),1) transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,2),1)],...
                          [transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,1),2) transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,2),2)],...
                          [transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,1),3) transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,2),3)],...
                          'color', [0, 0, 0, transparency(max(min(col_num,col_ind+1),1))],'LineWidth', 1.5)
                     
                grid on
                view(d)
                        
                end
                        
             end
                    
         end

         caxis([0 1])
                
         % Plot correspondence lines between source and transformed source
         for k = 1:length(updatedSource.connComp)
    
            for l = 1:size(updatedSource.connComp{k}.x,1)
                        
                if d == 2
                    plot([updatedSource.connComp{k}.x(l,1) transfSource.connComp{k}.x(l,1)],...
                         [updatedSource.connComp{k}.x(l,2) transfSource.connComp{k}.x(l,2)],'w--', 'LineWidth', 0.1) 
                         
                elseif d == 3
                    plot3([updatedSource.connComp{k}.x(l,1) transfSource.connComp{k}.x(l,1)],...
                          [updatedSource.connComp{k}.x(l,2) transfSource.connComp{k}.x(l,2)],...
                          [updatedSource.connComp{k}.x(l,3) transfSource.connComp{k}.x(l,3)],'w--', 'LineWidth', 0.1)
                     
                grid on
                view(d)
                          
                end
                        
            end
    
         end
                
        % Overlay geodesic
        Nt = splineData{1}.Nt;
        geodesic_transparency = linspace(transparency(1), transparency(end), Nt);

        for k = 1:length(updatedSource.connComp)

            cPath = optPath{k};
            Nj = splineData{k}.N;

            for t = 1:Nt-1
                            
            curve_t.x = splineData{k}.quadData.B_endCurveS * cPath((t-1)*Nj+1:t*Nj,:);
            curve_t.G = updatedSource.connComp{k}.G;
            
            if delete_mass(k)
                transparency_t = linspace(geodesic_transparency(end-t), geodesic_transparency(end), col_num);
            else
                transparency_t = linspace(geodesic_transparency(t), geodesic_transparency(end), col_num);
            end
                            
                for l = 1:size(curve_t.G,1)
                    
                    if delete_mass(k)  % use transformed source edge weights for progressive shading
                        col_ind = round(col_num*(1+transfSource.connComp{k}.rho(l)-minrho_endcurve)/(maxrho_endcurve-minrho_endcurve));
                    else            % use source edge weights for progressive shading
                        col_ind = round(col_num*(1+updatedSource.connComp{k}.rho(l)-minrho_source)/(maxrho_source-minrho_source));
                    end

                    if d == 2       
                        plot([curve_t.x(curve_t.G(l,1),1) curve_t.x(curve_t.G(l,2),1)],...
                             [curve_t.x(curve_t.G(l,1),2) curve_t.x(curve_t.G(l,2),2)],...
                             'color', [0, 0, 0, transparency_t(max(min(col_num,col_ind+1),1))], 'LineWidth', 0.6, 'LineStyle', '--') 

                    elseif d == 3    
                        plot3([curve_t.x(curve_t.G(l,1),1) curve_t.x(curve_t.G(l,2),1)],...
                              [curve_t.x(curve_t.G(l,1),2) curve_t.x(curve_t.G(l,2),2)],...
                              [curve_t.x(curve_t.G(l,1),3) curve_t.x(curve_t.G(l,2),3)],...
                              'color', [0, 0, 0, transparency_t(max(min(col_num,col_ind+1),1))], 'LineWidth', 0.6, 'LineStyle', '--') 
                         
                    grid on
                    view(d)
                                
                    end

                end

            end
            
        end
                
        caxis([0 1])
        axis equal off
        hold off
        set(gcf, 'Color', 'w')
        
end


%% Plot (d): Time Series Geodesic Plot

% Extract parameters for time series plot
Nt = splineData{1}.Nt;
time_pts = linspace(0,1,num_time_pts);
t_ind = min(Nt, max(1,ceil(time_pts*Nt)));
geodesic_transparency = geodesic_transparency(t_ind);

minplot=1.1*min([min(transfTarget.x,[],1); min(updatedSource.x,[],1)],[],1);
maxplot=1.1*max([max(transfTarget.x,[],1); max(updatedSource.x,[],1)],[],1);

plot_ind = 0;
figure(4)
hold on

% Plot geodesic
for tt = 1:length(time_pts)
    plot_ind = plot_ind+1;
    t = time_pts(tt);
    figure(4)
    clf
    hold on
    
    for k = 1:length(updatedSource.connComp)
        
        cPath = optPath{k};
        curve_t.x = splineData{k}.quadData.B_endCurveS * evalPath(cPath, t, splineData{k}); % transform control points to vertices
        curve_t.G = updatedSource.connComp{k}.G;
        
        if delete_mass(k)
            transparency_t = linspace(geodesic_transparency(max(min(end-tt,end),1)), geodesic_transparency(end), col_num);
        else
            transparency_t = linspace(geodesic_transparency(tt), geodesic_transparency(end), col_num);
        end
        
        for l = 1:size(curve_t.G,1)
            
            if delete_mass(k)  % use transformed source edge weights for progressive shading
                col_ind = round(col_num*(1+transfSource.connComp{k}.rho(l)-minrho_endcurve)/(maxrho_endcurve-minrho_endcurve));
            else            % use source edge weights for progressive shading
                col_ind = round(col_num*(1+updatedSource.connComp{k}.rho(l)-minrho_source)/(maxrho_source-minrho_source));
            end
            
            if d == 2
                plot([curve_t.x(curve_t.G(l,1),1) curve_t.x(curve_t.G(l,2),1)],...
                     [curve_t.x(curve_t.G(l,1),2) curve_t.x(curve_t.G(l,2),2)],...
                     'color', [0, 0, 1, transparency_t(max(min(col_num, col_ind+1), 1))], 'LineWidth', 2, 'LineStyle', '-')
                
            elseif d == 3
                plot3([curve_t.x(curve_t.G(l,1),1) curve_t.x(curve_t.G(l,2),1)],...
                      [curve_t.x(curve_t.G(l,1),2) curve_t.x(curve_t.G(l,2),2)],...
                      [curve_t.x(curve_t.G(l,1),3) curve_t.x(curve_t.G(l,2),3)],...
                      'color', [0, 0, 1, transparency_t(max(min(col_num,col_ind+1),1))], 'LineWidth', 2, 'LineStyle', '-')
                
                grid on
                view(d)
                
            end
            
        end
        
    end
    caxis([0 1])
    axis([minplot(1) maxplot(1) minplot(2) maxplot(2)]);
    axis equal off
    hold off
    set(gcf, 'Color', 'w')
    pause(0.5)
    
    if save_fig == true
        saveas(figure(4),['./Figures/' file_name '_t' num2str(plot_ind) '.png'],'png')  
    end
    
end

% Overlay target and transformed source
figure(4)
clf
hold on

% transformed source
for k = 1:length(transfSource.connComp)

    for l = 1:size(transfSource.connComp{k}.G,1)
        col_ind = round(col_num*(1+transfSource.connComp{k}.rho(l)-minrho_endcurve)/(maxrho_endcurve-minrho_endcurve)); % transformed source edge weight color        
        
        if d == 2
            plot([transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,1),1) transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,2),1)],...
                 [transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,1),2) transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,2),2)],...
                 'color', [0, 0, 1, transparency_t(max(min(col_num,col_ind+1),1))], 'LineWidth', 1.5, 'LineStyle', '--')

        elseif d == 3
            plot3([transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,1),1) transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,2),1)],...
                  [transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,1),2) transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,2),2)],...
                  [transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,1),3) transfSource.connComp{k}.x(transfSource.connComp{k}.G(l,2),3)],...
                  'color', [0, 0, 1, transparency_t(max(min(col_num,col_ind+1),1))],'LineWidth', 1.5, 'LineStyle', '--')

        grid on
        view(d)

        end

     end

 end


% target
for l = 1:size(transfTarget.G,1)
    col_ind = round(col_num*(1+transfTarget.rho(l)-minrho_target)/(maxrho_target-minrho_target)); % target edge weight color
                    
        if d == 2
            plot([transfTarget.x(transfTarget.G(l,1),1) transfTarget.x(transfTarget.G(l,2),1)],...
                 [transfTarget.x(transfTarget.G(l,1),2) transfTarget.x(transfTarget.G(l,2),2)],...
                 'color', [1, 0, 0, transparency(max(min(col_num,col_ind+1),1))], 'LineWidth', 2)
         
        elseif d == 3
            plot3([transfTarget.x(transfTarget.G(l,1),1) transfTarget.x(transfTarget.G(l,2),1)],...
                  [transfTarget.x(transfTarget.G(l,1),2) transfTarget.x(transfTarget.G(l,2),2)],...
                  [transfTarget.x(transfTarget.G(l,1),3) transfTarget.x(transfTarget.G(l,2),3)],...
                  'color', [1, 0, 0, transparency(max(min(col_num,col_ind+1),1))], 'LineWidth', 2)
                 
        grid on
        view(d)
         
        end
                    
end

caxis([0 1])
axis([minplot(1) maxplot(1) minplot(2) maxplot(2)]);
axis equal off
hold off
set(gcf, 'Color', 'w')

if save_fig == true
    saveas(figure(4),['./Figures/' file_name '_t' num2str(plot_ind) '.png'],'png')  
end

end

