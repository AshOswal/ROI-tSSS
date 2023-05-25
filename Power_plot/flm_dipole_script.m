function Flmr_sum = flm_dipole_script(S)
%rq_depths = [0.04;0.07];
rq_depths = [0.04];
%radii     = [0.10;0.12;0.15];
radii     = 0.08;
roi_rad   = 0.02:0.02:0.08; 
Lin = 11;

load dipoles.mat Q RQ;
% Q is primary current and RQ is the position of each dipole

for i = 1:length(rq_depths)
    for j = 1:size(RQ,2)
        RQscaled(:,j) = rq_depths(i)*RQ(:,j)/norm(RQ(:,j));
    end
       
    if  1
        [x,y,z] = sphere(100);
        figure;plot3(RQscaled(1,:),RQscaled(2,:),RQscaled(3,:),'r.','MarkerSize',9);hold on;
        quiver3(RQscaled(1,:),RQscaled(2,:),RQscaled(3,:),Q(1,:),Q(2,:),Q(3,:),'k');hold on;
        s = surf(rq_depths(i)*x,rq_depths(i)*y,rq_depths(i)*z);
        s.EdgeColor = 'none';
        s.FaceAlpha = 0.2;
        s.FaceColor = [0 0 1];
        xlabel('X');ylabel('Y');zlabel('Z');axis equal
    end
   
    for r = 1:numel(roi_rad)
        if ~isequal(roi_rad(r),rq_depths) && strcmp(S.ROI,'sphere')
            sroi = surf(roi_rad(r)*x,roi_rad(r)*y,roi_rad(r)*z);
            sroi.EdgeColor = 'none';
            sroi.FaceAlpha = 0.1;
            sroi.FaceColor = [0.2+0.1*r 0.2+0.1*r 0.2+0.1*r];hold on
            xlim([-0.08 0.08]);
            ylim([[-0.08 0.08]]);
            zlim([[-0.08 0.08]]);
        elseif strcmp(S.ROI,'cube')
            % cube side is determined through application of pythagoras'
            cube_side = 2*roi_rad(r)/sqrt(3);
            plotcube([cube_side cube_side cube_side],[ -cube_side/2  -cube_side/2  -cube_side/2],.1,[0.2+0.1*r 0.2+0.1*r 0.2+0.1*r]);
            hold on;% plot surface of cube
            xlim([-0.08 0.08]);
            ylim([[-0.08 0.08]]);
            zlim([[-0.08 0.08]]);
        end
        for j = 1:length(radii)
            %disp([i j]);
            clear Flm;
            for k = 1:size(RQscaled,2)
                Flm(:,k) = flm_dipole(radii(j),RQscaled(:,k),Q(:,k),Lin)';
            end
            
            % obtain ROI weights
            if strcmp(S.ROI,'sphere')
                G = obtain_weights(roi_rad(r),radii(j),Lin)';
                if isequal(roi_rad(r),radii(j))
                    G = 1;
                end
            elseif strcmp(S.ROI,'cube')
                % side of cube
                
                G = GsGd(Lin,Lin,[cube_side radii(j)],'inner');
                G = diag(G);
            end
            
            Flmr = real(Flm);
            Flmr_sum{i}(j,r,:) = G.*(sum(Flmr'))';
        end
    end
end


function [deep_weight,sup_weight] = obtain_weights(separating_radius,sensor_radius,Lin)


ii =1;
deep_weight = zeros(1,(Lin+1)^2 - 1 );
Ls = [];
for L=1:Lin
    for M=-L:L

        deep_weight(ii) = (sensor_radius^5 - separating_radius^5) .*  ( (separating_radius^(2*L-2)) ./ ...
            ((sensor_radius^(2*L+3)) - (separating_radius^(2*L+3)))  );

        % normalise and ensure that tends to unity as r > R
        Ls = [Ls,L];
        ii = ii + 1;
    end
end
%sup_weight  = (1./deep_weight)./(1/deep_weight(end)) .* ((2*Lin)+3) ./ (2.*Ls+3);
deep_weight = deep_weight .* ((2*Ls)+3)/5;

sup_weight = (separating_radius.^(2*Lin - 2*Ls) .* (sensor_radius.^(2*Ls +3)) - separating_radius^(2*Lin+3)) ./...
    (sensor_radius^(2*Lin +3) - separating_radius^(2*Lin +3));

sup_weight = sup_weight .* ((2*Lin)+3) ./ (2.*Ls+3);


%{
ii =1;
deep_weight = zeros(1,(Lin+1)^2 - 1 );
Ls = [];
for L=1:Lin
    for M=[-L:L]
        %         deep_weight(ii) = (sensor_radius^3 - separating_radius^3) .*  ( (separating_radius^(2*L)) ./ ...
        %             ((sensor_radius^(2*L+3)) - (separating_radius^(2*L+3)))  ) .* (2*L+3) ./ 3 ;
        deep_weight(ii) = (sensor_radius^3 - separating_radius^3) .*  ( (separating_radius^(2*L)) ./ ...
            ((sensor_radius^(2*L+3)) - (separating_radius^(2*L+3)))  );
        
        Ls = [Ls,L];

        %deep_weight(ii) = deep_weight(ii).* (2*L+3)./5;
        ii = ii + 1;
    end
end
deep_weight = deep_weight./deep_weight(1).* ((2*Ls)+3)/5;
sup_weight = 1./deep_weight;
end

end
%}