function plot_weight

sensor_radius = 0.5;
separating    = 0.4999;%0:0.1:0.5;
Lin           = 4;


for i = 1:numel(separating)
    [d s]         = obtain_weights(separating(i),sensor_radius,Lin)
end




function [deep_weight,sup_weight] = obtain_weights(separating_radius,sensor_radius,Lin)

i =1;
deep_weight = zeros(1,(Lin+1)^2 - 1 );
Ls = [];
for L=1:Lin
    for M=[-L:L]
          
     
         deep_weight(i) = (sensor_radius^5 - separating_radius^5) .*  ( (separating_radius^(2*L-2)) ./ ...
             ((sensor_radius^(2*L+3)) - (separating_radius^(2*L+3)))  );
        
        % normalise and ensure that tends to unity as r > R 
        Ls = [Ls,L];
        i = i + 1;
    end
end
%sup_weight  = (1./deep_weight)./(1/deep_weight(end)) .* ((2*Lin)+3) ./ (2.*Ls+3);
deep_weight = deep_weight .* ((2*Ls)+3)/5;

sup_weight = (separating_radius.^(2*Lin - 2*Ls) .* (sensor_radius.^(2*Ls +3)) - separating_radius^(2*Lin+3)) ./...
    (sensor_radius^(2*Lin +3) - separating_radius^(2*Lin +3));

sup_weight = sup_weight .* ((2*Lin)+3) ./ (2.*Ls+3);

function [deep_weight,sup_weight] = obtain_weights_unnormalised(separating_radius,sensor_radius,Lin)

i =1;
deep_weight = zeros(1,(Lin+1)^2 - 1 );

for L=1:Lin
    for M=[-L:L]
        deep_weight(i) = (sensor_radius^3 - separating_radius^3) .*  ( (separating_radius^(2*L)) ./ ...
            ((sensor_radius^(2*L+3)) - (separating_radius^(2*L+3)))  ) .* (2*L+3) ./ 3 ;
        i = i + 1;
    end
end
sup_weight = 1./deep_weight;