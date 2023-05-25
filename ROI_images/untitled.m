M = gifti(fullfile(spm('dir'), 'canonical', 'cortex_8196.surf.gii'));
ir = find(M.vertices(:,1)<0);
M.vertices(ir,:) = [];

A = spm_mesh_adjacency(M.vertices)
ii=[];
for i= 1:size(M.faces,1)
    if any(intersect(M.faces(i,:),ir))
        ii = [ii i];
    end
end
M.faces(ii,:) = [];

figure;spm_mesh_render('Disp', A);