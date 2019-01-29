function [mesh,data] = get_mesh(i,data)
% Decide what brain we're actually using, return it

try   mesh = i.g;
      fprintf('Using user provided mesh\n');
      
      if ischar(mesh) || isnumeric(mesh)
          [mesh,data] = aplot.convert_mesh(mesh,data);
      end
      
catch mesh = read_nv();
      fprintf('(Using template brain mesh)\n');
end

% if i.inflate
%     try fprintf('Trying to inflate mesh\n');
%         dmesh.vertices = mesh.vertices;
%         dmesh.faces    = mesh.faces;
%         dmesh = spm_mesh_inflate(dmesh,100);
%         mesh.vertices = dmesh.vertices;
%         mesh.faces    = dmesh.faces;
%     catch
%         fprintf('Couldnt find spm_mesh_inflate: is SPM installed?\n');
%     end
% end


end