function data = sort_sourcemodel(data,i)
% Sort out what source model vertices we're going to use

if      isfield(i,'pos')
            fprintf('Using supplied sourcemodel vertices\n');
            pos = i.pos;
        
elseif  isfield(i,'A') && ischar(i.A)
        fprintf('Using coords in node-file as sourcemodel\n');
        
        [~,pos] = rw_edgenode(i.A); 
        pos = pos(:,1:3);
        
%  elseif  isfield(i,'L') && ischar(i.L)
%          %IF USING A NIFTI OVERLAY AS SOURCEMODEL, NEED THIS *before* TRYING TO DO
%          %TEMPLATE SPACE CONVERSION!
%          
%          [mesh,data] = get_mesh(i,data);
%          data.mesh = mesh;
%          [~,data]  = parse_overlay(i.L,data);
%          pos       = data.sourcemodel.pos;
         
else
        fprintf('Assuming AAL90 source vertices by default\n');
        load('AAL_SOURCEMOD');
        pos  = template_sourcemodel.pos;
end