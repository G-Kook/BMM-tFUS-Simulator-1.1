% Made by G. Kook
% Reference: https://itectec.com/matlab/matlab-how-to-export-3d-image-from-matlab-and-import-it-into-maya/

function voxelToStl(grid_x, grid_y, grid_z, model_voxelized, dir_out, compression, voxelSmooth)
if ndims(model_voxelized)~=3
    error('Model must be 3-dimensional');
end
if size(model_voxelized,1)~=numel(grid_x) | size(model_voxelized,2)~=numel(grid_y) | size(model_voxelized,3)~=numel(grid_z)
	error('Grid and model must be in the same size');
end

%% (1/4) remove zeros to reduce matrix size
Nx=numel(grid_x); Ny=numel(grid_y); Nz=numel(grid_z);
[idx,idy,idz] = ind2sub([Nx Ny Nz], find(model_voxelized));
model_reduced = model_voxelized(min(idx,[],'all'):max(idx,[],'all'), min(idy,[],'all'):max(idy,[],'all'), min(idz,[],'all'):max(idz,[],'all'));
grid_x = grid_x(min(idx,[],'all'):max(idx,[],'all'));
grid_y = grid_y(min(idy,[],'all'):max(idy,[],'all'));
grid_z = grid_z(min(idz,[],'all'):max(idz,[],'all'));
%% (2/4) optimize the model before export
model_export = 2*(model_reduced-0.5); % convert 0/1 to -1/1
if voxelSmooth
    model_export = smooth3(model_export); % smooth
end
if compression > 1
    [grid_Y,grid_X,grid_Z] = meshgrid(grid_y,grid_x,grid_z);
    model_export = reducevolume(model_export,compression);
    grid_X = reducevolume(grid_X,compression);
    grid_Y = reducevolume(grid_Y,compression);
    grid_Z = reducevolume(grid_Z,compression);
    grid_x = reshape(grid_X(:,1,1),1,[]);
    grid_y = reshape(grid_Y(1,:,1),1,[]);
    grid_z = reshape(grid_Z(1,1,:),1,[]);
end
model_export = padarray(model_export,[3 3 3],0,'both'); % because if the model exists right at the grid edge, the surface is not closed in the stl file.
dx=grid_x(2)-grid_x(1); dy=grid_y(2)-grid_y(1); dz=grid_z(2)-grid_z(1); 
grid_x = padarray(grid_x,[0 3],0,'both'); grid_y = padarray(grid_y,[0 3],0,'both'); grid_z = padarray(grid_z,[0 3],0,'both');
Nx=numel(grid_x); Ny=numel(grid_y); Nz=numel(grid_z);
grid_x(1:3) = grid_x(4)-3*dx : dx : grid_x(4)-dx;
grid_x(Nx-2:Nx) = grid_x(Nx-3)+dx : dx : grid_x(Nx-3)+3*dx;
grid_y(1:3) = grid_y(4)-3*dy : dy : grid_y(4)-dy;
grid_y(Ny-2:Ny) = grid_y(Ny-3)+dy : dy : grid_y(Ny-3)+3*dy;
grid_z(1:3) = grid_z(4)-3*dz : dz : grid_z(4)-dz;
grid_z(Nz-2:Nz) = grid_z(Nz-3)+dz : dz : grid_z(Nz-3)+3*dz;

%% (3/4) convert voxel to mesh
[grid_Y,grid_X,grid_Z] = meshgrid(grid_y,grid_x,grid_z);
mesh = isosurface(grid_X,grid_Y,grid_Z,model_export,0);
TR = triangulation(mesh.faces,mesh.vertices);
% figure('color','w'), h=trimesh(TR), axis equal
% set(h,'FaceColor','w','EdgeColor','b')

% %% (4/4) convert mesh to stl
stlwrite(TR,dir_out);

% end