function gmsh_mesh2d_write ( gmsh_filename, m, node_num, node_x, ...
  element_order, element_num, element_node )

%*****************************************************************************80
%
%% GMSH_MESH2D_WRITE writes 2D mesh data as a Gmsh mesh file.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    07 October 2014
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Christophe Geuzaine, Jean-Francois Remacle,
%    Gmsh: a three-dimensional finite element mesh generator with
%    built-in pre- and post-processing facilities,
%    International Journal for Numerical Methods in Engineering,
%    Volume 79, Number 11, pages 1309-1331, 2009.
%
%  Parameters:
%
%    Input, string GMSH_FILENAME, the name of the Gmsh file.
%
%    Input, integer M, the spatial dimension.
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Input, real NODE_X(M,NODE_NUM), the node coordinates.
%
%    Input, integer ELEMENT_ORDER, the order of the elements.
%
%    Input, integer ELEMENT_NUM, the number of elements.
%
%    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), the nodes
%    that make up each element.
%


%  Open the file.
%
  gmsh = fopen ( gmsh_filename, 'wt' );

  if ( gmsh < 0 ) 
    fprintf ( 1, '\n' );
    fprintf ( 1, 'GMSH_MESH2D_WRITE - Error!\n' );
    fprintf ( 1, '  Could not open the output file.\n' );
    error ( 'GMSH_MESH2D_WRITE - Error!' );
  end
%
%  Write the data.
%
  fprintf ( gmsh, '$MeshFormat\n' );
  fprintf ( gmsh, '  2.2 0 8\n' );
  fprintf ( gmsh, '$EndMeshFormat\n' );

  fprintf ( gmsh, '$Nodes\n' );
  fprintf ( gmsh, '  %d\n', node_num );
  for node = 1 : node_num
    fprintf ( gmsh, '  %d %g %g 0.0\n', node, node_x(1:m,node) );
  end
  fprintf ( gmsh, '$EndNodes\n' );

  if ( element_order == 3 )
    element_type = 2;
  elseif ( element_order == 6 )
    element_type = 9;
  end

  tag_num = 2;
  tag1 = 0;
  fprintf ( gmsh, '$Elements\n' );
  fprintf ( gmsh, '  %d\n', element_num );
  for element = 1 : element_num
    fprintf ( gmsh, '  %d  %d  %d  %d  %d', ...
      element, element_type, tag_num, tag1, element );
    for vertex = 1 : element_order
      fprintf ( gmsh, '  %d', element_node(vertex,element) );
    end
    fprintf ( gmsh, '\n' );
  end
  fprintf ( gmsh, '$EndElements\n' );

  fclose ( gmsh );

  return
end