# vtkwrite
vtkwrite writes 3D Matlab array into VTK file format

Paraview is a powerful open-source software for visualization of large 3D dataset. It offers more options, details and much better performance than built-in Matlab 3D visualization modules. This function is an integration of several previous submissions regarding export of 3D data into VTK format. The function can save multiple vector and scalar field of the same size into a single VTK-formatted file to be viewed in ParaView. It can also export line or polygon objects. To maximize compatibility between different operating system, numerical data is by default saved in ascii format with precision of 3 digits behind decimal point. User can specify precision , e.g. vtkwrite('execute', 'structured_points', 'mri', D, 'precision, 5), instead of vtkwrite('execute', 'structured_points', 'mri', D). User can also choose to save numerical data in 'float' data type ( this option is not available for POLYDATA dataset type) by adding 'binary' to the command, e.g. vtkwrite('execute', 'structured_points', 'mri', D, binary). 
  
Here are some example usages. Just type in the following codes.  

Example 1 : export 3D array (typical of an image volume )

    load mri  
    D = squeeze(D);  
    vtkwrite('mri.vtk', 'structured_points', 'mri', D)  

if you've already setup system path to include the folder containing the ParaView binary, you can invoke ParaView directly by 

    vtkwrite('execute', 'structured_points', 'mri', D)

In this case, a file named 'matlab_export.vtk' is saved and passed on to ParaView.

Example 2 : export 3D vector and scalar field 

    load wind   
    [cu,cv,cw] = curl(x, y, z, u, v, w);       
    div = divergence(x, y, z, u, v, w);     
    vtkwrite('wind.vtk', 'structured_grid', x, y, z, ...   
    'vectors', 'vector_field', u, v, w, 'vectors', 'vorticity', cu, cv, cw, 'scalars', 'divergence', div); 

Usage is very similar for unstructured 3D data. Just change 'structured_grid' to 'unstructured_grid'. For example : 

    vtkwrite('random_vector.vtk', 'unstructured_grid',x,y,z, 'vectors','random_vector',u,v,w,) 

Example 3 : export 3D line 

    x = 1:100; 
    y = sin(x); 
    z = sqrt(x);     
    vtkwrite('execute','polydata','lines',x,y,z);     

optionally, user can specify precision of data output(default is 3). for example: 

    vtkwrite('line.vtk','polydata','lines',x,y,z,'precision',5);     

Example 4 : export triangular data 

    [x,y,z] = peaks(100);     
    z = .4*z;     
    tri = delaunay(x,y);     
    vtkwrite('peaks.vtk','polydata','triangle',x,y,z,tri);    

Example 5 : export tetrahedral data 

    d = [-1 1];     
    [x, y, z] = meshgrid(d, d, d);     
    DT = delaunayTriangulation(x(:), y(:), z(:));     
    vtkwrite('execute', 'polydata','tetrahedron', x, y, z, DT.ConnectivityList);    

The usage of vector and scalar field input is the same as built-in quiver3 and scatter3 function, where x,y,z specifies the coordinates of the grid points and u,v,w the vector components. For multiple data array entry, they must have the same number of grid points. x,y,z only need to be specified once, followed by combination of [keyword, title, data]. You can add however many data arrays as you want. In the example above, you'll end up with a single VTK file with three data arrays of the same number of grid points. In the sample screenshot, the color of the arrow represents the magnitude of divergence.

In Paraview, press 'ctrl+o' to open file. Then on properties page immediately below pipeline browser, click 'Apply'. For line and polygon data, you should already see the correct representation of data. For 3D vector field, you have to further click on 'glyph' in the common toolbar, then choose the glyph object in the pipeline browser and click 'Apply'.

See file for more details on usage.
