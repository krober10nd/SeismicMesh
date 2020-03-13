function p=protate3d(p,phi)

rotx =[1 0 0; 0 cos(phi) -sin(phi) ; 0 sin(phi) cos(phi)] ;
roty =[cos(phi) 0 sin(phi) ; 0 1 0 ; -sin(phi) 0  cos(phi)] ;
rotz =[cos(phi) -sin(phi) 0 ; sin(phi) cos(phi) 0 ; 0 0 1] ;

p=  p*rotx; 
