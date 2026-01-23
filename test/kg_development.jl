
using Ferrite, TriShellFiniteElement 


E = 200000.0
t = 1.0
ν = 0.30

######
grid = generate_grid(Triangle, (1, 1), Vec((0.0, 0.0)), Vec((1.0, 1.0)));                      

ip3 = TriShellFiniteElement.IP3()
qr1 = QuadratureRule{RefTriangle}(1)  


ip_geo = ip3 
ip_shape = ip3
qr = qr1

cv = CellValues(qr1, ip3, ip3) 

cell = first(CellIterator(grid))

x = getcoordinates(cell) 
reinit!(cv, x)


    ip_geo = TriShellFiniteElement.IP3() 
    ip_shape = TriShellFiniteElement.IP3()
    qr = QuadratureRule{RefTriangle}(1)  

    cv = CellValues(qr, ip_shape, ip_geo) 


    kgx, kgy, kgxy = calculate_element_geometric_stiffness_matrix(cv, ip_geo, ip_shape, qr, x)




i1=elems_row(1); i2=elems_row(2); i3=elems_row(3);
P1=nodes(i1,1:3); P2=nodes(i2,1:3); P3=nodes(i3,1:3);


norm_vec=cross(P2-P1,P3-P1);
norm_vec=norm_vec/norm(norm_vec);
j3=norm_vec
j1=(P2-P1)/norm(P2-P1)
j2=cross(j3,j1)
T=[j1' j2' j3']


P1a=T'*(P1-P1)';
P2a=T'*(P2-P1)';
P3a=T'*(P3-P1)';
X1=P1a(1); Y1=P1a(2);
X2=P2a(1); Y2=P2a(2);
X3=P3a(1); Y3=P3a(2);
