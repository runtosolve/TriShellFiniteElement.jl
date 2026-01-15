using Ferrite, TriShellFiniteElement 


E = 200000.0
t = 1.0
ν = 0.30

######
grid = generate_grid(Triangle, (1, 1), Vec((0.0, 0.0)), Vec((1.0, 1.0)));                      

ip6 = TriShellFiniteElement.IP6()
ip3 = TriShellFiniteElement.IP3()
qr1 = QuadratureRule{RefTriangle}(1)  
qr3 = QuadratureRule{RefTriangle}(2)  

dh = close!(add!(DofHandler(grid), :u, ip3))         


K = allocate_matrix(dh)
TriShellFiniteElement.assemble_global!(K, dh, qr1, qr3, ip3, ip6, E, ν, t)





