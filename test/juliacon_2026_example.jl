using TriShellFiniteElement, Ferrite, LinearAlgebra

grid = let grid_2D = generate_grid(Ferrite.Triangle, (1,1), Ferrite.Vec(0.0, 0.0), Ferrite.Vec(100.0, 1000.0))
    Grid(grid_2D.cells, [Node((node.x[1], node.x[2], 0.0)) for node in grid_2D.nodes])
end

ip = Lagrange{RefTriangle,1}()
ip6 = TriShellFiniteElement.IP6()
ip3 = TriShellFiniteElement.IP3()
qr1 = QuadratureRule{RefTriangle}(1)  
qr3 = QuadratureRule{RefTriangle}(2)  


dh = DofHandler(grid)
add!(dh, :u, ip^3)
add!(dh, :θ, ip^3)
close!(dh)


cv = CellValues(qr1, ip6, ip3)

cell = first(CellIterator(dh))


 x_global = getcoordinates(cell)


T = TriShellFiniteElement.calculation_rotation_matrix(x_global)


x_local = TriShellFiniteElement.global_nodal_coords_to_planar_coords(x_global, T)


reinit!(cv, x_local)


dump(cv.fun_values)
