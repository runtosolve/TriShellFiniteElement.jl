module TriShellFiniteElement

using Ferrite 

struct IP6 <: ScalarInterpolation{RefTriangle, 2}
end


function Ferrite.reference_shape_value(ip::IP6, ξ::Vec{2}, shape_number::Int)
    ξ₁ = ξ[1]
    ξ₂ = ξ[2]

    shape_number == 1 && return 1 - ξ₁ - ξ₂    
    shape_number == 2 && return ξ₁
    shape_number == 3 && return ξ₂ 
    shape_number == 4 && return 4* ξ₁ * (1 - ξ₁ - ξ₂)   
    shape_number == 5 && return 4 * ξ₁ * ξ₂ 
    shape_number == 6 && return 4 * ξ₂ * (1 - ξ₁ - ξ₂)  

    throw(ArgumentError("no shape function $shape_number for interpolation $ip"))
end


Ferrite.getnbasefunctions(::IP6) = 6

Ferrite.adjust_dofs_during_distribution(::IP6) = false


struct IP3 <: ScalarInterpolation{RefTriangle, 2}
end


function Ferrite.reference_shape_value(ip::IP3, ξ::Vec{2}, shape_number::Int)
    ξ₁ = ξ[1]
    ξ₂ = ξ[2]

    shape_number == 1 && return 1 - ξ₁ - ξ₂    
    shape_number == 2 && return ξ₁
    shape_number == 3 && return ξ₂ 

    throw(ArgumentError("no shape function $shape_number for interpolation $ip"))
end


Ferrite.vertexdof_indices(::IP3) = ((1,2,3,4,5), (6,7,8,9,10), (11,12,13,14,15))

Ferrite.getnbasefunctions(::IP3) = 3

Ferrite.adjust_dofs_during_distribution(::IP3) = false



function get_jacobian(ξ, ip, x)

    shape_number = 1
    dNdξ1 = Ferrite.reference_shape_gradient(ip, ξ, shape_number)

    shape_number = 2
    dNdξ2 = Ferrite.reference_shape_gradient(ip, ξ, shape_number)

    shape_number = 3
    dNdξ3 = Ferrite.reference_shape_gradient(ip, ξ, shape_number)


    J11 = [dNdξ1[1] dNdξ2[1] dNdξ3[1]] * [x[1][1], x[2][1], x[3][1]]
    J12 = [dNdξ1[1] dNdξ2[1] dNdξ3[1]] * [x[1][2], x[2][2], x[3][2]]
    J21 = [dNdξ1[2] dNdξ2[2] dNdξ3[2]] * [x[1][1], x[2][1], x[3][1]]
    J22 = [dNdξ1[2] dNdξ2[2] dNdξ3[2]] * [x[1][2], x[2][2], x[3][2]]


    J = [J11 J12
        J21 J22]

    return J

end


function calculate_membrane_constitutive_matrix(E, ν, t)

    G=E/(2*(1+ν))
    
    D=[    E/(1-ν^2) ν*E/(1-ν^2)    0  
        ν*E/(1-ν^2)    E/(1-ν^2)    0
                0             0    G]
    
    D = D .* t 

    return D 

end




function calculate_bending_constitutive_matrix(E, ν, t)

    D_const = E * t^3 / (12 * (1 - ν^2))
    D = D_const * [1.0  ν    0.0
                    ν    1.0  0.0
                    0.0  0.0  (1-ν)/2]

    return D 

end



function calculate_shear_constitutive_matrix(E, ν, t)

    G=E/(2*(1+ν))
    
    D=[  5/6*G*t    0.0    
        0.0        5/6*G*t]

    return D 

end

function calculate_element_membrane_stiffness_matrix(D, cv, ip_geo, ip_shape, qr, x)

    num_shape_functions = getnbasefunctions(ip_shape)

    ke = zeros(Float64, 6, 6)
    

    for q_point in 1:getnquadpoints(cv)

        ξ = qr.points[q_point]
        J = get_jacobian(ξ, ip_geo, x)
        Jinv = inv(J)
       
        B_node_all = []
    
        for i in 1:3  
       
            dNdξ1 = cv.fun_values.dNdξ[i + (q_point-1)*num_shape_functions][1]
            dNdξ2 = cv.fun_values.dNdξ[i + (q_point-1)*num_shape_functions][2]

            B_node = [dNdξ1*Jinv[1,1] + dNdξ2*Jinv[1,2]         0.0
                        0.0                                       dNdξ1*Jinv[2,1] + dNdξ2*Jinv[2,2]
                        dNdξ1*Jinv[2,1] + dNdξ2*Jinv[2,2]         dNdξ1*Jinv[1,1] + dNdξ2*Jinv[1,2]]
                
            push!(B_node_all, B_node)

        end

        B = hcat(B_node_all...)
    
        ke += B' * D * B .* det(J) .* qr.weights[q_point]

    end

    return ke 

end



function calculate_element_bending_stiffness_matrix(D, cv, ip_geo, ip_shape, qr, x)

    num_shape_functions = getnbasefunctions(ip_shape)
   
    ke = zeros(Float64, 18, 18)
   
    for q_point in 1:getnquadpoints(cv)

        ξ = qr.points[q_point]
        J = get_jacobian(ξ, ip_geo, x)
        Jinv = inv(J)
       
        B_node_all = []
    
        for i in 1:3  
       
            dNdξ1 = cv.fun_values.dNdξ[i + (q_point-1)*num_shape_functions][1]
            dNdξ2 = cv.fun_values.dNdξ[i + (q_point-1)*num_shape_functions][2]

            B_node = [0.0       0.0             dNdξ1*Jinv[1,1] + dNdξ2*Jinv[1,2]
                        0.0       -(dNdξ1*Jinv[2,1] + dNdξ2*Jinv[2,2])  0.0
                        0.0       -(dNdξ1*Jinv[1,1] + dNdξ2*Jinv[1,2])  dNdξ1*Jinv[2,1] + dNdξ2*Jinv[2,2]]
                
            push!(B_node_all, B_node)

        end

        push!(B_node_all, zeros(3, 9))
        B = hcat(B_node_all...)
  
        ke += B' * D * B .* det(J) .* qr.weights[q_point]

    end

    return ke 

end




function calculate_element_shear_stiffness_matrix(D, cv, ip_geo, ip_shape, qr, x)
   
    num_shape_functions = getnbasefunctions(ip_shape)
    
    ke = zeros(Float64, 18, 18)

    for q_point in 1:getnquadpoints(cv)

        ξ = qr.points[q_point]
        J = get_jacobian(ξ, ip_geo, x)
        Jinv = inv(J)
       
        B_node_all = []
    
        for i=1:num_shape_functions
  
            dNdξ1 = cv.fun_values.dNdξ[i + (q_point-1)*num_shape_functions][1]
            dNdξ2 = cv.fun_values.dNdξ[i + (q_point-1)*num_shape_functions][2]
                
            B_node_bending = [dNdξ1*Jinv[1,1] + dNdξ2*Jinv[1,2]         0.0     0.0
                                dNdξ1*Jinv[2,1] + dNdξ2*Jinv[2,2]         0.0     0.0]

            if i <= 3  
                
                N = Ferrite.reference_shape_value(ip_shape, ξ, i)

                B_node_shear = [0.0     0.0     N
                
                                0.0     -N      0.0]

            else 

                B_node_shear = zeros(Float64, 2, 3)
        
            end
        
            B_node = B_node_bending .+ B_node_shear

            push!(B_node_all, B_node)

        end
                  
        B = hcat(B_node_all...)
   
        ke += B' * D * B .* det(J) .* qr.weights[q_point]

    end

    return ke 

end



function elastic_stiffness_matrix!(qr1, qr3, ip3, ip6, E, ν, t, x)

    #####membrane
    cv = CellValues(qr1, ip3, ip3) 
    reinit!(cv, x)
    
    Dm = TriShellFiniteElement.calculate_membrane_constitutive_matrix(E, ν, t)

    D = Dm 
    ip_geo = ip3
    ip_shape = ip3 
    qr = qr1
    ke_m = TriShellFiniteElement.calculate_element_membrane_stiffness_matrix(D, cv, ip_geo, ip_shape, qr, x)

    ######bending
    cv = CellValues(qr1, ip3, ip3) 
    reinit!(cv, x)
    
    Db = TriShellFiniteElement.calculate_bending_constitutive_matrix(E, ν, t)

    D = Db 
    ip_geo = ip3
    ip_shape = ip3 
    qr = qr1
    ke_b = TriShellFiniteElement.calculate_element_bending_stiffness_matrix(D, cv, ip_geo, ip_shape, qr, x)

    ######shear
    cv = CellValues(qr3, ip6, ip3) 
    reinit!(cv, x)
    
    Ds = TriShellFiniteElement.calculate_shear_constitutive_matrix(E, ν, t)

    D = Ds
    ip_geo = ip3
    ip_shape = ip6 
    qr = qr3
    ke_s = TriShellFiniteElement.calculate_element_shear_stiffness_matrix(D, cv, ip_geo, ip_shape, qr, x)

    ke_bs = ke_b + ke_s

    indices = [1:10; 13; 16]
    ke_bs = ke_bs[indices, indices]

    inda=1:9
    indi=10:12
    ke_bs = ke_bs[inda,inda]-ke_bs[inda,indi]*inv(ke_bs[indi,indi])*ke_bs[indi,inda]


    ke = zeros(Float64, 15, 15)

    induv=[1 2 6 7 11 12]
    indwt=[3 4 5 8 9 10 13 14 15]
    ke[induv,induv]=ke_m
    ke[indwt,indwt]=ke_bs

    return ke 

end


function assemble_global!(K, dh, qr1, qr3, ip3, ip6, E, ν, t)
  
    assembler = start_assemble(K)
    for cell in CellIterator(dh)
        x = getcoordinates(cell) 
        ke = TriShellFiniteElement.elastic_stiffness_matrix!(qr1, qr3, ip3, ip6, E, ν, t, x)
        assemble!(assembler, celldofs(cell), ke)
    end
    return K
end



end # module TriShellFiniteElement
