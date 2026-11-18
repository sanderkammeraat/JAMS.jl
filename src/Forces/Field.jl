


@kwdef struct phoresis<:FieldForce
    ontypes::Union{Int64,Vector{Int64}}
    mobility::Float64
end

function contribute_field_force!(p_i, field_j, field_indices, t, dt, rngs_particles, system, force::phoresis)
    if p_i.type in force.ontypes && field_j.type in force.ontypes


        x_index = field_indices[1]
        y_index = field_indices[2]

        ∇C = kernel_grad(p_i, field_j, x_index, y_index)
        
        p_i.f += force.mobility * ∇C


    end
    return p_i, field_j
end

@kwdef struct align_with_gradC<:FieldForce
    ontypes::Union{Int64,Vector{Int64}}
    J::Float64
end

function contribute_field_force!(p_i,field_j,field_indices, t, dt, rngs_particles, system, force::align_with_gradC)
    if p_i.type[1] in force.ontypes && field_j.type in force.ontypes


        x_index = field_indices[1]
        y_index = field_indices[2]

        ∇C = kernel_grad(p_i,field_j, x_index, y_index)

        p_i.T+= force.J*cross(p_i.p, ∇C)
    end

    return p_i, field_j
end



@kwdef struct field_propulsion_interpolation{Tif} <:FieldForce
    ontypes::Union{Int64,Vector{Int64}}
    consumption::Float64
    interpolation_function::Tif
end

function contribute_field_force!(p_i,field_j,field_indices, t, dt, rngs_particles,system, force::field_propulsion_interpolation)
    if p_i.type in force.ontypes && field_j.type in force.ontypes
        x_index = field_indices[1]
        y_index = field_indices[2]
        #print(x_index)

        cavg = 0.0

        kernel_x = p_i.kernel_x
        kernel_y = p_i.kernel_y

        for ki in eachindex(kernel_x)

            xki = kernel_x[ki]
            yki = kernel_y[ki]
            xbin_ind, ybin_ind = kernel2index_polar(p_i,field_j, x_index, y_index, xki, yki)

            
            
            if field_j.C[xbin_ind, ybin_ind]>0
                cavg+= field_j.C[xbin_ind, ybin_ind]
                field_j.Cf[xbin_ind, ybin_ind]+=-force.consumption
            end

            if field_j.C[xbin_ind, ybin_ind]<0
                field_j.Cf[xbin_ind, ybin_ind]=0
            end
            
        end
        cavg/=length(kernel_x)
        if cavg>0
            v0fact = force.interpolation_function(cavg)
            p_i.f += p_i.zeta * v0fact * SVector{3, Float64}(p_i.p[1],p_i.p[2], 0.)
    
        end
    end

    return p_i, field_j

end









#Helper functions
function kernel2index_polar(p_i, field_j, x_index, y_index, x_kernel, y_kernel)

    xki = x_kernel
    yki = y_kernel
    xkirot = xki* p_i.p[1]  -   yki * p_i.p[2]
    ykirot = xki* p_i.p[2]  +  yki * p_i.p[1]

    
    xbin_ind = x_index + floor(Int64,xkirot/field_j.lbin)

    ybin_ind = y_index + floor(Int64,ykirot/field_j.lbin)

    xlen = size(field_j.Cf)[1]
    ylen  = size(field_j.Cf)[2]

    if xbin_ind <=1
        xbin_ind =  xlen + xbin_ind -2
    elseif xbin_ind >= xlen
        xbin_ind = 2 + (xbin_ind -xlen)
    end

    if ybin_ind <=1
        ybin_ind =  ylen + ybin_ind -2
    elseif ybin_ind >= ylen
        ybin_ind = 2 + (ybin_ind -ylen)
    end


    return xbin_ind, ybin_ind
    
end

function grad(field_j, x_index, y_index)



    ∇C_x = (field_j.C[x_index+1, y_index] - field_j.C[x_index-1, y_index])/(2*field_j.lbin)

    ∇C_y = (field_j.C[x_index, y_index+1] - field_j.C[x_index, y_index-1])/(2*field_j.lbin)

    
    return SVector{3, Float64}(∇C_x, ∇C_y, 0.)
end

function kernel_grad(p_i,field_j, x_index, y_index)

    ∇C_x = 0.
    ∇C_y = 0.

    kernel_x = p_i.kernel_x
    kernel_y = p_i.kernel_y

    for ki in eachindex(kernel_x)

        xki = kernel_x[ki]
        yki = kernel_y[ki]
        xbin_ind, ybin_ind = kernel2index_polar(p_i, field_j, x_index, y_index, xki, yki)
        
        dx, dy = grad(field_j, xbin_ind, ybin_ind)
        ∇C_x+=dx
        ∇C_y+=dy

    end
    ∇C_x/=length(kernel_x)
    ∇C_y/=length(kernel_x)

    return SVector{3, Float64}(∇C_x, ∇C_y, 0.)

end



