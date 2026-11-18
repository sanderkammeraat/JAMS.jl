@kwdef struct overdamped_CCvCf <: FieldDOFevolver
    ontypes::Union{Int64,Vector{Int64}}
end

function evolve_field!(field, t, dt, dofevolver::overdamped_CCvCf)

    if field.type in dofevolver.ontypes

        Threads.@threads for j in 2:size(field.C,2)-1

            @inbounds for i in 2:size(field.C,1)-1

                field.Cv[i,j] = field.Cf[i,j]
                field.C[i,j] += field.Cv[i,j]*dt
                #reinitalize
                field.Cf[i,j]= 0.0
            end
        end
    end
    return field
end