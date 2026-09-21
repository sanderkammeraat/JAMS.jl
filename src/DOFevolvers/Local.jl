
struct overdamped_xvf<:LocalDOFevolver
    ontypes::Union{Int64,Vector{Int64}}
end

function evolve_locally!(i, current_particle_state, t, dt, dofevolver::overdamped_xvf)
    p_i = current_particle_state[i]
    if p_i.type[1] in dofevolver.ontypes

        #evolve
        current_particle_state.v[i] = p_i.f ./p_i.zeta
        
        current_particle_state.x[i]+= p_i.v * dt
        current_particle_state.xuw[i]+= p_i.v * dt
        
        #reinitialize
        current_particle_state.f[i]*=0.
    end
    return p_i
end


struct overdamped_pqT<:LocalDOFevolver
    ontypes::Union{Int64,Vector{Int64}}
end

function evolve_locally!(i, current_particle_state, t, dt, dofevolver::overdamped_pqT)
    p_i = current_particle_state[i]
    if p_i.type[1] in dofevolver.ontypes

        #evolve
        current_particle_state.q[i] = p_i.T ./p_i.zeta_R
        
        current_particle_state.p[i] += cross(current_particle_state.q[i],p_i.p) * dt
        current_particle_state.p[i] = normalize(current_particle_state.p[i])
        
        #reinitialize
        current_particle_state.T[i]*=0.
    end
    return p_i
end

struct overdamped_pqT_xyc<:LocalDOFevolver
    ontypes::Union{Int64,Vector{Int64}}
end
function evolve_locally!(i, current_particle_state, t, dt, dofevolver::overdamped_pqT_xyc)
    p_i = current_particle_state[i]
    if p_i.type[1] in dofevolver.ontypes

        dθ = p_i.T[3] ./p_i.zeta_R * dt
        #evolve
        pxc = p_i.p[1]
        pyc = p_i.p[2]
        pxn = cos(dθ) * pxc  - sin(dθ) *  pyc
        pyn = sin(dθ) * pxc  + cos(dθ) *  pyc
        current_particle_state.p[i] = SVector{3, Float64}(pxn, pyn, 0.)
        current_particle_state.p[i] = normalize(current_particle_state.p[i])

        current_particle_state.q[i] = copy(p_i.T)./p_i.zeta_R
        #reinitialize
        current_particle_state.T[i]*=0.
    end
    return p_i
end