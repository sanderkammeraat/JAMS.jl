
struct overdamped_xvf<:LocalDOFevolver
    ontypes::Union{Int64,Vector{Int64}}
end

function evolve_locally!(p_i, t, dt, dofevolver::overdamped_xvf)
    if p_i.type in dofevolver.ontypes

        #evolve
        p_i.v = p_i.f ./p_i.zeta
        
        p_i.x += p_i.v * dt
        p_i.xuw += p_i.v * dt
        
        #reinitialize
        p_i.f *=0.
    end
    return p_i
end


struct overdamped_pqT<:LocalDOFevolver
    ontypes::Union{Int64,Vector{Int64}}
end

function evolve_locally!(p_i, t, dt, dofevolver::overdamped_pqT)
    if p_i.type in dofevolver.ontypes

        #evolve
        p_i.q = p_i.T ./p_i.zeta_R
        
        p_i.p += cross(p_i.q,p_i.p) * dt
        p_i.p = normalize(p_i.p)
        
        #reinitialize
        p_i.T*=0.
    end
    return p_i
end

struct overdamped_pqT_xyc<:LocalDOFevolver
    ontypes::Union{Int64,Vector{Int64}}
end
function evolve_locally!(p_i, t, dt, dofevolver::overdamped_pqT_xyc)
    if p_i.type in dofevolver.ontypes

        dθ = p_i.T[3] ./p_i.zeta_R * dt
        #evolve
        pxc = p_i.p[1]
        pyc = p_i.p[2]
        pxn = cos(dθ) * pxc  - sin(dθ) *  pyc
        pyn = sin(dθ) * pxc  + cos(dθ) *  pyc
        p_i.p = SVector{3, Float64}(pxn, pyn, 0.)
        p_i.p = normalize(p_i.p )

        p_i.q = copy(p_i.T)./p_i.zeta_R
        #reinitialize
        p_i.T*=0.
    end
    return p_i
end