
"""
    overdamped_xvf(ontypes)

Overdamped (Euler–Maruyama) update of the position:
`v = f / zeta`, then `x += v * dt` and `xuw += v * dt`, and finally `f` is reset to zero.

# Fields

- `ontypes`: particle type (`Int`) or types (`Vector{Int}`) this evolver acts on

Requires particle fields `type`, `zeta`, `x`, `xuw`, `v` and `f`.
"""
@kwdef struct overdamped_xvf<:LocalDOFevolver
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


"""
    overdamped_pqT(ontypes)

Overdamped update of the polarity in 3D: the angular velocity is `q = T / zeta_R`, the
polarity is rotated as `p += (q × p) * dt` and normalized again, and `T` is reset to zero.

For particles that move in the xy-plane, use [`overdamped_pqT_xyc`](@ref) instead.

# Fields

- `ontypes`: particle type (`Int`) or types (`Vector{Int}`) this evolver acts on

Requires particle fields `type`, `zeta_R`, `p`, `q` and `T`.
"""
@kwdef struct overdamped_pqT<:LocalDOFevolver
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

"""
    overdamped_pqT_xyc(ontypes)

Overdamped update of the polarity for particles in the xy-plane. The polarity is rotated
exactly about the z-axis by the angle `dθ = T[3] / zeta_R * dt`, which keeps `p` in the
xy-plane with unit length. Then `q = T / zeta_R` and `T` is reset to zero.

# Fields

- `ontypes`: particle type (`Int`) or types (`Vector{Int}`) this evolver acts on

Requires particle fields `type`, `zeta_R`, `p`, `q` and `T`.
"""
@kwdef struct overdamped_pqT_xyc<:LocalDOFevolver
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




@kwdef struct overdamped_deform<:LocalDOFevolver
    ontypes::Union{Int64,Vector{Int64}}
    incompressible::Bool
end




function evolve_locally!(p_i, t, dt, dofevolver::overdamped_deform)
    if p_i.type in dofevolver.ontypes


        #Find new Lambda matrix 


        #Define shorthands
        Lambda = p_i.Lambda
        L11 = Lambda[1,1]
        L12 = Lambda[1,2]
        L21 = L12
        L22 = Lambda[2,2]

        co = p_i.p[1]
        si = p_i.p[2]

        #Note that we will not be updating Lambda0, instead we will use the rotated version instead

        semi_minor0 = p_i.semi_minor0
        semi_major0 = p_i.semi_major0

        L011 =  co^2 * semi_major0 + si^2 * semi_minor0
        L012 = si * co * (semi_major0 - semi_minor0)
        L021 = L012
        L022 =  si^2 * semi_major0 + co^2 * semi_minor0
        TL = L11 + L22

        sigma = p_i.stress
        sigma11 = sigma[1,1]
        sigma12 = sigma[1,2]
        sigma21 = sigma[2,1]
        sigma22 = sigma[2, 2]

        tau = p_i.tau

        K = p_i.K

        R0 = p_i.R

        tau_inv = 1/tau

        mu_dynamic = p_i.mu #unless mu_alpha=0
        mu_alpha = p_i.mu_alpha
        #unless
        if mu_alpha>0.
            delta_lambda2 = ( (TL)^2 - 4 * (L11*L22 - L21 * L12) ) / (TL^2)

            mu_dynamic = p_i.mu * (1 + mu_alpha * delta_lambda2 )
        end

        
        if dofevolver.incompressible

            dLmda_dt11 = -tau_inv * (L11 - L011) + tau_inv * R0/(4*mu_dynamic) * ( sigma11  - sigma22 )

            dLmda_dt22 = -tau_inv * (L22 - L022) + tau_inv * R0/(4*mu_dynamic) * ( sigma22 - sigma11 )

        else
            dLmda_dt11 = -tau_inv * (L11 - L011) + tau_inv * R0/(4*mu_dynamic * K) * ( (mu_dynamic + K ) * sigma11 + (mu_dynamic - K ) *sigma22)

            dLmda_dt22 = -tau_inv * (L22 - L022) + tau_inv * R0/(4*mu_dynamic * K) * ( (mu_dynamic + K ) * sigma22 + (mu_dynamic - K ) *sigma11)

        end

        dLmda_dt12 = -tau_inv * (L12 - L012) + tau_inv * R0/(2*mu_dynamic) *sigma12

        dLmda_dt21 = -tau_inv * (L21 - L021) + tau_inv * R0/(2*mu_dynamic) *sigma21

    


        temp_Lambda =  Lambda + SMatrix{2,2, Float64, 4}(dLmda_dt11, dLmda_dt21, dLmda_dt12, dLmda_dt22) * dt

        theta_def = 0.5*atan( temp_Lambda[1,2] + temp_Lambda[2,1] , temp_Lambda[1,1] - temp_Lambda[2,2])

        theta_current = angle.(p_i.p[1] + 1im * p_i.p[2])

        theta_def = theta_def - pi* round(Int64, (theta_def - theta_current)/pi )

        d_theta_def = theta_def - theta_current


        dθ = p_i.T[3] ./p_i.zeta_R * dt + d_theta_def
        #evolve
        pxc = p_i.p[1]
        pyc = p_i.p[2]
        pxn = cos(dθ) * pxc  - sin(dθ) *  pyc
        pyn = sin(dθ) * pxc  + cos(dθ) *  pyc
        p_i.p = SVector{3, Float64}(pxn, pyn, 0.)
        p_i.p = normalize(p_i.p )

        pxn = p_i.p[1]
        pyn = p_i.p[2]

        p_i.q = copy(p_i.T)./p_i.zeta_R
        #reinitialize
        p_i.T*=0.


        new_T = temp_Lambda[1,1]+temp_Lambda[2,2]
        Ldif = L11 - L22
        root = sqrt( 0.25*(Ldif)^2 + temp_Lambda[1,2 ] * temp_Lambda[2,1]  )

        new_semi_minor = new_T/2 - root

        new_semi_major = new_T/2 + root

        new_L11 =  pxn^2*new_semi_major+pyn^2*new_semi_minor

        new_L12 = pyn*pxn*(new_semi_major-new_semi_minor) 

        new_L21 = new_L12

        new_L22 =  pyn^2*new_semi_major+pxn^2*new_semi_minor

        p_i.Lambda = SMatrix{2,2, Float64, 4}(new_L11, new_L21, new_L12, new_L22)

        p_i.semi_major = new_semi_major

        p_i.semi_minor = new_semi_minor

        p_i.stress *=0.0

    #Update aspect ratio?
    end
    return p_i
end