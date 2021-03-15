using LinearAlgebra

function energy(sim::System)
    E=0.0
    β=1/sim.temp
    for (i,j,K,r0) in sim.bonds
        r=sqrt(sum(abs2,sim.coords[i]-sim.coords[j]))
        E+=harmonic_bond(β,K,r,r0)
    end
    for (i,j,k,K,ϕ0) in sim.angles
        a=sim.coords[i].-sim.coords[j]
        b=sim.coords[j].-sim.coords[k]
        ab=dot(a,b)
        ϕ=acos(ab/(norm(a)*norm(b)))
        E+=harmonic_angle(β,K,ϕ,ϕ0)
    end
return E
end

function harmonic_bond(β,K,r,r0)
    return 0.5*K/β*(r-r0)^2
end

function harmonic_angle(β,K,ϕ,ϕ0)
    return 0.5*K/β*(ϕ-ϕ0)^2
end
