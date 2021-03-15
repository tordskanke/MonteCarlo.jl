struct System{A,X,AN,B,T}
    atoms::A
    coords::X
    bonds::B
    angles::AN
    temp::T
    box_size::T
    energy::Vector{T}
    step_size::T
    permittivity::T
    n_steps::Int
    n_steps_made::Vector{Int}
    n_steps_accepted::Vector{Int}
    n_steps_rejected::Vector{Int}
end


function System(atoms,
                    coords,
                    bonds,
                    angles,
                    temp=300,
                    box_size=1.0,
                    energy=[0.0],
                    step_size=0.1,
                    permittivity=80.0,
                    n_steps=100,
                    n_steps_made=[0],
                    n_steps_accepted=[0],
                    n_steps_rejected=[0])
    A=typeof(atoms)
    X=typeof(coords)
    B=typeof(bonds)
    AN=typeof(angles)
    T=typeof(step_size)
    return System{A,X,AN,B,T}(atoms, coords, bonds , angles, temp, box_size,
                            energy, step_size, permittivity,
                            n_steps, n_steps_made, n_steps_accepted,
                            n_steps_rejected)
end
