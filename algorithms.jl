using Plots
using StaticArrays
function simulate!(sim::System)
  MonteCarlo(sim)
end

function write_xyz(pre,iter::Int,freq::Int)
  if iter%freq==0
    file = open(pre*string(trunc(Int,floor(iter/freq)))*".xyz", "w")
    write(file,string(length(atoms),"\n"))
    write(file,"Monte-carlo.jl generated\n")
    for i in sim.coords
      write(file, "C   "*string(round(i[1],digits=2),"  ",round(i[2],digits=2),"  ",round(i[3],digits=2)));
      write(file,"\n")
    end
    close(file)
  end
end

function MonteCarlo!(sim::System)
  β=1/sim.temp
  while sim.n_steps_made[1] < sim.n_steps
      sim.n_steps_made[1]+=1
      particle=rand(1:length(sim.coords))
      q=rand(1:length(sim.coords[1]))
      step=sim.step_size.*(2.0.*rand(3).-1.0)
      sim.coords[particle]+=step
      trial_energy=energy(sim)
      if trial_energy < sim.energy[1] || exp(-(trial_energy-sim.energy[1])*β) > rand(1)[1]
        sim.energy[1] = trial_energy
        sim.n_steps_accepted[1]+=1
      else
        sim.coords[particle]-=step
        sim.n_steps_rejected[1]+=1
      end
      write_xyz("temp",sim.n_steps_made[1],100000)
    end
#println(sum(lengths)/sim.n_steps)
#println([sqrt(sum(abs2,sim.coords[i]-sim.coords[i+1])) for i in 1:length(sim.coords)-1])
end
atoms=[i for i in 1:64]
X=[ @SVector [i*1.5,0.0,0.0] for i in atoms]
bonds=[(i,i+1,1,1.0) for i in atoms[1:end-1]]
angles=[(i,i+1,i+2,10,0.0) for i in atoms[1:end-2]]
sim=System(atoms,X,bonds,angles,50,1.0,[0.0],0.31,80.0,100000)

@time MonteCarlo!(sim)

println(sim.n_steps_accepted/sim.n_steps)
plot([sim.coords[i][1] for i in 1:length(sim.coords)],[sim.coords[i][2] for i in 1:length(sim.coords)],seriestype=scatter,
ylims=(-20.0,20.0))
