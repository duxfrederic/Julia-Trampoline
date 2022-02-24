###
println("precompiling GL libraries ...")
using GLMakie, Makie, GeometryBasics
println("done")
###
function norm(a::Vector{Float64})
"""
  L2 norm of a vector a
"""
   n = 0.
   @inbounds @simd for i=1:length(a)
   n += a[i]^2
   end
   return n^0.5
end


function springForce(a::Vector{Float64}, b::Vector{Float64}, x0::Float64, k::Float64)
"""
   returns the force exerced on particle a by particle b (given that they are 
   linke by a spring). 
"""
   diff = b .- a 
   dist = norm(diff)
   dir  = diff ./ dist 
   return  k * (dist - x0) .* dir
end

function allForces(masses, links, g=-0.2, b=-0.02)
"""
   calculates the sum of the forces acting on each particle
"""
   forces = [[0.,0.,0.] for i = 1:length(masses)]
   
   @inbounds for i = 1:length(masses)
       # springs:
       others = links[i]
       @inbounds  for j = 1:length(others)
           forces[i] .+= springForce(masses[i].x, masses[others[j]].x, x0::Float64, k::Float64)
       end
       # gravity:
       forces[i] .+= [0., 0., g*masses[i].m]
       # friction:
       forces[i] .-= b .* masses[i].v
   end

   return forces
end


function VelocityVerletStep!(masses, links, dt)
"""
  Time evolution using Velocity Verlet
"""
   forces = allForces(masses, links)
   @inbounds for i = 1:length(masses)
       if  ! masses[i].fixed
          masses[i].x .+= dt .* masses[i].v + (0.5 * dt^2) .* forces[i] / masses[i].m
       end
   end
   forcesnew = allForces(masses, links)
   @inbounds for i = 1:length(masses)
       masses[i].v .+= (0.5 * dt) .* (forces[i] .+ forcesnew[i]) / masses[i].m
   end
end


function square(N, cm, sep1, sep2)
"""
 initializes the trampoline. For now it's a square of
 particles linked by springs (nearest neighbours are linked),
 outermost particles are fixed.
"""
    orig = cm .- sep1 .* (N / 2.)
    masses = Vector{NamedTuple{(:x, :v, :m, :fixed), 
                     Tuple{Vector{Float64}, 
                     Vector{Float64},
                     Float64,
                     Bool}}}(undef, 0)

    links = links = [ Vector{Int64}(undef, 0) for i=1:N^2 ]
    # create the particles
    for i=1:N
        for j=1:N
            if i == 1 || i == N || j == 1 || j == N
                fixed = true 
            else
                fixed = false
            end
            push!(masses, (x=orig.+i.*sep1.+j*sep2, v=[0.,0.,0.], m=1., fixed=fixed))
        end
    end
    # go row by row, 
    for i=1:N
        for j=1:N-1
            our = (i-1)*N + j
            other = (i-1)*N + j + 1
            
            push!(links[our], other)
            push!(links[other], our)   
        end
    end

    for i=1:N-1
        for j=1:N
            our = (i-1)*N + j
            other2 = (i*N) + j
            push!(links[our], other2)
            push!(links[other2], our)
        end
    end
    (masses, links)
end
#  creating the trampoline
N = 50 # size of the square
k = 5. # stiffness of the springs
x0 = 0. # rest length of the springs
masses, links = square(N, [0.,0.,0.], [1.,0, 0.], [0., 1., 0.])
# preparing the visualization
GLMakie.activate!()
scene = Scene(backgroundcolor=:black)
display(scene)
plotpoints = Observable([Point3f0(m.x...) for m in masses])
cam3d!(scene)
meshscatter!(scene, plotpoints, 
             color=:white, 
             markersize=0.4)
center!(scene)
# parameters
fps = 100
nframes = 1300
dt = 0.1
# main loop
for i in 1:nframes
    VelocityVerletStep!(masses, links, dt)
    if i%3 == 0
        plotpoints[] = [Point3f0(m.x...) for m in masses]
        sleep(1/fps)
    end
end
