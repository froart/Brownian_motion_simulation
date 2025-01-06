using GLMakie, GeometryBasics, Observables, LinearAlgebra, Distributions

# TODO create function draw and zoom, which would draw objects on the screen 
# TODO make structure for solid objects

mutable struct Container
   x_min::Float64
   x_max::Float64
   y_min::Float64
   y_max::Float64
end

mutable struct Particles
   diam::Float64
   num::Int64
   pos::Observable{Array{Point2f0, 1}}
   vel::Array{Vec2f0, 1}
   collided::Array{Bool, 2}
   # Constructor
   function Particles(diam::Float64, num::Int64, boundaries::Container, maxspeed::Float64)
      @assert diam > 0 "Radius must be positive"
      @assert num > 0 "Number of particles must be positive"
      pos = Observable([Point2f0(rand(Uniform(boundaries.x_min + diam/2, boundaries.x_max - diam/2)), 
                                 rand(Uniform(boundaries.y_min + diam/2, boundaries.y_max - diam/2)))
                        for i in 1:num])
      vel = [Vec2f0(rand(Uniform(-maxspeed, maxspeed)),
                    rand(Uniform(-maxspeed, maxspeed)))
             for i in 1:num]
      col = falses(num, num + 4)
      return new(diam, num, pos, vel, col)
   end
end

function step!(pts::Particles, container::Container, dt::Float64)

    # TODO make computation via law of momentum conservation
    # FIXME cumulative velocity is "jumping", which it shouldn't
    for i in 1:pts.num
        # Check for collisions with other particles and exchange velocities if needed
        for j in (i+1):pts.num
            dist = sqrt((pts.pos[][j][1] - pts.pos[][i][1])^2 + (pts.pos[][j][2] - pts.pos[][i][2])^2)
      
            # if collision
            if dist <= 2pts.diam/2 && pts.collided[i,j] == false
               # unit vector from centre to centre
               dij      =  Vec2f0( pts.pos[][j][1] - pts.pos[][i][1], pts.pos[][j][2] - pts.pos[][i][2] )
               dij_unit =  dij / sqrt( dij[1]^2 + dij[2]^2 )
               dji_unit = -dij_unit
               # projection of vi & vj unto dij_unit & dji_unit accordingly
               proj_d_vi = dij_unit * dot(pts.vel[i], dij_unit) / ( dij_unit[1]^2 + dij_unit[2]^2 ) # formally speaking
               proj_d_vj = dji_unit * dot(pts.vel[j], dji_unit) / ( dji_unit[1]^2 + dji_unit[2]^2 ) # formally speaking
               ###   ELASTIC COLLISION
               pts.vel[j] = Vec2f0( pts.vel[j][1] + proj_d_vi[1] - proj_d_vj[1],
                                    pts.vel[j][2] + proj_d_vi[2] - proj_d_vj[2] )

               pts.vel[i] = Vec2f0( pts.vel[i][1] + proj_d_vj[1] - proj_d_vi[1],
                                    pts.vel[i][2] + proj_d_vj[2] - proj_d_vi[2] )
               pts.collided[i,j] = true
            end

            # When a particle escape another particle
            if dist > pts.diam && pts.collided[i,j] == true
               pts.collided[i,j] = false
            end
        end

       # Check for collisions with the walls and reverse velocity if needed
       # Left
       if (pts.pos[][i][1] - pts.diam/2) <= container.x_min && pts.collided[i, pts.num+1] == false
          pts.vel[i] = Vec2f0( -pts.vel[i][1], pts.vel[i][2] )
          pts.collided[i, pts.num+1] = true
       end
       if (pts.pos[][i][1] - pts.diam/2) > container.x_min && pts.collided[i, pts.num+1] == true
          pts.collided[i, pts.num+1] = false
       end
       # Right
       if (pts.pos[][i][1] + pts.diam/2) >= container.x_max && pts.collided[i, pts.num+2] == false
          pts.vel[i] = Vec2f0( -pts.vel[i][1], pts.vel[i][2] )
          pts.collided[i, pts.num+2] = true
       end
       if (pts.pos[][i][1] + pts.diam/2) < container.x_max && pts.collided[i, pts.num+2] == true
          pts.collided[i, pts.num+2] = false
       end
       # Bottom
       if (pts.pos[][i][2] - pts.diam/2) <= container.y_min && pts.collided[i, pts.num+3] == false
          pts.vel[i] = Vec2f0( pts.vel[i][1], -pts.vel[i][2] )
          pts.collided[i, pts.num+3] = true
       end
       if (pts.pos[][i][2] - pts.diam/2) > container.y_min && pts.collided[i, pts.num+3] == true
          pts.collided[i, pts.num+3] = false
       end
       # Upper
       if (pts.pos[][i][2] + pts.diam/2) >= container.y_max && pts.collided[i, pts.num+4] == false
          pts.vel[i] = Vec2f0( pts.vel[i][1], -pts.vel[i][2] )
          pts.collided[i, pts.num+4] = true
       end
       if (pts.pos[][i][2] + pts.diam/2) < container.y_max && pts.collided[i, pts.num+4] == true
          pts.collided[i, pts.num+4] = false
       end

       # Update positions based on velocities
       pts.pos[][i] = Point2f0( pts.pos[][i][1] + dt * pts.vel[i][1], pts.pos[][i][2] + dt * pts.vel[i][2] )

    end
end

function to_pixels(value_in_cm::Float64)
   dpi = sqrt(1920^2 + 1080^2) / 15.6 # should be fixed with each new monitor
   return (value_in_cm / 2.54) * dpi
end

# Configuring GLMakie
GLMakie.activate!()
GLMakie.closeall()

# Defining parameters of the simulation
particles_num   = 1000
speed_max       = 1.5 # (m/s)
dt              = 0.1 # (s)
particle_diameter = 0.25 # (cm)
window_width    = 25.0 # (cm)
window_height   = 25.0 # (cm)
container_size  = 10.0 # (cm)

# Coordinates for the borders
container = Container(0.0, container_size, 0.0, container_size)

# Generating particles
pts = Particles(particle_diameter, particles_num, container, speed_max)

# Plot configuration
fig = Figure(size = (to_pixels(window_width), to_pixels(window_height)))

# Axis of paticle box
ax  = Axis(fig[1,1][1,1], aspect = 1, limits = (container.x_min, container.x_max, container.y_min, container.y_max), width = to_pixels(container_size), height = to_pixels(container_size))

hidexdecorations!(ax)
hideydecorations!(ax)
# Bottom
lines!(ax, [container.x_min, container.x_max], [container.y_min, container.y_min], color=:black, linewidth = 3 )
# Top
lines!(ax, [container.x_min, container.x_max], [container.y_max, container.y_max], color=:black, linewidth = 3 )
# Left
lines!(ax, [container.x_min, container.x_min], [container.y_min, container.y_max], color=:black, linewidth = 3 )
# Right
lines!(ax, [container.x_max, container.x_max], [container.y_min, container.y_max], color=:black, linewidth = 3 )
# FIXME change the size of the axis
ax2 = Axis(fig[1,1][1,2], aspect = 1, limits = (0, particles_num + 1, 0, 2.0speed_max), width = to_pixels(container_size), height = to_pixels(container_size), xlabel = "Velocity of each particle")

# Draw particles
particle_diameter_px = to_pixels(particle_diameter)
scat = scatter!(ax, pts.pos, markersize = particle_diameter_px, markerspace = :pixel)

quit = Observable(false)
on(events(fig).keyboardbutton) do event
   if event.action == Keyboard.press && event.key == Keyboard.escape # ESC key
      quit[] = true
      notify(quit)
   end
end

fps = 60
frame_count = 0
display(fig)
xs = 1:particles_num

while !quit[] # live simulation loop

      step!(pts, container, dt)
      notify(pts.pos)

      ys = [ sqrt(pts.vel[i][1]^2 + pts.vel[i][2]^2) for i in xs ]
      empty!(ax2)
      barplot!(ax2, xs, ys, color = :blue)

      # Calculate the cumulative velocity
      ax2.subtitle = "Cumulative velocity: $(sum(ys))"

      sleep(1/fps)
      frame_count = frame_count + 1
end

GLMakie.closeall()

