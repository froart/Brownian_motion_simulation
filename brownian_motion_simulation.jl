using GLMakie, GeometryBasics, Observables, LinearAlgebra, Distributions

# TODO make values to be depicted precisely in cm, m, etc on the screen
# TODO create function draw and zoom, which would draw objects on the screen 
# TODO make structure for solid objects

mutable struct Container
   x_min::Float64
   x_max::Float64
   y_min::Float64
   y_max::Float64
end

mutable struct Particles
   rad::Float64
   num::Int64
   pos::Observable{Array{Point2f0, 1}}
   vel::Array{Vec2f0, 1}
   collided::Array{Bool, 2}
   # Constructor
   function Particles(rad::Float64, num::Int64, boundaries::Container, maxspeed::Float64)
      @assert rad > 0 "Radius must be positive"
      @assert num > 0 "Number of particles must be positive"
      pos = Observable([Point2f0(rand(Uniform(boundaries.x_min + rad, boundaries.x_max - rad)), 
                                 rand(Uniform(boundaries.y_min + rad, boundaries.y_max - rad)))
                        for i in 1:num])
      vel = [Vec2f0(rand(Uniform(-maxspeed, maxspeed)), 
                    rand(Uniform(-maxspeed, maxspeed)))
             for i in 1:num]
      col = falses(num, num + 4)
      return new(rad, num, pos, vel, col)
   end
end

function step!(pts::Particles, container::Container, dt::Float64)

    # FIXME cumulative velocity is "jumping", which it shouldn't
    for i in 1:pts.num
        # Check for collisions with other particles and exchange velocities if needed
        for j in (i+1):pts.num
            dist = sqrt((pts.pos[][j][1] - pts.pos[][i][1])^2 + (pts.pos[][j][2] - pts.pos[][i][2])^2)
      
            # if collision
            if dist <= 2pts.rad && pts.collided[i,j] == false 
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
            if dist > 2pts.rad && pts.collided[i,j] == true
               pts.collided[i,j] = false
            end
        end

       # Check for collisions with the walls and reverse velocity if needed
       # Left
       if (pts.pos[][i][1] - pts.rad) <= container.x_min && pts.collided[i, pts.num+1] == false
          pts.vel[i] = Vec2f0( -pts.vel[i][1], pts.vel[i][2] )
          pts.collided[i, pts.num+1] = true
       end
       if (pts.pos[][i][1] - pts.rad) > container.x_min && pts.collided[i, pts.num+1] == true
          pts.collided[i, pts.num+1] = false
       end
       # Right
       if (pts.pos[][i][1] + pts.rad) >= container.x_max && pts.collided[i, pts.num+2] == false
          pts.vel[i] = Vec2f0( -pts.vel[i][1], pts.vel[i][2] )
          pts.collided[i, pts.num+2] = true
       end
       if (pts.pos[][i][1] + pts.rad) < container.x_max && pts.collided[i, pts.num+2] == true
          pts.collided[i, pts.num+2] = false
       end
       # Bottom
       if (pts.pos[][i][2] - pts.rad) <= container.y_min && pts.collided[i, pts.num+3] == false
          pts.vel[i] = Vec2f0( pts.vel[i][1], -pts.vel[i][2] )
          pts.collided[i, pts.num+3] = true
       end
       if (pts.pos[][i][2] - pts.rad) > container.y_min && pts.collided[i, pts.num+3] == true
          pts.collided[i, pts.num+3] = false
       end
       # Upper
       if (pts.pos[][i][2] + pts.rad) >= container.y_max && pts.collided[i, pts.num+4] == false
          pts.vel[i] = Vec2f0( pts.vel[i][1], -pts.vel[i][2] )
          pts.collided[i, pts.num+4] = true
       end
       if (pts.pos[][i][2] + pts.rad) < container.y_max && pts.collided[i, pts.num+4] == true
          pts.collided[i, pts.num+4] = false
       end

       # Update positions based on velocities
       pts.pos[][i] = Point2f0( pts.pos[][i][1] + dt * pts.vel[i][1], pts.pos[][i][2] + dt * pts.vel[i][2] )

    end
end

# Configuring GLMakie
GLMakie.activate!()
GLMakie.closeall()

# Defining parameters of the simulation
particles_num   = 30
speed_max       = 0.5 # (m/s)
dt              = 1.0 # (s)
particle_radius = 1.0 # (cm)
dpi = sqrt(1920^2 + 1080^2) / 15.6 # should be fixed with each new monitor
particle_radius_pxs = (particle_radius / 2.54) * dpi
# Coordinates for the borders
container = Container(0.0, 20.0, 0.0, 20.0)

# Generating particles
pts = Particles(particle_radius, particles_num, container, speed_max)

# Plot configuration
fig = Figure(size = (1200, 600))
# colsize!(fig.layout, 1, Relative(1))
# Axis of paticle box
ax  = Axis(fig[1,1][1,1], aspect = 1, limits = (container.x_min, container.x_max, container.y_min, container.y_max), width = 500, height = 500)
hidexdecorations!(ax)
hideydecorations!(ax)
# Axis of speed bar plot
ax2 = Axis(fig[1,1][1,2], aspect = 1, limits = (0, particles_num + 1, 0, 2.0speed_max), width = 500, height = 500, xlabel = "Velocity of each particle")
xs = 1:particles_num
# ESC key
quit = Observable(false)
on(events(fig).keyboardbutton) do event
   if event.action == Keyboard.press && event.key == Keyboard.escape
      quit[] = true
      notify(quit)
   end
end
# Draw lines to form a border
# Bottom
lines!(ax, [container.x_min, container.x_max], [container.y_min, container.y_min], color=:black, linewidth = 8 )
# Top
lines!(ax, [container.x_min, container.x_max], [container.y_max, container.y_max], color=:black, linewidth = 8 )
# Left
lines!(ax, [container.x_min, container.x_min], [container.y_min, container.y_max], color=:black, linewidth = 8 )
# Right
lines!(ax, [container.x_max, container.x_max], [container.y_min, container.y_max], color=:black, linewidth = 8 )
# Draw particles
scat = scatter!(ax, pts.pos, markersize = particle_radius_pxs)

fps = 60
frame_count = 0
display(fig)

while !quit[] # simulation loop

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

