using GLMakie, GeometryBasics, Observables, LinearAlgebra, Distributions

function step!(positions, velocities, particle_diam, num_particles, containerside, dt, already_collided, walls_coord)

    # FIX cumulative velocity is "jumping", which it shouldn't
    for i in 1:num_particles
        # Check for collisions with other particles and exchange velocities if needed
        for j in (i+1):num_particles
            dist = sqrt(( positions[][j][1] - positions[][i][1] )^2 + ( positions[][j][2] - positions[][i][2] )^2)
      
            # if collision
            if dist <= particle_diam && already_collided[i,j] == false 
               # unit vector from centre to centre
               dij      =  Vec2f0( positions[][j][1] - positions[][i][1], positions[][j][2] - positions[][i][2] )
               dij_unit =  dij / sqrt( dij[1]^2 + dij[2]^2 )
               dji_unit = -dij_unit
               # projection of vi & vj unto dij_unit & dji_unit accordingly
               proj_d_vi = dij_unit * dot(velocities[i], dij_unit) / ( dij_unit[1]^2 + dij_unit[2]^2 ) # formally speaking
               proj_d_vj = dji_unit * dot(velocities[j], dji_unit) / ( dji_unit[1]^2 + dji_unit[2]^2 ) # formally speaking
               ###   ELASTIC COLLISION
               velocities[j] = Vec2f0( velocities[j][1] + proj_d_vi[1] - proj_d_vj[1],
                                       velocities[j][2] + proj_d_vi[2] - proj_d_vj[2] )

               velocities[i] = Vec2f0( velocities[i][1] + proj_d_vj[1] - proj_d_vi[1],
                                       velocities[i][2] + proj_d_vj[2] - proj_d_vi[2] )
               already_collided[i,j] = true
            end

            # if a particle escape another particle
            if dist > particle_diam && already_collided[i,j] == true
               already_collided[i,j] = false
            end
        end

       # Check for collisions with the walls and reverse velocity if needed
       # Left
       if (positions[][i][1] - particle_diam[]/2) <= walls_coord[:xmin] && already_collided[i, num_particles+1] == false
          velocities[i] = Vec2f0( -velocities[i][1], velocities[i][2] )
          already_collided[i, num_particles+1] = true
       end
       if (positions[][i][1] - particle_diam[]/2) > walls_coord[:xmin] && already_collided[i, num_particles+1] == true
          already_collided[i, num_particles+1] = false
       end
       # Right
       if (positions[][i][1] + particle_diam[]/2) >= walls_coord[:xmax] && already_collided[i, num_particles+2] == false
          velocities[i] = Vec2f0( -velocities[i][1], velocities[i][2] )
          already_collided[i, num_particles+2] = true
       end
       if (positions[][i][1] + particle_diam[]/2) < walls_coord[:xmax] && already_collided[i, num_particles+2] == true
          already_collided[i, num_particles+2] = false
       end
       # Bottom
       if (positions[][i][2] - particle_diam[]/2) <= walls_coord[:ymin] && already_collided[i, num_particles+3] == false
          velocities[i] = Vec2f0( velocities[i][1], -velocities[i][2] )
          already_collided[i, num_particles+3] = true
       end
       if (positions[][i][2] - particle_diam[]/2) > walls_coord[:ymin] && already_collided[i, num_particles+3] == true
          already_collided[i, num_particles+3] = false
       end
       # Upper
       if (positions[][i][2] + particle_diam[]/2) >= walls_coord[:ymax] && already_collided[i, num_particles+4] == false
          velocities[i] = Vec2f0( velocities[i][1], -velocities[i][2] )
          already_collided[i, num_particles+4] = true
       end
       if (positions[][i][2] + particle_diam[]/2) < walls_coord[:ymax] && already_collided[i, num_particles+4] == true
          already_collided[i, num_particles+4] = false
       end

       # Update positions based on velocities
       positions[][i] = Point2f0( positions[][i][1] + dt * velocities[i][1], positions[][i][2] + dt * velocities[i][2] )

    end
end

function run_simulation()

    GLMakie.activate!()
    GLMakie.closeall()

    num_particles = 30
    zoomframeside = Observable(0.4)
    containerside = 2
    particle_diam = 0.1
    # FIX make scattered particles and axis of fixed sizes for real measure representation
    diam_for_scat = Observable( particle_diam * 270 )
    maxspeed      = 0.05
    dt            = 1
    # Coordinates for the border
    walls_coord = ( xmin = 0, xmax = containerside, ymin = 0, ymax = containerside)

    # Define initial positions and velocities for the particles
    positions  = Observable([ Point2f0( 
                              rand(Uniform(walls_coord[:xmin] + particle_diam/2, walls_coord[:xmax] - particle_diam/2)), 
                              rand(Uniform(walls_coord[:ymin] + particle_diam/2, walls_coord[:ymax] - particle_diam/2)))
		                               for _ in 1:num_particles ])
    velocities = [ Vec2f0( maxspeed * rand() * rand([-1,1]), 
	                         maxspeed * rand() * rand([-1,1])) 
                           for _ in 1:num_particles ]

    # Plot configuration
    fig = Figure(size = (1200, 600))
    # colsize!(fig.layout, 1, Relative(1))
    # Axis of paticle box
    ax  = Axis(fig[1,1][1,1], aspect = 1, limits = (0, containerside, 0, containerside), width = 500, height = 500)
    hidexdecorations!(ax)
    hideydecorations!(ax)
    # Axis of speed bar plot
    ax2 = Axis(fig[1,1][1,2], aspect = 1, limits = (0, num_particles+1, 0, 2.0maxspeed), width = 500, height = 500, xlabel = "Velocity of each particle")
    xs = 1:num_particles
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
    lines!(ax, [walls_coord[:xmin], walls_coord[:xmax]], [walls_coord[:ymin], walls_coord[:ymin]], color=:black, linewidth = 8 )
    # Top
    lines!(ax, [walls_coord[:xmin], walls_coord[:xmax]], [walls_coord[:ymax], walls_coord[:ymax]], color=:black, linewidth = 8 )
    # Left
    lines!(ax, [walls_coord[:xmin], walls_coord[:xmin]], [walls_coord[:ymin], walls_coord[:ymax]], color=:black, linewidth = 8 )
    # Right
    lines!(ax, [walls_coord[:xmax], walls_coord[:xmax]], [walls_coord[:ymin], walls_coord[:ymax]], color=:black, linewidth = 8 )
    # Draw particles
    scat = scatter!(ax, positions, markersize = diam_for_scat)

    fps = 60
    frame_count = 0
    display(fig)
    already_collided = falses(num_particles, num_particles+4)

    while !quit[] # simulation loop

          step!(positions, velocities, particle_diam, num_particles, containerside, dt, already_collided, walls_coord)
          notify(positions)

          ys = [ sqrt(velocities[i][1]^2 + velocities[i][2]^2) for i in xs ]
          empty!(ax2)
          barplot!(ax2, xs, ys, color = :blue)

          # Calculate the cumulative velocity
          ax2.subtitle = "Cumulative velocity: $(sum(ys))"

          sleep(1/fps)
          frame_count = frame_count + 1
    end

    GLMakie.closeall()
end

