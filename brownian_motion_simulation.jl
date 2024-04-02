using GLMakie, GeometryBasics, Observables, LinearAlgebra

function run_simulation()

    GLMakie.activate!()
    GLMakie.closeall()

    num_particles = 30
    zoomframeside = Observable(0.4)
    containerside = 2 # cm
    particle_diam = 0.1 # cm
    # FIX make scattered particles and axis of fixed sizes for real measure representation
    diam_for_scat = Observable( particle_diam * 270 )
    maxspeed      = 0.05
    dt            = 1

    # Define initial positions and velocities for the particles
    positions  = Observable([ Point2f0( ( -particle_diam/2 + containerside ) * rand() + particle_diam/2, 
		    	                              ( -particle_diam/2 + containerside ) * rand() + particle_diam/2 ) 
		                                       for _ in 1:num_particles ])
    velocities = [ Vec2f0( maxspeed * rand() * rand([-1,1]), 
	                         maxspeed * rand() * rand([-1,1])) 
                           for _ in 1:num_particles ]

    # Plot configuration
    fig = Figure(size = (600, 600))
    colsize!(fig.layout, 1, Relative(1))
    fig[1, 1] = grid1 = GridLayout()
    # Axis of paticle box
    ax  = Axis(grid1[1, 1], aspect = 1, limits = (0, containerside, 0, containerside))
    hidexdecorations!(ax)
    hideydecorations!(ax)
    fig[2, 1] = grid2 = GridLayout()
    # Axis of speed bar plot
    ax2 = Axis(grid1[1,2], aspect = 1, limits = (0, num_particles+1, 0, maxspeed))
    xs = 1:num_particles
    # Velocity sum label
    label = Label(grid2[1, 1], "")
    # Escape button
    quit   = Observable(false)
    on(events(fig).keyboardbutton) do event
       if event.action == Keyboard.press && event.key == Keyboard.escape
          quit[] = true
          notify(quit)
       end
    end
    # Coordinates for the border
    xmin, xmax, ymin, ymax = 0, containerside, 0, containerside
    # Draw lines to form a border
    lines!(ax, [xmin, xmax], [ymin, ymin], color=:black, linewidth = 8 ) # Bottom
    lines!(ax, [xmin, xmax], [ymax, ymax], color=:black, linewidth = 8 ) # Top
    lines!(ax, [xmin, xmin], [ymin, ymax], color=:black, linewidth = 8 ) # Left
    lines!(ax, [xmax, xmax], [ymin, ymax], color=:black, linewidth = 8 ) # Right
    # Draw particles
    scat = scatter!(ax, positions, markersize = diam_for_scat)

    fps = 60
    frame_count = 0
    display(fig)
    pair_check_prev = []

    while !quit[] # simulation loop
          step!(positions, velocities, particle_diam, num_particles, containerside, dt, pair_check_prev)
          scat[1][] = positions[]
          notify(positions)

          # Calculate the cumulative velocity
          tot_vel = 0
          for i in 1:num_particles
              tot_vel = tot_vel + sqrt(velocities[i][1] ^ 2 + velocities[i][2] ^ 2 )
          end
          label.text = "Sum: $(tot_vel)"

          ys = [ sqrt(velocities[i][1]^2 + velocities[i][2]^2) for i in xs ]
          empty!(ax2)
          barplot!(ax2, xs, ys, color = :blue)

          #=
          if frame_count > 520
             posx = positions[][1][1]
             posy = positions[][1][2]
             limits!(ax, posx - zoomframeside[]/2, posx + zoomframeside[]/2, posy - zoomframeside[]/2, posy + zoomframeside[]/2)
          end
          if frame_count == 521
             global diam_for_scat[] = diam_for_scat[] * containerside / zoomframeside[] 
          end
          =#

          sleep(1/fps)
          frame_count = frame_count + 1
    end

    GLMakie.closeall()
end

function step!(positions, velocities, particle_diam, num_particles, containerside, dt, pair_check_prev)

    pair_check = []
    for i in 1:num_particles
        # Check for collisions with other particles and exchange velocities if needed
        for j in 1:num_particles
            # if the same particle or if we have already balanced the particular pair of particles in this step
            if j == i || ((j, i) in pair_check) continue end
            # if collision
            if sqrt(( positions[][j][1] - positions[][i][1] )^2 + ( positions[][j][2] - positions[][i][2] )^2) <= particle_diam && !((i, j) in pair_check_prev) && !((j, i) in pair_check_prev)
               # unit vector from centre to centre
               dij      = Vec2f0( positions[][j][1] - positions[][i][1], positions[][j][2] - positions[][i][2] )
               dij_unit =  dij / sqrt( dij[1]^2 + dij[2]^2 )
               dji_unit = -dij_unit
               # projection of vi & vj unto dij_unit & dji_unit accordingly
               proj_d_vi_magn = dot(velocities[i], dij_unit) / ( dij_unit[1]^2 + dij_unit[2]^2 ) # formally speaking
               proj_d_vi      = proj_d_vi_magn * dij_unit
               proj_d_vj_magn = dot(velocities[j], dji_unit) / ( dji_unit[1]^2 + dji_unit[2]^2 ) # formally speaking
               proj_d_vj      = proj_d_vj_magn * dji_unit
               ###   ELASTIC COLLISION
               velocities[j] = Vec2f0( velocities[j][1] + proj_d_vi[1] - proj_d_vj[1],
                                       velocities[j][2] + proj_d_vi[2] - proj_d_vj[2] )

               velocities[i] = Vec2f0( velocities[i][1] + proj_d_vj[1] - proj_d_vi[1],
                                       velocities[i][2] + proj_d_vj[2] - proj_d_vi[2] )
               push!(pair_check, (i, j))
            end
        end

       # Check for collisions with the walls and reverse velocity if needed
       if (positions[][i][1] - particle_diam[]/2) < 0 || (positions[][i][1] + particle_diam[]/2) > containerside
          velocities[i] = Vec2f0( -velocities[i][1], velocities[i][2] )
       end
       if (positions[][i][2] - particle_diam[]/2) < 0 || (positions[][i][2] + particle_diam[]/2) > containerside
          velocities[i] = Vec2f0( velocities[i][1], -velocities[i][2] )
       end

       # Update positions based on velocities
       positions[][i] = Point2f0( positions[][i][1] + dt * velocities[i][1], positions[][i][2] + dt * velocities[i][2] )

    end
  
    # Update the drawing
    empty!(pair_check_prev)
    append!(pair_check_prev, pair_check)
end
