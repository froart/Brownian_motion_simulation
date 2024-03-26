using GLMakie, GeometryBasics, Observables, LinearAlgebra

GLMakie.activate!()
GLMakie.closeall()

num_particles = 2
zoomframeside = Observable(0.4)
containerside = 2 # cm
particle_diam = 0.1 # cm
diam_for_scat = Observable( particle_diam * 270 )
maxspeed      = 0.05
dt            = 1

# Define initial positions and velocities for the circles
#positions  = Observable([ Point2f0( ( -particle_diam/2 + containerside ) * rand() + particle_diam/2, 
#				    ( -particle_diam/2 + containerside ) * rand() + particle_diam/2 ) 
#			                 for _ in 1:num_particles ])
#velocities = [ Vec2f0( maxspeed * rand() * rand([-1,1]), 
#		       maxspeed * rand() * rand([-1,1])) 
#                           for _ in 1:num_particles ]
positions = Observable( [ Point2f0( 0.5, 1.0 ), Point2f0( 1.5, 1.0 ) ] )
velocities = [ Vec2f0( maxspeed, 0.0 ), Vec2f0( -maxspeed, 0.0 ) ] 

fig = Figure(size = (600, 600))
colsize!(fig.layout, 1, Relative(1))

fig[1, 1] = grid1 = GridLayout()

ax  = Axis(grid1[1, 1], aspect = 1, limits = (0, containerside, 0, containerside))
#hidexdecorations!(ax)
#hideydecorations!(ax)

fig[2, 1] = grid2 = GridLayout()

# Exit button
button = grid2[1, 2] = Button(fig, label = "Exit")
exit   = Observable(false)
on(button.clicks) do n
   exit[] = true
   notify(exit)
end

# Velocity sum label
label = Label(grid2[1, 1], "")

ax2 = Axis(grid1[1,2], aspect = 1, limits = (0, num_particles+1, 0, maxspeed))
xs  = 1:num_particles

# Coordinates for the border
xmin, xmax, ymin, ymax = 0, containerside, 
                         0, containerside
# Draw lines to form a border
lines!(ax, [xmin, xmax], [ymin, ymin], color=:black, linewidth = 8 ) # Bottom
lines!(ax, [xmin, xmax], [ymax, ymax], color=:black, linewidth = 8 ) # Top
lines!(ax, [xmin, xmin], [ymin, ymax], color=:black, linewidth = 8 ) # Left
lines!(ax, [xmax, xmax], [ymin, ymax], color=:black, linewidth = 8 ) # Right

# Draw particles
scat = scatter!(ax, positions, markersize = diam_for_scat)

# Animation function
function update(positions, velocities)

    for i in 1:num_particles
	# Check for collisions with other particles and exchange velocities if needed
        for j in 1:num_particles
	    if j == i continue end
	    # if collision
	    if sqrt( ( positions[][j][1] - positions[][i][1] )^2 + ( positions[][j][2] - positions[][i][2] )^2 ) <= particle_diam
               # unit vector from centre to centre
	       dij      = Vec2f0( positions[][j][1] - positions[][i][1], positions[][j][2] - positions[][i][2] )
	       dij_unit = dij / sqrt( dij[1]^2 + dij[2]^2 )
	       # projection of vi unto dij_unit
	       proj_d_vi_magn = dot(velocities[i], dij_unit) / ( dij_unit[1]^2 + dij_unit[2]^2 ) # formally speaking
               proj_d_vi      = proj_d_vi_magn * dij_unit
	       if proj_d_vi_magn > 0
		  ###   ELASTIC COLLISION
		  # FIXME make momemtum exchange in the same cycle
                  velocities[j] = Vec2f0( velocities[j][1] + proj_d_vi[1], velocities[j][2] + proj_d_vi[2])
                  velocities[i] = Vec2f0( velocities[i][1] - proj_d_vi[1], velocities[i][2] - proj_d_vi[2])
	       end
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
    scat[1][] = positions[]
end

# Run the animation
fps = 60
frame_count = 0
while !exit[]
    display(fig)
    update(positions, velocities)
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
    #println(fig.window.size)
    sleep(1/fps)
    global frame_count = frame_count + 1
end

GLMakie.closeall()
