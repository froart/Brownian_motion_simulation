using GLMakie, GeometryBasics, Observables, LinearAlgebra

GLMakie.activate!()
GLMakie.closeall()

num_particles = 10 
zoomframeside = Observable(0.4)
containerside = 2 # cm
particle_diam = 0.1 # cm
diam_for_scat = Observable( particle_diam * 240 )
maxspeed      = 0.05

# Define initial positions and velocities for the circles
positions  = Observable([ Point2f0( ( -particle_diam/2 + containerside ) * rand() + particle_diam/2, 
				    ( -particle_diam/2 + containerside ) * rand() + particle_diam/2 ) 
			                 for _ in 1:num_particles ])
velocities = [ Vec2f0( maxspeed * rand() * rand([-1,1]), 
		       maxspeed * rand() * rand([-1,1])) 
                           for _ in 1:num_particles ]

fig = Figure(size = (1240, 860))
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

    pair_check = [] # clear

    for i in 1:num_particles
	# Check for collisions with other particles and exchange velocities if needed
        for j in 1:num_particles
            if j == i #=|| (j, i) in pair_check=# continue end
	    if sqrt(( positions[][i][1] - positions[][j][1] )^2 + ( positions[][i][2] - positions[][j][2] )^2 ) < particle_diam
	       d = [ positions[][j][1] - positions[][i][1], positions[][j][2] - positions[][i][2] ]
	       magn_d_v1 = dot( d, velocities[i] ) / sum( d .^ 2 )
               proj_d_v1 = magn_d_v1 * d
	       if magn_d_v1 > 0
                  velocities[j] = Vec2f0( velocities[j][1] + proj_d_v1[1], velocities[j][2] + proj_d_v1[2])
                  velocities[i] = Vec2f0( velocities[i][1] - proj_d_v1[1], velocities[i][2] - proj_d_v1[2])
               end
	    end
	    #push!(pair_check, (i, j))
	end

        # Check for collisions with the walls and reverse velocity if needed
        if (positions[][i][1] - particle_diam[]/2) < 0 || (positions[][i][1] + particle_diam[]/2) > containerside
           velocities[i] = Vec2f0( velocities[i][1] * (-1), velocities[i][2] )
        end
        if (positions[][i][2] - particle_diam[]/2) < 0 || (positions[][i][2] + particle_diam[]/2) > containerside
           velocities[i] = Vec2f0( velocities[i][1], velocities[i][2] * (-1) )
        end
        #positions[][i] = Point2f0(positions[][i][1], positions[][i][2])
    end

    # Update positions based on velocities
    for i in 1:num_particles
        positions[][i] = Point2f0( positions[][i][1] + velocities[i][1], positions[][i][2] + velocities[i][2] )
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
    println(fig.window.size)
    sleep(1/fps)
    global frame_count = frame_count + 1
end

GLMakie.closeall()
