using GLMakie, GeometryBasics, Observables

GLMakie.activate!()
GLMakie.closeall()

num_particles = 50
particle_diam = Observable(20)
zoomframeside = Observable(0.4)
containerside = 2
maxspeed      = 0.1

# Define initial positions and velocities for the circles
positions  = Observable([ Point2f0(containerside * rand(), 
				   containerside * rand()) 
			               for _ in 1:num_particles ])
velocities = [ Vec2f0( maxspeed * rand() * rand([-1,1]), 
		       maxspeed * rand() * rand([-1,1])) 
                           for _ in 1:num_particles ]

fig = Figure()
ax  = Axis(fig[1, 1], aspect = 1, limits = (0, containerside, 0, containerside))
#hidexdecorations!(ax)
#hideydecorations!(ax)

# Coordinates for the border
xmin, xmax, ymin, ymax = 0, containerside, 
                         0, containerside

# Draw lines to form a border
lines!(ax, [xmin, xmax], [ymin, ymin], color=:black, linewidth = 8 ) # Bottom
lines!(ax, [xmin, xmax], [ymax, ymax], color=:black, linewidth = 8 ) # Top
lines!(ax, [xmin, xmin], [ymin, ymax], color=:black, linewidth = 8 ) # Left
lines!(ax, [xmax, xmax], [ymin, ymax], color=:black, linewidth = 8 ) # Right

scat = scatter!(ax, positions, markersize = particle_diam)

# Animation function
function update(positions, velocities)
    for i in 1:num_particles
        # Update positions based on velocities
        new_x = positions[][i][1] + velocities[i][1]
        new_y = positions[][i][2] + velocities[i][2]

        # Check for collisions with the walls and reverse velocity if needed
	# FIXME conditions don't work correctly
	if (new_x - particle_diam[]/2) < 0 || (new_x + particle_diam[]/2) > containerside 
           velocities[i] = Vec2f0( velocities[i][1] * (-1), velocities[i][2] )
        end
	if (new_y - particle_diam[]/2) < 0 || (new_y + particle_diam[]/2) > containerside
           velocities[i] = Vec2f0( velocities[i][1], velocities[i][2] * (-1) )
        end
        positions[][i] = Point2f0(new_x, new_y)
    end
    scat[1][] = positions[]
end

# Run the animation
fps = 60
frame_count = 0
while true
    display(fig)
    update(positions, velocities)
    notify(positions)
    if frame_count > 121
       posx = positions[][1][1]
       posy = positions[][1][2]
       limits!(ax, posx - zoomframeside[]/2, posx + zoomframeside[]/2, posy - zoomframeside[]/2, posy + zoomframeside[]/2)
       global particle_diam[] = 50
    end
    sleep(1/fps)
    global frame_count = frame_count + 1
end

