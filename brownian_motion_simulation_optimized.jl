using GLMakie, GeometryBasics, Observables, LinearAlgebra, Distributions

# ======================
# Structures
# ======================
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

    function Particles(diam::Float64, num::Int64, boundaries::Container, maxspeed::Float64)
        @assert diam > 0 "Diameter must be positive"
        @assert num > 0 "Number of particles must be positive"

        pos = Observable([Point2f0(rand(Uniform(boundaries.x_min + diam/2, boundaries.x_max - diam/2)),
                                    rand(Uniform(boundaries.y_min + diam/2, boundaries.y_max - diam/2)))
                          for i in 1:num])

        vel = [Vec2f0(rand(Uniform(-maxspeed, maxspeed)),
                      rand(Uniform(-maxspeed, maxspeed))) for i in 1:num]

        collided = falses(num, num + 4)
        new(diam, num, pos, vel, collided)
    end
end

# ======================
# Simulation step
# ======================
function step!(pts::Particles, container::Container, dt::Float64)
    for i in 1:pts.num
        # Particle-particle collisions
        for j in (i+1):pts.num
            dist = sqrt((pts.pos[][j][1] - pts.pos[][i][1])^2 + (pts.pos[][j][2] - pts.pos[][i][2])^2)

            if dist <= pts.diam && !pts.collided[i,j]
                dij      = Vec2f0(pts.pos[][j][1] - pts.pos[][i][1], pts.pos[][j][2] - pts.pos[][i][2])
                dij_unit = dij / sqrt(dij[1]^2 + dij[2]^2)
                dji_unit = -dij_unit

                proj_d_vi = dij_unit * dot(pts.vel[i], dij_unit) / (dij_unit[1]^2 + dij_unit[2]^2)
                proj_d_vj = dji_unit * dot(pts.vel[j], dji_unit) / (dji_unit[1]^2 + dji_unit[2]^2)

                pts.vel[j] = Vec2f0(pts.vel[j][1] + proj_d_vi[1] - proj_d_vj[1],
                                     pts.vel[j][2] + proj_d_vi[2] - proj_d_vj[2])

                pts.vel[i] = Vec2f0(pts.vel[i][1] + proj_d_vj[1] - proj_d_vi[1],
                                     pts.vel[i][2] + proj_d_vj[2] - proj_d_vi[2])

                pts.collided[i,j] = true
            elseif dist > pts.diam && pts.collided[i,j]
                pts.collided[i,j] = false
            end
        end

        # Wall collisions
        # Left
        if (pts.pos[][i][1] - pts.diam/2) <= container.x_min && !pts.collided[i, pts.num+1]
            pts.vel[i] = Vec2f0(-pts.vel[i][1], pts.vel[i][2])
            pts.collided[i, pts.num+1] = true
        elseif (pts.pos[][i][1] - pts.diam/2) > container.x_min && pts.collided[i, pts.num+1]
            pts.collided[i, pts.num+1] = false
        end

        # Right
        if (pts.pos[][i][1] + pts.diam/2) >= container.x_max && !pts.collided[i, pts.num+2]
            pts.vel[i] = Vec2f0(-pts.vel[i][1], pts.vel[i][2])
            pts.collided[i, pts.num+2] = true
        elseif (pts.pos[][i][1] + pts.diam/2) < container.x_max && pts.collided[i, pts.num+2]
            pts.collided[i, pts.num+2] = false
        end

        # Bottom
        if (pts.pos[][i][2] - pts.diam/2) <= container.y_min && !pts.collided[i, pts.num+3]
            pts.vel[i] = Vec2f0(pts.vel[i][1], -pts.vel[i][2])
            pts.collided[i, pts.num+3] = true
        elseif (pts.pos[][i][2] - pts.diam/2) > container.y_min && pts.collided[i, pts.num+3]
            pts.collided[i, pts.num+3] = false
        end

        # Top
        if (pts.pos[][i][2] + pts.diam/2) >= container.y_max && !pts.collided[i, pts.num+4]
            pts.vel[i] = Vec2f0(pts.vel[i][1], -pts.vel[i][2])
            pts.collided[i, pts.num+4] = true
        elseif (pts.pos[][i][2] + pts.diam/2) < container.y_max && pts.collided[i, pts.num+4]
            pts.collided[i, pts.num+4] = false
        end

        # Update position
        pts.pos[][i] = Point2f0(pts.pos[][i][1] + dt * pts.vel[i][1],
                                pts.pos[][i][2] + dt * pts.vel[i][2])
    end
end

function step_optimized!(pts::Particles, container::Container, dt::Float64)
    # ----------------------
    # Build spatial grid
    # ----------------------
    cell_size = pts.diam
    nx = ceil(Int, (container.x_max - container.x_min) / cell_size)
    ny = ceil(Int, (container.y_max - container.y_min) / cell_size)

    grid = [Int[] for _ in 1:(nx*ny)]
    cell_coords = Dict{Int, Tuple{Int,Int}}()

    for i in 1:pts.num
        x_idx = clamp(floor(Int, (pts.pos[][i][1] - container.x_min) / cell_size) + 1, 1, nx)
        y_idx = clamp(floor(Int, (pts.pos[][i][2] - container.y_min) / cell_size) + 1, 1, ny)
        cell_id = (y_idx - 1) * nx + x_idx
        push!(grid[cell_id], i)
        cell_coords[i] = (x_idx, y_idx)
    end

    # ----------------------
    # Helper: neighboring cells
    # ----------------------
    neighboring_cells(x::Int, y::Int) = begin
        xs = clamp.(x-1:x+1, 1, nx)
        ys = clamp.(y-1:y+1, 1, ny)
        cells = Int[]
        for xi in xs, yi in ys
            push!(cells, (yi-1)*nx + xi)
        end
        return cells
    end

    # ----------------------
    # Particle-particle collisions
    # ----------------------
    for i in 1:pts.num
        xi, yi = cell_coords[i]
        for cell_id in neighboring_cells(xi, yi)
            for j in grid[cell_id]
                if j > i
                    dist = sqrt((pts.pos[][j][1] - pts.pos[][i][1])^2 +
                                (pts.pos[][j][2] - pts.pos[][i][2])^2)

                    if dist <= pts.diam && !pts.collided[i,j]
                        dij      = Vec2f0(pts.pos[][j][1] - pts.pos[][i][1],
                                           pts.pos[][j][2] - pts.pos[][i][2])
                        dij_unit = dij / sqrt(dij[1]^2 + dij[2]^2)
                        dji_unit = -dij_unit

                        proj_d_vi = dij_unit * dot(pts.vel[i], dij_unit) / (dij_unit[1]^2 + dij_unit[2]^2)
                        proj_d_vj = dji_unit * dot(pts.vel[j], dji_unit) / (dji_unit[1]^2 + dji_unit[2]^2)

                        pts.vel[j] = Vec2f0(pts.vel[j][1] + proj_d_vi[1] - proj_d_vj[1],
                                             pts.vel[j][2] + proj_d_vi[2] - proj_d_vj[2])

                        pts.vel[i] = Vec2f0(pts.vel[i][1] + proj_d_vj[1] - proj_d_vi[1],
                                             pts.vel[i][2] + proj_d_vj[2] - proj_d_vi[2])

                        pts.collided[i,j] = true
                    elseif dist > pts.diam && pts.collided[i,j]
                        pts.collided[i,j] = false
                    end
                end
            end
        end
    end

    # ----------------------
    # Wall collisions
    # ----------------------
    for i in 1:pts.num
        # Left wall
        if (pts.pos[][i][1] - pts.diam/2) <= container.x_min && !pts.collided[i, pts.num+1]
            pts.vel[i] = Vec2f0(-pts.vel[i][1], pts.vel[i][2])
            pts.collided[i, pts.num+1] = true
        elseif (pts.pos[][i][1] - pts.diam/2) > container.x_min && pts.collided[i, pts.num+1]
            pts.collided[i, pts.num+1] = false
        end

        # Right wall
        if (pts.pos[][i][1] + pts.diam/2) >= container.x_max && !pts.collided[i, pts.num+2]
            pts.vel[i] = Vec2f0(-pts.vel[i][1], pts.vel[i][2])
            pts.collided[i, pts.num+2] = true
        elseif (pts.pos[][i][1] + pts.diam/2) < container.x_max && pts.collided[i, pts.num+2]
            pts.collided[i, pts.num+2] = false
        end

        # Bottom wall
        if (pts.pos[][i][2] - pts.diam/2) <= container.y_min && !pts.collided[i, pts.num+3]
            pts.vel[i] = Vec2f0(pts.vel[i][1], -pts.vel[i][2])
            pts.collided[i, pts.num+3] = true
        elseif (pts.pos[][i][2] - pts.diam/2) > container.y_min && pts.collided[i, pts.num+3]
            pts.collided[i, pts.num+3] = false
        end

        # Top wall
        if (pts.pos[][i][2] + pts.diam/2) >= container.y_max && !pts.collided[i, pts.num+4]
            pts.vel[i] = Vec2f0(pts.vel[i][1], -pts.vel[i][2])
            pts.collided[i, pts.num+4] = true
        elseif (pts.pos[][i][2] + pts.diam/2) < container.y_max && pts.collided[i, pts.num+4]
            pts.collided[i, pts.num+4] = false
        end

        # Update positions
        pts.pos[][i] = Point2f0(pts.pos[][i][1] + dt * pts.vel[i][1],
                                pts.pos[][i][2] + dt * pts.vel[i][2])
    end
end

# ======================
# Utilities
# ======================
to_pixels(value_in_cm::Float64) = (value_in_cm / 2.54) * sqrt(1920^2 + 1080^2) / 15.6

# ======================
# Main
# ======================
GLMakie.activate!()
GLMakie.closeall()

# Parameters
particles_num    = 500
speed_max        = 1.5
dt               = 0.1
particle_diameter = 0.23
window_width     = 25.0
window_height    = 25.0
container_size   = 10.0

container = Container(0.0, container_size, 0.0, container_size)
pts = Particles(particle_diameter, particles_num, container, speed_max)

# Figure & Axis
fig = Figure(size=(to_pixels(window_width), to_pixels(window_height)))
ax  = Axis(fig[1,1][1,1], aspect=1,
           limits=(container.x_min, container.x_max, container.y_min, container.y_max))

hidexdecorations!(ax)
hideydecorations!(ax)

scat = scatter!(ax, pts.pos, markersize=to_pixels(particle_diameter), markerspace=:pixel)

quit = Observable(false)
on(events(fig).keyboardbutton) do event
    if event.action == Keyboard.press && event.key == Keyboard.escape
        quit[] = true
        notify(quit)
    end
end

fps = 60
display(fig)




# ======================
# Simulation loop
# ======================
#n = 0
#t = 0.0
while !quit[]
    #n += 1
    #t += @elapsed step!(pts, container, dt)
    t += @elapsed step_optimized!(pts, container, dt)
    notify(pts.pos)
    sleep(1/fps)
    #n == 2000 && break
end

#avg = t/n

GLMakie.closeall()
