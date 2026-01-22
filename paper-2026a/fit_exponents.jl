using GLMakie, JLD2, FourierFlows

@load (@__DIR__) * "/circulation_moments.jld2" hshs orders moments
include("../core/circulation.jl")

grid = TwoDGrid(CPU(), nx=1024, Lx=2π)
fig = Figure()
ax = Axis(fig[1, 1], xscale=log10, yscale=log10)
tb = Textbox(fig[2, 1], placeholder="Enter a exponent",
    validator=Float64, tellwidth=false)


x = Observable(map(hshs) do hsh
    periarea(hsh, 0, grid)
end)

on(tb.stored_string) do str
    α = parse(Float64, str)
    x[] = map(hshs) do hsh
        periarea(hsh, α, grid)
    end
    autolimits!(ax)
end

for y in moments
    scatter!(ax, x, y)
end