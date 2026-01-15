using GLMakie, JLD2, FourierFlows

@load (@__DIR__) * "/circulation_moments.jld2" hshs orders moments
include("../core/circulation.jl")

grid = TwoDGrid(GPU(), nx = 1024, Lx = 2π)
fig = Figure()
ax = Axis(fig[1, 1])
tb = Textbox(fig[2, 1], placeholder = "Enter a exponent",
	validator = Float64, tellwidth = false)

tb2 = Textbox(fig[3, 1], placeholder = "Enter a integer",
	validator = Int, tellwidth = false)


x = Observable(map(hshs) do hsh
	periarea(hsh, 0, grid)
end)

y = Observable(moments[1])
on(tb.stored_string) do str
	α = parse(Float64, str)
	x[] = map(hshs) do hsh
		periarea(hsh, α, grid)
	end
    autolimits!(ax)
end
on(tb2.stored_string) do str
	n = parse(Int, str)
	y[] = moments[n]
    autolimits!(ax)
end
scatter!(ax, x, y)
