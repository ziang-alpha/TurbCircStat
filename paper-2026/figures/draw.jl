using CairoMakie, HDF5, LinearAlgebra
using FourierFlows, CUDA
mytheme = merge(
	Theme(
		Axis = (
			xgridvisible = false,
			xminorticksvisible = true,
			xminorticksize = 2,
			xticksize = 3,
			xticklabelsize = 8,
			xlabelsize = 8,
			xlabelpadding = 2,
			xtickalign = 1,
			xminortickalign = 1,
			ygridvisible = false,
			yminorticksvisible = true,
			yminorticksize = 2,
			yticksize = 3,
			yticklabelsize = 8,
			ylabelsize = 8,
			ylabelpadding = 2,
			ytickalign = 1,
			yminortickalign = 1,
		),
		Lines = (
			joinstyle = :round,
			linewidth = 1,
		),
		Scatter = (
			markersize = 5,
		),
		Legend = (
			labelsize = 8,
			padding = 6,
			patchsize = (10, 5),
			colgap = 7,
			patchlabelgap = 3,
			framevisible = false,
		),
		Errorbars = (
			linewidth = 0.5,
			whiskerwidth = 2,
		),
		Text = (
			fontsize = 8,
		),
	),
	theme_latexfonts(),
)

## Global constants
pfile = h5open((@__DIR__) * "/postdata.h5", "r")
cm = 96/2.54
pt = 4/3
Reynolds = [3.3; 3.5; 3.7; 3.9; 4.0; 4.02; 4.04; 4.042; 4.05]
αβs = map(Reynolds) do Re
	data = pfile["wavenumber-energy_spectra"]["$(Re)"]
	x, y = read(data["x"]), read(data["y"])
	αβ = hcat(y[1:300], x[1:300] .^ 2 .* y[1:300])\(2π .* x[1:300])
	Re => (αβ[1], αβ[2])
end
αβs = Dict(αβs)

colormaps = cgrad(:roma, length(Reynolds); categorical = true)
eqSpec(x, Re) = @. 2π * x / (αβs[Re][1] + αβs[Re][2] * x^2)

## Figure 1
with_theme(mytheme) do
	fig = Figure(size = (8.5cm, 8.2cm))

	# Layout
	begin
		ga = GridLayout(fig[1, 1])
		axa = Axis(ga[1, 1];
			xlabel = L"k",
			ylabel = L"E(k)",
			xscale = log10,
			yscale = log10,
			limits = (1, 512, 1e-6, 1),
			xticks = 2 .^ (0:2:8),
			xminorticks = IntervalsBetween(3),
			yminorticks = IntervalsBetween(9),
		)

		gb = GridLayout(fig[1, 2])
		axb = Axis(gb[1, 1];
			xlabel = L"u/\sqrt{⟨u^2⟩}",
			ylabel = "PDF",
			yscale = log10,
			limits = (-7, 7, 1e-10, 1),
			yminorticksvisible = false,
		)

		gc = GridLayout(fig[2, 1])
		axct = Axis(gc[1, 1];
			ylabel = L"\Pi_E/ε",
			xscale = log10,
			limits = (1, 512, -0.6, 0.4),
			xticklabelsvisible = false,
			xticks = 2 .^ (0:2:8),
			xminorticks = IntervalsBetween(3),
		)
		axcb = Axis(gc[2, 1];
			xlabel = L"k",
			ylabel = L"\Pi_\Omega/η",
			xscale = log10,
			limits = (1, 512, -0.3, 0.6),
			xticks = 2 .^ (0:2:8),
			xminorticks = IntervalsBetween(3),
		)
		rowgap!(gc, 8)

		gd = GridLayout(fig[2, 2])
		axd = Axis(gd[1, 1];
			xlabel = L"Re",
			ylabel = L"l_\text{eq}k_f",
			yscale = log10,
			limits = (3.2, 4.1, 1, 10^3),
			yminorticksvisible = false,
		)

		for (gl, lb) in zip([ga, gb, gc, gd], ["(a)", "(b)", "(c)", "(d)"])
			Label(gl[1, 1, TopLeft()], lb;
				fontsize = 10,
				font = :bold,
				padding = (0, 22, -10, 0),
			)
		end

		rowgap!(fig.layout, 10)
		colgap!(fig.layout, 10)
	end

	for (Re, c) in zip(Reynolds, colormaps)
		data = pfile["wavenumber-energy_spectra"]["$(Re)"]
		α, β = αβs[Re]
		x, y = read(data["x"]), read(data["y"])
		lines!(axa, x, y; color = c)
		if αβs[Re][1] > 0 && αβs[Re][2] > 0
			xref = x[4:300]
			lines!(axa, xref, eqSpec(xref, Re); linestyle = :dashdot, color = :black)
		end
	end
	lines!(axa, 2:64, 2e-6 .* (2:64); linestyle = :dash, color = :red)
	lines!(axa, 2:64, 1e-1 .* (2:64) .^ (-1); linestyle = :dash, color = :blue)
	text!(axa, 16, 1e-5; text = L"k", color = :red)
	text!(axa, 16, 1e-2; text = L"k^{-1}", color = :blue)


	for (Re, c) in zip(Reynolds, colormaps)
		data = pfile["velocity-pdf:filtered"]["$(Re)"]
		x, y, Δy = read(data["x"]), read(data["y"]), read(data["Δy"])
		std = sqrt(sum(@. y*x^2) * (x[2]-x[1]))
		scatter!(axb, x ./ std, y .* std; color = c)
		errorbars!(axb, x ./ std, y .* std, Δy .* std; color = c)
	end
	xref = -7:0.05:7
	yref = @. exp(-xref^2/2)/sqrt(2π)
	lines!(axb, xref, yref; linestyle = :dashdot, color = :black)

	for (Re, c) in zip(Reynolds, colormaps)
		efdata, wfdata = pfile["wavenumber-energy_fluxes"]["$(Re)"], pfile["wavenumber-enstrophy_fluxes"]["$(Re)"]
		x1, y1 = read(efdata["x"]), read(efdata["y"])
		x2, y2 = read(wfdata["x"]), read(wfdata["y"])
		lines!(axct, x1, y1; color = c)
		lines!(axcb, x2, y2; color = c)
	end

	x = [Re for Re in Reynolds if αβs[Re][1] > 0 && αβs[Re][2] > 0]
	y = [sqrt(αβs[Re][2]/αβs[Re][1])*1024 / 3sqrt(Re) for Re in Reynolds if αβs[Re][1] > 0 && αβs[Re][2] > 0]
	lines!(axd, x, y; color = :black, alpha = 0.2)
	for (Re, c) in zip(Reynolds, colormaps)
		if αβs[Re][1] > 0 && αβs[Re][2] > 0
			scatter!(axd, Re, sqrt(αβs[Re][2]/αβs[Re][1])*1024 / 3sqrt(Re); color = c)
		end
	end
	Legend(gd[1, 1],
		[LineElement(color = c, linestyle = nothing) for c in colormaps],
		["$(Re)" for Re in Reynolds];
		height = Relative(0.55),
		width = Relative(0.75),
		halign = 0,
		valign = 1,
		nbanks = 2,
	)
	save((@__DIR__) * "/diag.pdf", fig)
end

## Figure 2
with_theme(mytheme) do
	fig = Figure(size = (8.5cm, 8.2cm))

	# Layout
	begin
		ga = GridLayout(fig[1, 1])
		axa = Axis(ga[1, 1];
			xlabel = L"\sqrt{l_1l_2}/2\pi l_\text{eq}",
			ylabel = L"\sqrt{⟨Γ_C^2⟩/⟨Γ_□^2⟩}",
			xscale = log10,
			limits = (1e-2, 1e2, 0.95, 1.2),
			xminorticks = IntervalsBetween(9),
		)

		gb = GridLayout(fig[1, 2])
		axb = Axis(gb[1, 1];
			xlabel = L"l_1/l_2",
			ylabel = L"(⟨Γ_C^p⟩/⟨Γ_□^p⟩)^{1/p}",
			xscale = log10,
			xticks = 10.0 .^ (-2:0),
			xminorticks = IntervalsBetween(9),
			limits = (1e-2, 1, 0.98, 1.02),
		)

		gc = GridLayout(fig[2, 1])
		axc = Axis(gc[1, 1];
			xlabel = L"(l_1 + l_2)/4\pi l_\text{eq}",
			ylabel = L"\sqrt{⟨Γ_C^2⟩/⟨Γ_□^2⟩}",
			xscale = log10,
			limits = (1e-2, 1e2, 0.75, 1.05),
			xminorticks = IntervalsBetween(9),
		)

		gd = GridLayout(fig[2, 2])
		axd = Axis(gd[1, 1];
			xlabel = L"l_1/l_2",
			ylabel = L"(⟨Γ_C^p⟩/⟨Γ_□^p⟩)^{1/p}",
			xscale = log10,
			xticks = 10.0 .^ (-2:0),
			xminorticks = IntervalsBetween(9),
			limits = (1e-2, 1, 0.98, 1.02),
		)

		for (gl, lb) in zip([ga, gb, gc, gd], ["(a)", "(b)", "(c)", "(d)"])
			Label(gl[1, 1, TopLeft()], lb;
				fontsize = 10,
				font = :bold,
				padding = (0, 22, -10, 0),
			)
		end

		rowgap!(fig.layout, 10)
		colgap!(fig.layout, 10)
	end

	for (Re, c) in zip(Reynolds, colormaps)
		data = pfile["loop_sizes-square_moments"]["$(Re)"]
		x_rect = read(data["2"]["5×20"]["x"])
		y_rect = read(data["2"]["5×20"]["y"])
		Δy_rect = read(data["2"]["5×20"]["Δy"])
		npoints = length(x_rect)
		x_sq = read(data["2"]["10×10"]["x"])[1:npoints]
		y_sq = read(data["2"]["10×10"]["y"])[1:npoints]
		Δy_sq = read(data["2"]["10×10"]["Δy"])[1:npoints]

		if αβs[Re][1] > 0 && αβs[Re][2] > 0
			leq = sqrt(αβs[Re][2]/αβs[Re][1])
			x, y = 10x_rect ./ 1024leq, sqrt.(y_rect ./ y_sq)
			lines!(axa, x, y; color = c)
			# scatter!(axa, x, y; color = c)
		end
	end

	for (Re, c) in zip(Reynolds, colormaps)
		data = pfile["loop_sizes-square_moments"]["$(Re)"]
		x_rect = read(data["2"]["4×16"]["x"])
		y_rect = read(data["2"]["4×16"]["y"])
		Δy_rect = read(data["2"]["4×16"]["Δy"])
		npoints = length(x_rect)
		x_sq = read(data["2"]["10×10"]["x"])[1:npoints]
		y_sq = read(data["2"]["10×10"]["y"])[1:npoints]
		Δy_sq = read(data["2"]["10×10"]["Δy"])[1:npoints]

		if αβs[Re][1] > 0 && αβs[Re][2] > 0
			leq = sqrt(αβs[Re][2]/αβs[Re][1])
			x, y = 10x_rect ./ 1024leq, sqrt.(y_rect ./ y_sq)
			lines!(axc, x, y; color = c)
			# scatter!(axc, x, y; color = c)
		end
	end
	markers = [:circle, :utriangle, :dtriangle, :rect, :diamond]
	orders = 2:2:10

	for (n, mk) in zip(orders, markers)
		data = pfile["aspect_ratio-moments"]["area8100"]["4.04"]["$(n)"]
		x = read(data["x"])
		y = read(data["y"])
		Δy = read(data["Δy"])
		scatter!(axb, x, (y ./ y[end]) .^ (1/n); marker = mk, color = :transparent, strokecolor = :black, strokewidth = 0.5)
	end

	for (n, mk) in zip(orders, markers)
		data = pfile["aspect_ratio-moments"]["perimeter180"]["3.5"]["$(n)"]
		x = read(data["x"])
		y = read(data["y"])
		Δy = read(data["Δy"])
		scatter!(axd, x, (y ./ y[end]) .^ (1/n); marker = mk, color = :transparent, strokecolor = :black, strokewidth = 0.5)
	end


	save((@__DIR__) * "/areaperi.pdf", fig)
end


## Figure3
with_theme(mytheme) do
	fig = Figure(size = (8.5cm, 4.1cm))

	# Layout
	begin
		ga = GridLayout(fig[1, 1])
		axa = Axis(ga[1, 1];
			xlabel = L"p",
			ylabel = L"\zeta_p",
			limits = (0, 11, 0, 11),
			xticks = 1:2:11,
			yticks = 1:2:11,
		)

		gb = GridLayout(fig[1, 2])
		axb = Axis(gb[1, 1];
			xlabel = L"l/2\pi l_\text{eq}",
			ylabel = L"\alpha^{p/2} ⟨|Γ_C|^p⟩",
			xscale = log10,
			yscale = log10,
			xticks = 10.0 .^ (-2:0),
			xminorticks = IntervalsBetween(9),
			limits = (1e-2, 1, 0.98, 1.02),
		)

		for (gl, lb) in zip([ga, gb], ["(a)", "(b)"])
			Label(gl[1, 1, TopLeft()], lb;
				fontsize = 10,
				font = :bold,
				padding = (0, 22, -10, 0),
			)
		end

		rowgap!(fig.layout, 10)
		colgap!(fig.layout, 10)
	end

	for (Re, c) in zip(Reynolds, colormaps)
		if αβs[Re][2] > 0 && αβs[Re][1] > 0
			leq = sqrt(αβs[Re][2]/αβs[Re][1])
			α = αβs[Re][1]
			for n in 1:10
				data = pfile["loop_sizes-square_moments"]["$(Re)"]["$(n)"]["10×10"]
				x, y = read(data["x"]), read(data["y"])
				lines!(axb, 10x ./ 1024leq, y .* α ^ (n/2); color = c)
				scatter!(axb, 10x ./ 1024leq, y .* α ^ (n/2); color = :white, strokecolor = c, strokewidth = 0.3)
			end
		end
	end

	sc_area = map(1:10) do n
		data = pfile["loop_sizes-square_moments"]["4.04"]["$(n)"]["10×10"]
		leq = sqrt(αβs[4.04][2]/αβs[4.04][1])
		α = αβs[4.04][1]
		x, y = read(data["x"])[2:16], read(data["y"])[2:16]
		fit = hcat(log.(x), ones(length(x)))\log.(y)
		yref = exp.(hcat(log.(x), ones(length(x))) * fit)
		lines!(axb, 10x ./ 1024leq, yref .* α ^ (n/2); color = :blue)
		fit[1]
	end

	sc_peri = map(1:10) do n
		data = pfile["loop_sizes-square_moments"]["3.5"]["$(n)"]["10×10"]
		leq = sqrt(αβs[3.5][2]/αβs[3.5][1])
		α = αβs[3.5][1]
		x, y = read(data["x"])[2:100], read(data["y"])[2:100]
		fit = hcat(log.(x), ones(length(x)))\log.(y)
		yref = exp.(hcat(log.(x), ones(length(x))) * fit)
		lines!(axb, 10x ./ 1024leq, yref .* α ^ (n/2); color = :red)
		fit[1]
	end
	scatter!(axa, 1:10,sc_area,color = :blue)
	scatter!(axa, 1:10,sc_peri,color = :red)
	autolimits!(axb)


	save((@__DIR__) * "/scaling.pdf", fig)
end

