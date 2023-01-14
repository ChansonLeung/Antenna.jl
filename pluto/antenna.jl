### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 8c4ceaf6-4a7c-4c35-b0d3-46e44e952177
begin
	using Pkg
	Pkg.activate("v1.8", shared=true)
	using Revise
	using Plots
	using PlutoUI
	using Antenna
	plotlyjs()
end

# ╔═╡ 714dae6e-8c3a-11ed-31ce-85ea53b0c14e
begin
	md"""
	- Nx $(@bind Nx Slider(1:10, default=8, show_value=true))
	- Ny $(@bind Ny Slider(1:10, default=8, show_value=true))
	- θₜ $(@bind θₜ Slider(0:0.1:pi, show_value=true) ) 
	- ϕₜ $(@bind ϕₜ Slider(0:0.1:2pi, show_value=true))
	"""
end

# ╔═╡ bcd8f860-ea2c-4cc6-90d6-286674940753
begin
	set_param(freq=2.7e9)
	point = point_rectangle(Nx = Nx,Ny=Ny, dx=λ/2, dy=λ/2,pattern=pattern_identity)
	
	res_pattern = cal_pattern(point, θₜ ,ϕₜ)
	plot(res_pattern, zlims=(-40,30))
end

# ╔═╡ Cell order:
# ╠═8c4ceaf6-4a7c-4c35-b0d3-46e44e952177
# ╠═714dae6e-8c3a-11ed-31ce-85ea53b0c14e
# ╠═bcd8f860-ea2c-4cc6-90d6-286674940753
