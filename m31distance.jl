### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 5a96edd0-a432-11eb-119b-bb814232c52a
begin
	using PlutoUI
	using FITSIO
	using DataFrames
	using KernelDensity
	using Makie
	using Unitful
	using UnitfulAstro
	using Polynomials
	using Statistics
end

# ╔═╡ 112f6c01-5ef7-4980-81ce-282b242a9cf2
# Change the width of the notebook
html"""<style>
main {
    max-width: 1200px;
}
"""

# ╔═╡ 5e95a0fd-9f0b-40a8-9cc8-f945160d28e8
md"""
# Andmetöötlus

### Andmete sisselugemine.

Gaias olevate tsefeiidide andmed on saadud tabelsit gaiadr2.vari_cepheid. See on liidetud gaia EDR3 ja Pan-STARRS-1 andmetega, et saada vastavalt parallaksid ja nähtavatad heledused õigetes filtrites.
"""

# ╔═╡ 26956c96-480e-4e7b-926c-3202b7bce3df
df_gaia_raw = DataFrame(FITS("gaia_cepheid_panstarrs2.fits")[2])

# ╔═╡ 0d08c76d-48f1-4ccf-813b-cf59e12273b1


# ╔═╡ b0e70be5-17c8-4d84-a23c-84e0526d9e7b
df_m31_raw = DataFrame(FITS("M31_cepheid.fits")[2])

# ╔═╡ f557afac-405f-4ef2-85de-15bd2ab61297
md"Sisendväärtus: $(@bind poe PlutoUI.Slider(1:50, show_value=true, default=10))"

# ╔═╡ 41c3b585-032b-4d50-bf02-9a4aa5034062
combine(groupby(df_gaia_raw, ["type_best_classification"]), nrow)

# ╔═╡ 9cced706-0fcb-4d28-ad3b-78f2b16c79b6
begin
	parallax_over_error_lim = 5
	mask_gaia_wrangling = (df_gaia_raw.mode_best_classification .== "FUNDAMENTAL") .& 
		.!isnan.(df_gaia_raw.pf) .&
		.!isnan.(df_gaia_raw.parallax) .&
		(0 .< df_gaia_raw.parallax .< 4) .&
		(df_gaia_raw.parallax_over_error .>= parallax_over_error_lim)
	df_gaia_wrangled = df_gaia_raw[mask_gaia_wrangling, setdiff(names(df_gaia_raw), ("source_id", "type_best_classification"))]
end

# ╔═╡ b2341b40-71d4-4bee-afe6-41d2cd7baac3
filter(x -> x.ID ∈ df_gaia_wrangled.obj_id, df_m31_raw)

# ╔═╡ f8002ac9-476d-470e-9a72-3d0656f58412
md"""
### Andmete puhastamine

Eemaldan Gaia andmetest read, kus on kas periood või parallaks ei ole määratud. Samuti kasutan ainult tulemusi, mille prallaksi viga on väiksem kui $(round(Int, 100/parallax_over_error_lim))% ja kasutan ainult klassikalisi (delta) tsefeiide. Tekitan ka eraldi maskid iga värvifiltri jaoks (g,r,i,z).

"""

# ╔═╡ 11e6f16c-7789-4e3d-bada-8aa39f02924f
mean(filter(!isnan, filter(row -> row.mode_best_classification == "FUNDAMENTAL", df_gaia_wrangled).g_mean_psf_mag).g_mean_psf_mag)

# ╔═╡ a39559cf-1dc4-4409-9572-f33f73dbaa32
begin
	maskg = .!isnan.(df_gaia_wrangled.g_mean_psf_mag)
	maskr = .!isnan.(df_gaia_wrangled.r_mean_psf_mag)
	maski = .!isnan.(df_gaia_wrangled.i_mean_psf_mag)
	maskz = .!isnan.(df_gaia_wrangled.z_mean_psf_mag)
end;

# ╔═╡ 926dc1d8-5eaf-4a8b-8f31-a40a5180941b
combine(groupby(df_m31_raw, ["Type"]), nrow)

# ╔═╡ 29e0bf05-def9-4373-8ee1-47116faba7fb
begin
	mask_m31_wrangling = df_m31_raw.Type .∈ Ref(("FM", "FO"))
	df_m31_wrangled = df_m31_raw[mask_m31_wrangling, :]
	# for colname in names(df_m31_wrangled)
	# 	replace!(df_m31_wrangled[!, colname], NaN => missing)
	# end
	# df_m31_wrangled
end

# ╔═╡ 6aabe140-4745-4b34-9733-c723340ad683


# ╔═╡ 49e95cb1-7faa-454d-8c02-1d486222eef0
md"""
# Analüüs
"""

# ╔═╡ 42ece2d8-6cba-47ae-bdb3-bc52707f858d
begin
	dist_from_parallax_mas(parallax) = u"kpc"/parallax
	absmag(appmag::Real, parallax_in_mas::Real) = appmag - 5. * log10(dist_from_parallax_mas(parallax_in_mas)/u"pc") + 5
	distmod(appmag::Real, absmag::Real) = appmag - absmag
	distance(distmod::Real) = 10^(distmod/5 + 1)
end

# ╔═╡ 8e1c5586-11b5-4194-8df1-3844e875e4bb
absmag(12, 1)

# ╔═╡ cf16eb08-69e2-4570-80cd-26d49d93f4ac
row = df_gaia_wrangled[1,:]

# ╔═╡ 10523041-ad3d-4a66-b647-f7c6a2746430
row.g_mean_psf_mag

# ╔═╡ cba03e79-eaa4-430b-9855-f348a36b87db
absmag(row.g_mean_psf_mag, row.parallax)

# ╔═╡ cce43589-c8c2-44c4-aee5-a783053c51cb
begin
	absmagg = absmag.(df_gaia_wrangled[maskg, :g_mean_psf_mag], df_gaia_wrangled[maskg, :parallax])
	logpg = log10.(df_gaia_wrangled[maskg, :pf])
	
	absmagr = absmag.(df_gaia_wrangled[maskr, :r_mean_psf_mag], df_gaia_wrangled[maskr, :parallax])
	logpr = log10.(df_gaia_wrangled[maskr, :pf])

	absmagi = absmag.(df_gaia_wrangled[maski, :i_mean_psf_mag], df_gaia_wrangled[maski, :parallax])
	logpi = log10.(df_gaia_wrangled[maski, :pf])
	
	absmagz = absmag.(df_gaia_wrangled[maskz, :z_mean_psf_mag], df_gaia_wrangled[maskz, :parallax])
	logpz = log10.(df_gaia_wrangled[maskz, :pf])
end;

# ╔═╡ d94e9303-2955-420e-afbe-5451c224f018
begin
	maskicut = absmagi .< 0;
	plumicut = fit(logpi[maskicut], absmagi[maskicut], 1)
end

# ╔═╡ e3643abf-d0c5-41e7-a1a5-2b97cd93a4e7
begin
	plumg = fit(logpg, absmagg, 1)
	plumr = fit(logpr, absmagr, 1)
	plumi = fit(logpi, absmagi, 1)
	plumz = fit(logpz, absmagz, 1)
end;

# ╔═╡ ce737155-b540-4513-8aef-538916a349cf
begin
	ff = Figure(resolution = (1200,600))
	a1 = ff[1,1] = Axis(ff, xlabel = "Absoluutne heledus")
	a2 = ff[1,2] = Axis(ff, xlabel = "Näiline heledus")
	a3 = ff[1,3] = Axis(ff, xlabel = "Parallaks")
	hist!(a1, absmagi, bins=20, color=:lightblue)
	hist!(a2, df_gaia_wrangled[maski, :i_mean_psf_mag], bins=20, color=:red)
	hist!(a3, df_gaia_wrangled[maski, :parallax], bins=20, color=:orange)
	ff
end

# ╔═╡ e29bf196-1701-4180-b4ff-ae6c25d6b9b3
begin

	f = Figure(resolution = (1200, 800))
	ax = Axis(f[1,1], xlabel = "log (P(d))", ylabel = "g-band magnitude")
	ax2 = Axis(f[1,1], yaxisposition = :right, xaxisposition = :top)
	ax.xminorticksvisible = true
	ax.yminorticksvisible = true
	ax.xgridvisible = false
	ax.xminorticks = IntervalsBetween(4)
	ax.xminortickalign = 1
	ax.xtickalign = 1
	ax.xminorticksize = 10
	ax.xticksize = 15
	ax.yreversed = true
	ax2.yreversed = true

	scatter!(ax, logpi, absmagi, color=:lightblue)
	scatter!(ax2, logpi, absmagi, markersize = 0)
	line_ends = [minimum(logpi), maximum(logpi)]
	lines!(ax, line_ends, plumi.(line_ends), color=:red, linewidth=2)
	
	f
	
end

# ╔═╡ e53b9d50-0bcb-4932-89cf-2d848d938bf6
begin

	fff = Figure(resolution = (1200, 800))
	ax1 = Axis(fff[1,1], xlabel = "log (P(d))", ylabel = "g-band magnitude")
	ax1.xminorticksvisible = true
	ax1.yminorticksvisible = true
	ax1.xgridvisible = false
	ax1.xminorticks = IntervalsBetween(4)
	ax1.xminortickalign = 1
	ax1.xtickalign = 1
	ax1.xminorticksize = 10
	ax1.xticksize = 15
	ax1.yreversed = true

	scatter!(ax1, logpi[maskicut], absmagi[maskicut], color=:lightblue)
	lines!(ax1, line_ends, plumicut.(line_ends), color=:red, linewidth=2)
	
	fff
	
end

# ╔═╡ 09212331-cd90-4b49-bcab-adeb6d70d0bc
length(df_m31_raw.rP1mag) - count(isequal(-99), df_m31_raw.rP1mag)

# ╔═╡ 0f792688-b1e5-4aa2-9150-465852abf72e
length(df_m31_raw.Type) - count(isequal(" "), df_m31_raw.Type)

# ╔═╡ 6e9b4c42-e521-4974-8ae6-908682b6bf2f
count(isequal("FM"), df_m31_raw.Type)

# ╔═╡ 8ffe0786-863a-4fbe-befa-cbb89c4e547c
df_m31_raw.ID

# ╔═╡ 21d68135-f5b5-4444-b351-ff5b0015dded
mean(filter(!isnan, plumg.(df_m31_wrangled.gP1mag)))

# ╔═╡ Cell order:
# ╠═5a96edd0-a432-11eb-119b-bb814232c52a
# ╟─112f6c01-5ef7-4980-81ce-282b242a9cf2
# ╟─5e95a0fd-9f0b-40a8-9cc8-f945160d28e8
# ╠═26956c96-480e-4e7b-926c-3202b7bce3df
# ╠═0d08c76d-48f1-4ccf-813b-cf59e12273b1
# ╠═b0e70be5-17c8-4d84-a23c-84e0526d9e7b
# ╠═b2341b40-71d4-4bee-afe6-41d2cd7baac3
# ╠═f8002ac9-476d-470e-9a72-3d0656f58412
# ╟─f557afac-405f-4ef2-85de-15bd2ab61297
# ╠═41c3b585-032b-4d50-bf02-9a4aa5034062
# ╠═9cced706-0fcb-4d28-ad3b-78f2b16c79b6
# ╠═11e6f16c-7789-4e3d-bada-8aa39f02924f
# ╠═a39559cf-1dc4-4409-9572-f33f73dbaa32
# ╠═926dc1d8-5eaf-4a8b-8f31-a40a5180941b
# ╠═29e0bf05-def9-4373-8ee1-47116faba7fb
# ╠═6aabe140-4745-4b34-9733-c723340ad683
# ╠═49e95cb1-7faa-454d-8c02-1d486222eef0
# ╠═42ece2d8-6cba-47ae-bdb3-bc52707f858d
# ╠═8e1c5586-11b5-4194-8df1-3844e875e4bb
# ╠═cf16eb08-69e2-4570-80cd-26d49d93f4ac
# ╠═10523041-ad3d-4a66-b647-f7c6a2746430
# ╠═cba03e79-eaa4-430b-9855-f348a36b87db
# ╠═cce43589-c8c2-44c4-aee5-a783053c51cb
# ╠═d94e9303-2955-420e-afbe-5451c224f018
# ╠═e3643abf-d0c5-41e7-a1a5-2b97cd93a4e7
# ╠═ce737155-b540-4513-8aef-538916a349cf
# ╠═e29bf196-1701-4180-b4ff-ae6c25d6b9b3
# ╠═e53b9d50-0bcb-4932-89cf-2d848d938bf6
# ╠═09212331-cd90-4b49-bcab-adeb6d70d0bc
# ╠═0f792688-b1e5-4aa2-9150-465852abf72e
# ╠═6e9b4c42-e521-4974-8ae6-908682b6bf2f
# ╠═8ffe0786-863a-4fbe-befa-cbb89c4e547c
# ╠═21d68135-f5b5-4444-b351-ff5b0015dded
