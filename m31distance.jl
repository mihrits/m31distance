### A Pluto.jl notebook ###
# v0.14.5

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

# ╔═╡ 9e747c0c-ab6d-4a15-aa6d-1d17716f06f3
md"""
# Andromeda (M31) galaktika kauguse määramine tsefeiidide periood-heledus seose põhjal

Käesoleva töö eesmärk on määrata Andromeda kaugus muutlike tähtede, [tsefeiidide](https://en.wikipedia.org/wiki/Cepheid_variable) perioodide ja heleduste põhjal. Tsefeiididel esineb seos heleduse muutuse perioodi ja absoluutse heleduse vahel ([LPR](https://en.wikipedia.org/wiki/Period-luminosity_relation)), mille me esmalt kalibreerime ja seejärel rakendame saadud seost Andromeda galaktikas asuvate tsefeiidide andmetele, et leida galaktika kaugus.
"""

# ╔═╡ 5e95a0fd-9f0b-40a8-9cc8-f945160d28e8
md"""
# Andmetöötlus

## Andmete sisselugemine.

Tsefeiidide perioodi ja heleduse seose määramiseks on tarvis tsefeiidide perioode, näivat heledust ja kaugust. Perioodid ja näivad heledused, koos neeldumisparandiga, saame [Groenewegen (2018)](https://ui.adsabs.harvard.edu/abs/2018A%26A...619A...8G/abstract) (edaspidi G18) artikli andmestikust. Kaugused, parallaksid, saame [Gaia Early DR3](https://www.cosmos.esa.int/web/gaia/early-data-release-3) andmestikust ja tsefeiidide tüübid [Gaia DR2 tsefeiidide](https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_variability_tables/ssec_dm_vari_cepheid.html) andmestikust. Edaspidi nimetame neid andmeid **kalibratsiooniandmestikuks**.

Andromeda tsefeiidide kohta sain andmed [Kodric et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018AJ....156..130K/abstract) artikli andmestikust nimega PAndromeda. Sealt kasutan perioodi ja näilise heleduse väärtusi. Edaspidi nimetame neid andmeid **M31 andmestikuks**.

Andmete allalaadimiseks kasutan [TOPCAT](http://www.star.bris.ac.uk/~mbt/topcat/) tarkvara. Kõigepealt esitasin päringu G18 andmetele Vizieri tabelist "J/A+A/619/A8/table1"
```
SELECT Name, Period, Vmag, "E(B-V)" as reddening, "_RA", "_DE"
FROM "J/A+A/619/A8/table1"
```

Seejärel leidsin igale tähele lähima vaste 1 kaaresekundi kauguselt Gaia EDR3 kataloogist (gaiaedr3.gaia\_source), kasutades selle jaoks TOPCATi kasutajaliidest. Kui vasted olid leitud, siis ühendsin kataloogi Gaia DR2 tsefeiidide kataloogiga (gaiadr2.vari\_cepheid), et saada sealt tsefeiidide tüübid, päringuga
```
SELECT cep.Name, cep2.mode_best_classification as gaia2_mode, cep.Period, cep2.pf,
cep.Vmag, cep.reddening, cep.parallax, cep.parallax_over_error, cep.ruwe,
cep.ra_epoch2000, cep.dec_epoch2000, cep.angDist
FROM TAP_UPLOAD.t2 as cep
JOIN gaiaedr3.dr2_neighbourhood as nbr2 ON nbr2.dr3_source_id = cep.source_id
JOIN gaiadr2.vari_cepheid as cep2 ON nbr2.dr2_source_id = cep2.source_id
```

M31 andmed sain Vizieri tabelis "J/AJ/156/130/main" päringuga
```
SELECT Pr, gP1mag, rP1mag, Type, RAJ2000, DEJ2000
FROM "J/AJ/156/130/main"
```
"""

# ╔═╡ 5ef8ac87-f667-4ecc-8d57-d7dbf71a3e5f
df_calibration_raw = DataFrame(FITS("cal_gaia_g18.fits")[2])

# ╔═╡ b0e70be5-17c8-4d84-a23c-84e0526d9e7b
df_m31_raw = DataFrame(FITS("m31.fits")[2])

# ╔═╡ f8002ac9-476d-470e-9a72-3d0656f58412
md"""
## Andmete puhastamine

Analüüsis kasutame ainult [1. tüüpi põhisagedusega tsefeiide](https://en.wikipedia.org/wiki/Classical_Cepheid_variable). Selleks jaoks tuleb nii kalibratsiooni- kui ka M31 andmestikust välja filtreerida vastavalt "FUNDAMENTAL" ja "FM" tüüpi tähed.

Kalibratsiooniandmestiku tsefeiidide tüübid:
"""

# ╔═╡ 8ccdaf15-5af0-441e-8eee-42c5caa00a9c
combine(groupby(df_calibration_raw, :gaia2_mode), nrow)

# ╔═╡ 70c39032-3255-4f46-97a0-70e53957b6d3
md"M31 andmestiku tsefeiidide tüübid:"

# ╔═╡ 926dc1d8-5eaf-4a8b-8f31-a40a5180941b
combine(groupby(df_m31_raw, ["Type"]), nrow)

# ╔═╡ 0a334697-f3a6-4b28-8f6e-3d0add984e05
md"""
Kalibratsiooniandmestikust kasutame ainult neid sissekandeid, mille parallaksi väärtused on suuremad kui 0. Lisaks paneme veel peale täiendavad piirangud Gaia poolt antud vigadele nende arvutatud andmetele ([RUWE](https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_ruwe.html) < 1,4), mis tuginevad sellel, kui hästi tähemudel sobib mõõtmistulemustega. Lõpuks kasutame ainult neid tähti, mille parallaksi suhteline viga on väiksem kui $(@bind parallax_relative_error NumberField(1:1:50, default=20))%.
"""

# ╔═╡ dc8bfc19-86c7-4d78-a480-cf4590f12dc2
begin
	
	function mask_calibration(df::DataFrame, typecol::String, types::Union{Symbol, Tuple{Vararg{String}}}, ruwe_limit::Float64, parallax_over_error_limit::Real)
		if types == :All
			type_mask = true
		else
			type_mask = df[:, typecol] .∈ Ref(types)
		end
		type_mask .& .!isnan.(df.parallax) .& (df.parallax .> 0) .&
			(df.parallax_over_error .> parallax_over_error_limit) .&
			(df.ruwe .< ruwe_limit) .&
			(df.period .> 10)
	end 
	
	parallax_over_error_lim = 100. / parallax_relative_error
	mask_calibration_wrangling = mask_calibration(df_calibration_raw, "gaia2_mode", ("FUNDAMENTAL",), 1.4, parallax_over_error_lim)
	
	df_calibration = df_calibration_raw[mask_calibration_wrangling, [:name, :period, :vmag, :reddening, :parallax, :ra_epoch2000, :dec_epoch2000]]
end ;

# ╔═╡ 22d57db1-8374-41f0-bd7a-ed1aa4eb4889
md"Sellisel juhul on kalibratsiooniandmestikus $(size(df_calibration, 1)) kirjet."

# ╔═╡ c6d59e16-70f8-47b3-93c8-75247f6dad6f
begin

	fig_cal_pos = Figure(resolution = (1200, 800))
	ax_cal_pos = Axis(fig_cal_pos[1,1],
		title = "Kalibratsiooniandmestik",
		xlabel = "RA (kraad)",
		ylabel = "DEC (kraad)",
		xminorticksvisible = true,
		yminorticksvisible = true,
		xgridvisible = false,
		ygridvisible = false,
		xminorticks = IntervalsBetween(4),
		yminorticks = IntervalsBetween(4),
		xminortickalign = 1,
		yminortickalign = 1,
		yminorgridvisible = false,
		xtickalign = 1,
		ytickalign = 1,
		xminorticksize = 10,
		xticksize = 15,
		xreversed = true)
	
	scatter!(ax_cal_pos, df_calibration_raw.ra_epoch2000, df_calibration_raw.dec_epoch2000, color=(:blue, 1), markersize = 5, marker = :star5, label="Välja jäänud tsefeiidid")
	scatter!(ax_cal_pos, df_calibration.ra_epoch2000, df_calibration.dec_epoch2000, color=(:lightblue, 1), markersize = 9, marker = :star5, label="Kalibratsiooniandmestiku tsefeiidid")
	
	axislegend(position=:ct)
	
	fig_cal_pos
	
end

# ╔═╡ 29e0bf05-def9-4373-8ee1-47116faba7fb
begin
	mask_m31_wrangling = (df_m31_raw.Type .∈ Ref(("FM",))) .&
		.!isnan.(df_m31_raw.Pr) .&
		.!isnan.(df_m31_raw.gP1mag) .&
		.!isnan.(df_m31_raw.rP1mag) .& (df_m31_raw.gP1mag .> 16)
	
	df_m31 = df_m31_raw[mask_m31_wrangling, [:Pr, :gP1mag, :rP1mag, :RAJ2000, :DEJ2000]]
end ;

# ╔═╡ 57c63a7d-2fd6-400b-aebf-98ce7e022330
md"M31 andmestikust filtreerime veel välja sissekanded, millel ei ole perioodi või näiva heleduse väärtust ja eemaldame mõned põhigrupist kaugemelae jäävad tähed. Sellesse andmestikku jääb $(size(df_m31, 1)) kirjet."

# ╔═╡ e53b9d50-0bcb-4932-89cf-2d848d938bf6
begin

	fig_appmag_logperiod = Figure(resolution = (1200, 800))
	ax_appmag_logperiod = Axis(fig_appmag_logperiod[1,1],
		title = "PAndromeda",
		xlabel = "log (P(päev))",
		ylabel = "g-riba näivheledus",
		xminorticksvisible = true,
		yminorticksvisible = true,
		xgridvisible = true,
		xminorticks = IntervalsBetween(5),
		yminorticks = IntervalsBetween(3),
		xminortickalign = 1,
		yminortickalign = 1,
		yminorgridvisible = true,
		xtickalign = 1,
		ytickalign = 1,
		xminorticksize = 10,
		xticksize = 15,
		yreversed = true)
	
	scatter!(ax_appmag_logperiod, log10.(df_m31.Pr), df_m31.gP1mag, color=(:lightblue, 1), markersize = 9, marker = :diamond)
	
	fig_appmag_logperiod
	
end

# ╔═╡ ce4d8304-eb47-4ea3-8720-7633ba7f7d43
begin

	fig_m31_pos = Figure(resolution = (1200, 1400))
	ax_m31_pos = Axis(fig_m31_pos[1,1],
		title = "PAndromeda",
		ylabel = "DEC (kraad)",
		xminorticksvisible = true,
		yminorticksvisible = true,
		xgridvisible = false,
		ygridvisible = false,
		xminorticks = IntervalsBetween(4),
		yminorticks = IntervalsBetween(4),
		xminortickalign = 1,
		yminortickalign = 1,
		yminorgridvisible = false,
		xtickalign = 1,
		ytickalign = 1,
		xminorticksize = 10,
		xticksize = 15,
		xreversed = true)
	
	ax_only_m31_pos = Axis(fig_m31_pos[2,1],
		xlabel = "RA (kraad)",
		ylabel = "DEC (kraad)",
		xminorticksvisible = true,
		yminorticksvisible = true,
		xgridvisible = false,
		ygridvisible = false,
		# xminorticks = IntervalsBetween(4),
		# yminorticks = IntervalsBetween(4),
		xminortickalign = 1,
		yminortickalign = 1,
		yminorgridvisible = false,
		xtickalign = 1,
		ytickalign = 1,
		xminorticksize = 10,
		xticksize = 15,
		xreversed = true)
	
	scatter!(ax_m31_pos, df_calibration.ra_epoch2000, df_calibration.dec_epoch2000, color=(:lightblue, 0.6), markersize = 9, marker = :star5)
	scatter!(ax_m31_pos, df_m31.RAJ2000, df_m31.DEJ2000, color=(:blue, 0.5), markersize = 2, marker = :star4)
	
	scatter!(ax_only_m31_pos, df_m31.RAJ2000, df_m31.DEJ2000, color=(:blue, 1), markersize = 9, marker = :star4)


	
	fig_m31_pos
	
end

# ╔═╡ 49e95cb1-7faa-454d-8c02-1d486222eef0
md"""
# Analüüs

Andromeda kauguse leidmiseks on vaja esmalt kalibreerida seos tsefeiidide perioodi ja heleduse vahel.
"""

# ╔═╡ 3de9f4a7-2e79-465a-8405-2428bb2fb9a8
md"Kõigepealt tuleb teisendada parallaks parsekiteks."

# ╔═╡ e663cc9b-aab8-4569-a60e-04f7fd0b04c8
parallax_to_parsecs(p::Real) = u"kpc"/p ;

# ╔═╡ 6a68e1fe-5203-4c70-ab43-e5461fa3f070
md"""Seejärel leiame näivast heledusest absoluutse heleduse kasutades tähe kaugust 
```math
M = m - 5 ( \log_{10} d_{pc} - 1 ) \ .
```
"""

# ╔═╡ 7fa96134-23f2-4c8e-bd20-1097d672401e
abs_mag(app_mag::Real, distance::Unitful.Length) = app_mag - 5. * ( log10.(distance/u"pc") - 1 ) ;

# ╔═╡ 37486888-e4f2-488a-815b-b3ca78152b05
md"""
Näiva ja absoluutse heleduse teisendusel tuleb ka arvestada neeldumisega ehk heledus ei ole vähenenud mitte ainult kaugusest tulenevalt, vaid ka neeldumise tõttu.

V-riba heleduse neeldumise saab leida seosest 
```math
A_\mathrm{V} = 3,\!1 \cdot \mathrm{E(B-V)}
```
[(Wang & Chen 2019)](https://ui.adsabs.harvard.edu/abs/2019ApJ...877..116W/abstract). Suurust E(B-V) nimetatakse punanemiseks ja see on antud G18 andmestikus.
"""

# ╔═╡ aa84b64d-bf5d-42d2-b0d1-ea1fe3770e8a
app_mag_extincion_corrected(app_mag::Real, reddening::Real) = app_mag - 3.1 * reddening ;

# ╔═╡ 8b419528-6c7c-47ea-9460-1bd7f989c5ab
md"Seos tsefeiidide perioodi ja absoluutse heleduse vahel on logaritmiline, seega võtame perioodist kümnendlogaritmi."

# ╔═╡ 94bd487d-3efe-47eb-bf15-90494cecb2ac
log_period(period::Real) = log10(period) ;

# ╔═╡ 0c5a7ffa-0c4d-4089-b036-b0e4aadd9060
md"Leian lineaarse seose absoluutse heleduse ja perioodi logaritmi vahel."

# ╔═╡ 661993c0-0370-4de4-8161-de5e76b8ca4c
function LPR(df::DataFrame, parallax_col::T, app_mag_col::T, reddening_col::T, period_col::T) where T<:Union{String, Symbol}
	distances = parallax_to_parsecs.(df[:, parallax_col])
	app_mags = app_mag_extincion_corrected.(df[:, app_mag_col], df[:, reddening_col])
	abs_mags = abs_mag.(app_mags, distances)
	logperiods = log_period.(df[:, period_col])
	(fit = fit(logperiods, abs_mags, 1), abs_mags = abs_mags, logperiods = logperiods)
end ;

# ╔═╡ aea97c39-ad55-4b84-849c-2ae42689f113
LPR_cal = LPR(df_calibration, :parallax, :vmag, :reddening, :period) ;

# ╔═╡ 94a4ac61-6ee4-46e1-a3f2-9583006aad0c
LPR_cal.fit

# ╔═╡ ef84dc6a-1c73-4fc8-96ca-155de6b7546c
begin

	fig_absmag_logperiod = Figure(resolution = (1200, 800))
	ax_absmag_logperiod = Axis(fig_absmag_logperiod[1,1],
		title = "Kalibratsiooniandmestik",
		xlabel = "log (P(päev))",
		ylabel = "V-riba absoluutne heledus",
		xminorticksvisible = true,
		yminorticksvisible = true,
		xgridvisible = true,
		xminorticks = IntervalsBetween(4),
		yminorticks = IntervalsBetween(3),
		xminortickalign = 1,
		yminortickalign = 1,
		yminorgridvisible = true,
		xtickalign = 1,
		ytickalign = 1,
		xminorticksize = 10,
		xticksize = 15,
		yreversed = true)
	
	
	scatter!(ax_absmag_logperiod, LPR_cal.logperiods, LPR_cal.abs_mags, color=(:lightgreen, 0.5), markersize = 11, marker = :diamond)
	line_ends = [minimum(LPR_cal.logperiods), maximum(LPR_cal.logperiods)]
	lines!(ax_absmag_logperiod, line_ends, LPR_cal.fit.(line_ends), color=:black, linewidth=2)
	
	fig_absmag_logperiod
	
end

# ╔═╡ 220c936f-6c30-4cb6-a6b4-885d8a8c00d5
md"""
PAndromeda andmestikus ei oleks kahjuks V-riba heledusi, seal on Pan-STARRSi g-, r- ja i-riba heledused.
"""

# ╔═╡ 5ccdca31-929a-48a6-bc0f-b6b952482f0b
Show(MIME"image/png"(), read("panstarrs_vs_BVRI.png"))

# ╔═╡ 06fe4a26-5e5a-490e-b94a-8601a4ed654d
md"""Allikas: [Kostov & Bonev (2018)](https://ui.adsabs.harvard.edu/abs/2018BlgAJ..28....3K/abstract)

Õnneks Pan-STARRSi filtrid sarnanevad SDSSi g-, r- ja i-riba filtritele."""

# ╔═╡ cf9aa7c4-f4be-4e71-b3b5-c78427ed3f70
Show(MIME"image/png"(), read("sdss_vs_panstarrs.png"))

# ╔═╡ b231b847-13f4-475a-8739-93f567d8add1
md"Allikas: [Pan-STARRS in the Era of Wide-Field Surveys, Armin Rest, Harvard University](https://slideplayer.com/slide/5931271/)"

# ╔═╡ a0be3cd4-e255-4b90-97be-9cbdf46c45a6
md"""SDSS g- ja r-riba filtreid kasutades on võimalik leida V-riba filtri heledus seosega
```math
\mathrm{V} = \mathrm{g} - 0,\!59 (\mathrm{g} - \mathrm{r}) - 0,\!01
```
[(Transformations between SDSS magnitudes and other systems)](http://www.sdss3.org/dr8/algorithms/sdssUBVRITransform.php). Õnneks PAndromeda andmestiku heledustel on neeldumine juba arvesse võetud.
"""

# ╔═╡ 1a945d9f-4183-468c-aa2c-6c955a742127
gr_to_V(g::Real, r::Real) = g - 0.59*(g-r) - 0.01 ;

# ╔═╡ f36c6585-fb50-4f51-8d91-97bb091a607d
md"""
Eelnevalt leitud V-riba absoluutse heleduse ja perioodi logaritmi vahelise seose abil on võimalik leida Andromedas asuvate tsefeiidide absoluutsed heledused. Näivat ja absoluutset heledust kasutades saame leida kaugusmooduli, mis on nende vahe.
"""

# ╔═╡ 29d2006d-21fb-45f1-b755-d1e85b1afa28
distance_modulus(app_mag::Real, abs_mag::Real) = app_mag - abs_mag ;

# ╔═╡ acb8e054-6403-416d-a2e5-18c11a2cb3f6
md"""Kaugusmooduli μ abil on saame leida kauguse parsekites seosest
```math
d_\mathrm{pc} = 10^{\mu/5 + 1} \ .
```
"""

# ╔═╡ dff6e148-48bf-4b79-bfad-288de7b4604f
distance_from_modulus(distance_modulus::Real) = 10^(distance_modulus/5 + 1) ;

# ╔═╡ 38391fc5-d7bf-43e0-927b-da4edd2d2919
md"Leian Andromeda kauguse."

# ╔═╡ aa6ee1ef-e3ef-45c6-bdbb-0e336447bde6
function m31_distance(df::DataFrame, app_gmag_col::T, app_rmag_col::T, period_col::T, LPR::Polynomial{<:Real}) where T<:Union{String, Symbol}
	app_mags = gr_to_V.(df[:, app_gmag_col], df[:, app_rmag_col])
	abs_mags = LPR.(log_period.(df[:, period_col]))
	dist_modulus = distance_modulus.(app_mags, abs_mags)
	mean_distance_modulus = mean(dist_modulus)
	(distance = round(typeof(1.0u"kpc"), distance_from_modulus.(mean_distance_modulus) * u"pc" |> u"kpc", sigdigits=4), modulus = dist_modulus)
end ;

# ╔═╡ 7e38a11a-4abf-4040-858e-810524cf62ea
function m31_distance(df::DataFrame, app_mag_col::T, period_col::T, LPR::Polynomial{<:Real}) where T<:Union{String, Symbol}
	abs_mags = LPR.(log_period.(df[:, period_col]))
	dist_modulus = distance_modulus.(df[:, app_mag_col], abs_mags)
	mean_distance_modulus = mean(dist_modulus)
	(distance = round(typeof(1.0u"kpc"), distance_from_modulus.(mean_distance_modulus) * u"pc" |> u"kpc", sigdigits=4), modulus = dist_modulus)
end ;

# ╔═╡ cdd17294-77fb-4e11-9725-b889a1fa48a3
md"# Andromeda kauguse arvutamine"

# ╔═╡ d6bbca13-374b-4966-827c-1fb8c42aff44
md"""
Kasutades 
$(@bind band Select(["(:gP1mag, :rP1mag)" => "V-riba tuletatud heledusi", "(:gP1mag,)" => "g-riba heledusi", "(:rP1mag,)" => "r-riba heledusi"]))
saame Andromeda kauguseks
"""

# ╔═╡ 8ef83720-958a-4804-beb1-898d829a1a48
md"**$(m31_distance(df_m31, eval(Meta.parse(band))..., :Pr, LPR_cal.fit).distance)**"

# ╔═╡ 54d57b23-8b0e-405e-90d0-d10f64d0882a
md"Wikipediast võetud [Andromeda](https://en.wikipedia.org/wiki/Andromeda_Galaxy) kaugus on 765 kpc."

# ╔═╡ Cell order:
# ╟─5a96edd0-a432-11eb-119b-bb814232c52a
# ╟─112f6c01-5ef7-4980-81ce-282b242a9cf2
# ╟─9e747c0c-ab6d-4a15-aa6d-1d17716f06f3
# ╟─5e95a0fd-9f0b-40a8-9cc8-f945160d28e8
# ╟─5ef8ac87-f667-4ecc-8d57-d7dbf71a3e5f
# ╟─b0e70be5-17c8-4d84-a23c-84e0526d9e7b
# ╟─f8002ac9-476d-470e-9a72-3d0656f58412
# ╟─8ccdaf15-5af0-441e-8eee-42c5caa00a9c
# ╟─70c39032-3255-4f46-97a0-70e53957b6d3
# ╟─926dc1d8-5eaf-4a8b-8f31-a40a5180941b
# ╟─0a334697-f3a6-4b28-8f6e-3d0add984e05
# ╟─22d57db1-8374-41f0-bd7a-ed1aa4eb4889
# ╠═dc8bfc19-86c7-4d78-a480-cf4590f12dc2
# ╟─c6d59e16-70f8-47b3-93c8-75247f6dad6f
# ╟─57c63a7d-2fd6-400b-aebf-98ce7e022330
# ╟─29e0bf05-def9-4373-8ee1-47116faba7fb
# ╟─e53b9d50-0bcb-4932-89cf-2d848d938bf6
# ╟─ce4d8304-eb47-4ea3-8720-7633ba7f7d43
# ╟─49e95cb1-7faa-454d-8c02-1d486222eef0
# ╟─3de9f4a7-2e79-465a-8405-2428bb2fb9a8
# ╠═e663cc9b-aab8-4569-a60e-04f7fd0b04c8
# ╟─6a68e1fe-5203-4c70-ab43-e5461fa3f070
# ╠═7fa96134-23f2-4c8e-bd20-1097d672401e
# ╟─37486888-e4f2-488a-815b-b3ca78152b05
# ╠═aa84b64d-bf5d-42d2-b0d1-ea1fe3770e8a
# ╟─8b419528-6c7c-47ea-9460-1bd7f989c5ab
# ╠═94bd487d-3efe-47eb-bf15-90494cecb2ac
# ╟─0c5a7ffa-0c4d-4089-b036-b0e4aadd9060
# ╠═661993c0-0370-4de4-8161-de5e76b8ca4c
# ╠═aea97c39-ad55-4b84-849c-2ae42689f113
# ╟─94a4ac61-6ee4-46e1-a3f2-9583006aad0c
# ╟─ef84dc6a-1c73-4fc8-96ca-155de6b7546c
# ╟─220c936f-6c30-4cb6-a6b4-885d8a8c00d5
# ╟─5ccdca31-929a-48a6-bc0f-b6b952482f0b
# ╟─06fe4a26-5e5a-490e-b94a-8601a4ed654d
# ╟─cf9aa7c4-f4be-4e71-b3b5-c78427ed3f70
# ╟─b231b847-13f4-475a-8739-93f567d8add1
# ╟─a0be3cd4-e255-4b90-97be-9cbdf46c45a6
# ╠═1a945d9f-4183-468c-aa2c-6c955a742127
# ╟─f36c6585-fb50-4f51-8d91-97bb091a607d
# ╠═29d2006d-21fb-45f1-b755-d1e85b1afa28
# ╟─acb8e054-6403-416d-a2e5-18c11a2cb3f6
# ╠═dff6e148-48bf-4b79-bfad-288de7b4604f
# ╟─38391fc5-d7bf-43e0-927b-da4edd2d2919
# ╠═aa6ee1ef-e3ef-45c6-bdbb-0e336447bde6
# ╟─7e38a11a-4abf-4040-858e-810524cf62ea
# ╟─cdd17294-77fb-4e11-9725-b889a1fa48a3
# ╟─d6bbca13-374b-4966-827c-1fb8c42aff44
# ╟─8ef83720-958a-4804-beb1-898d829a1a48
# ╟─54d57b23-8b0e-405e-90d0-d10f64d0882a
