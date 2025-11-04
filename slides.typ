#import "@preview/mannot:0.3.0": *
#import "touying/lib.typ": *
#import "@preview/pinit:0.1.4": *
#import "@preview/xarrow:0.3.0": xarrow
#import "@preview/cetz:0.4.1"
#import "psi-slides-0.6.1.typ": *
#import "@preview/grayness:0.3.0": image-grayscale
#import "@preview/algorithmic:1.0.3"
#import algorithmic: algorithm

// color-scheme can be navy-red, blue-green, or pink-yellow
// #let s = psi-slides.register(aspect-ratio: "16-9", color-scheme: "pink-yellow")
#show: psi-theme.with(aspect-ratio: "16-9",
                      color-scheme: "pink-yellow",
                             config-info(
                                title: [Correcting the failures of DFT],
                                subtitle: [Koopmans functionals, DFT+_U_, and more],
                                author: [Edward Linscott],
                                date: datetime(year: 2025, month: 10, day: 30),
                                location: [PsiQuantum],
                                references: [references.bib],
                             ))

#set footnote.entry(clearance: 0em)
#show bibliography: set text(0.6em)

#let primary = rgb("#dc005a")
#let secondary = rgb("#f0f500")

#let blcite(reference) = {
  text(fill: white, cite(reference))
}

#let delayedmark(start, content, tag: none, color: red, mark: mark, color-before: black, alternatives: none) = {
   let entries = (mark(content, tag: tag, color: color-before),)*(start - 1) + (mark(content, tag: tag, color: color),)
   alternatives(repeat-last: true, ..entries)
}

#let methods-with-marks(self) = {
  let (uncover, only, alternatives) = utils.methods(self)
  let dm = delayedmark.with(alternatives: alternatives)
  (uncover, only, alternatives, dm, dm.with(mark: markhl, color-before: white))
}

#title-slide()
#matrix-slide()[
  #image("media/photos/otago.jpg", height: 100%)
  #place(center + bottom, text(weight: "bold", fill: white, [Otago, NZ]) + v(0.5em))
  #pause
][
  #image("media/photos/cambridge.jpg", height: 100%)
  #place(center + bottom, text(weight: "bold", fill: white, [Cambridge, UK]) + v(0.5em))
  #pause
][
  #image("media/photos/epfl.jpg", height: 100%)
  #place(center + bottom, text(weight: "bold", fill: white, [EPFL, CH]) + v(0.5em))
  #pause
][
  #image("media/photos/psi_in_mist_large.jpg", height: 100%)
  #place(center + bottom, text(weight: "bold", fill: white, [PSI, CH]) + v(0.5em))
  #pause
][
  #image("media/photos/brisbane.png", height: 100%)
  #place(center + bottom, text(weight: "bold", fill: white, [... PsiQuantum \ AU?]) + v(0.5em))
]


== Predicting electronic excitations

Spectral properties are fundamental to understanding molecules and materials...

#align(center, 
grid(columns: 2, column-gutter:  1em,
//cetz.canvas({
//  import cetz.draw: *
//
//  // grid((0,-5), (8,5), stroke: gray + .5pt)
//
//  // Valence
//  rect((-1, -1), (1, 1), stroke: none, fill: secondary, alpha: 0.5)
//  content((1.75, 1), [$E_F$], align: left)
//  circle((0, 0), radius: 0.2, fill: white, stroke: (dash: "dashed", paint: primary))
//
//  // Vacuum
//  circle((0, 4), radius: 0.2, fill: primary, stroke: none)
//  line((-1, 3.5), (1, 3.5), stroke: (dash: "dashed", paint: primary))
//  content((1.75, 3.5), [$E_"vac"$], align: left)
//
//  // Arrow
//  arc((0,0), start: -30deg, stop: 30deg, radius: 4, mark: (end: ">", fill: black))
//  
//  let photon(amplitude: 1, phases: 2, scale: 8, samples: 1000, angle: 0, start-x: 0, start-y: 0, ..args) = {
//    line(..(for x in range(0, samples + 1) {
//      let x = x / samples
//      // A Gaussian envelope with sigma = 1/4 and mean = 1/2 and height = amplitude
//      let envelope = amplitude * calc.exp(-calc.pow(((x - 0.5) / (0.25)), 2))
//
//      let phase = (2 * phases * calc.pi) * x
//
//      // Rotate the output by angle
//      let xval = x * scale
//      let yval = calc.sin(phase) * envelope
//
//      let rotated-x = xval * calc.cos(angle) - yval * calc.sin(angle)
//      let rotated-y = xval * calc.sin(angle) + yval * calc.cos(angle)
//      ((start-x + rotated-x, start-y + rotated-y),)
//    }), ..args)
//
//    let subdivs = 8
//    for phase in range(0, phases) {
//      let x = phase / phases
//      for div in range(1, subdivs + 1) {
//        let p = 2 * calc.pi * (div / subdivs)
//        let y = calc.sin(p) * amplitude
//        let x = x * scale + div / subdivs * scale / phases
//      }
//    }
//  }
//  photon(amplitude: 0.8, phases: 9, start-x: -0.25, start-y: 0.25, scale: 3, fill: none, angle: 2.5, mark: (start: ">", fill: black))
//}),
  image("figures/nguyen_2015_buckyball_pes.png", height: 40%),
  image("figures/arpes_puppin.png", height: 40%),
))

#blcite(<Nguyen2015>)#blcite(<Puppin2020>)
#v(-2.5em)
#pause ... but how can we routinely compute them? #pause

- quantum chemistry: gold standard, but scales prohibitively (for now...)
- GW: accurate but expensive, often ill-behaved, diagrammatic
- DFT: plagued by intrinsic errors
#pause

*My work: understanding and correcting these intrinsic errors*

= The failures of DFT
// == Learning nothing from Icarus
// 
// Kohn-Sham DFT introduces an auxiliary non-interacting system such that...
// 
// $ rho^"DFT" (bold(r)) = rho^"exact" (bold(r)) $
// 
// and
// 
// $ E^"DFT" [rho(bold(r))] = E^"exact" $
// 
// What does the single-particle spectrum of this auxiliary system look like?

== Learning nothing from Icarus
From the auxiliary, non-interacting system:

#align(center, image("figures/onida_2002_cu_bs.png", height: 65%))

#pause
#align(right, [... the temptation is too much!])

#v(-1em)
#blcite(<Onida2002>)

== Learning nothing from Icarus
- single-particle excitation energies can be related to total energy differences
$ - epsilon_"HOMO" = I = E(N - 1) - E(N) $ #pause
- the long-range decay of the density is proportional to $epsilon_"HOMO"$@Almbladh1985 #pause
- the KS xc potential is the best local approx. to the xc self-energy@Casida1995 #pause

#slide(repeat: 2, self => [
  #let (uncover, only, alternatives, delayedmark, delayedmarkhl) = methods-with-marks(self)


#align(horizon + center,
alternatives()[
#image("figures/ip_1.svg", width: 90%)
][
#image("figures/ip_2.svg", width: 90%)
][
]
)

#blcite(<Dabo2010>)
])

== What's going wrong?

#align(horizon,
grid(align: horizon, columns: (1fr, 1fr), column-gutter: 1em,
  [
// The exact Green's function has poles that correspond to total energy differences

$
  ε_i =^? cases(E(N) - E_i (N-1) & "if" i in "occ", E_i (N+1) - E(N) & "if" i in "emp")
$

#pause
_cf._ Janak's theorem:
$
  ε_i^"DFT" = frac(dif E^"DFT", dif f_i)
$

#pause


// but DFT does #emph[not]
],[
  #image(width: 20em, "figures/fig_en_curve_gradients_zoom.svg")
]
))

= Aside: piecewise linearity
== Aside: piecewise linearity

#pause
- Perdew, Parr, Levy, and Balduz@Perdew1982a showed that the exact total energy is piecewise linear in electron number (using ensembles) #pause
- Yang, Zhang, and Ayers@Yang2000 provided an alternative proof without invoking ensembles #pause
- We showed@Burgess2023b that the "convexity condition" follows _i.e._
  $ 2 E(N) <= E(N + 1) + E(N - 1) $
  for all DFTs that are
  - exact for all $v$-representable densities
  - size-consistent
  - translationally invariant #pause
- ... and that similar reasoning applies to total energy as a function of total magnetisation@Burgess2024a

#focus-slide()[Core idea: enforce piecewise linearity]

== Imposing generalised piecewise linearity
#slide(repeat: 3, self => [
  #let (uncover, only, alternatives, delayedmark, delayedmarkhl) = methods-with-marks(self)

#align(horizon,
grid(align: horizon, columns: (1fr, auto), column-gutter: 1em,
[
  $
  E^"KI" &[rho, {rho_i}] =
  E^"DFT" [rho]
  \ & +
  sum_i (
    - delayedmarkhl(#2, integral_0^f_i lr(chevron.l phi_i mid(|) hat(h)^"DFT" (f) mid(|) phi_i chevron.r) dif f, tag: #<remove_nonlin>, color: primary)
  \ &
    + delayedmarkhl(#3, f_i integral_0^1 lr(chevron.l phi_i mid(|) hat(h)^"DFT" (f) mid(|) phi_i chevron.r) dif f, tag: #<restore_linear>, color: primary)
  )
$
// Bakes the total energy differences $E^"DFT" [rho^(f_i arrow.r 1)] - E^"DFT" [rho^(f_i arrow.r 0)]$ into the functional
#uncover("2-")[#annot(<remove_nonlin>, pos: bottom)[#align(center, [removes dependence on $f_i$])]]
#uncover("3-")[#annot(<restore_linear>, pos: bottom)[#align(center, [restores linear dependence on $f_i$])]]

],[
  #image(width: 20em, "figures/fig_en_curve_gradients_zoom.svg")
]))

])

== Details for the experts
$
  E^"KI"_bold(alpha) [rho, {rho_i}] =
  E^"DFT" [rho] +
  sum_i alpha_i { &
    - (E^"DFT" [rho] - E^"DFT" [rho - rho_i])
  \ &
    + f_i (E^"DFT" [rho - rho_i + n_i] - E^"DFT" [rho - rho_i])
  }
$

- orbital-density-dependence
- screening parameters
- total energy at integer occupations unchanged!

= Results

== Atoms
#slide(repeat: 2, self => [
  #let (uncover, only, alternatives, delayedmark, delayedmarkhl) = methods-with-marks(self)

#align(center + horizon,
  alternatives[
    #image("figures/ip_2.svg", width: 100%)
  ][
    #image("figures/ip_3.svg", width: 100%)
  ]
)

#blcite(<Dabo2010>)
])
== Molecules

=== Ionisation potentials@Colonna2019
#align(center + horizon,
image("figures/colonna_2019_gw100_ip.jpeg", width: 100%)
)

=== UV photoemission spectra@Nguyen2015
#align(center + horizon,
image("figures/fig_nguyen_prl_spectra_pink.png", width: 100%)
)


== Materials
#slide[
=== Prototypical semiconductors and insulators @Nguyen2018

#show table.cell: it => {
  if it.x == 3 or it.x == 4 {
    set text(fill: primary, weight: "semibold")
    it
  } else {
    it
  }
}

#grid(align: center + horizon, columns: 2, column-gutter: 1em,
image("figures/scatter_plot.png", height: 80%),
table(columns: (auto, 1fr, 1fr, 1fr, 1fr, 1fr), inset: 0.5em, stroke: none,
table.header([], [PBE], [G#sub[0]W#sub[0]], [KI], [KIPZ], [QSGW̃]),
table.hline(),
[$E_"gap"$], [2.54], [0.56], [0.27], [0.22], [0.18],
[IP], [1.09], [0.39], [0.19], [0.21], [0.49]
))
  
]

  
== Photocatalysis
#slide(repeat: 2, self => [

  #let (uncover, only, alternatives, delayedmark, delayedmarkhl) = methods-with-marks(self)


#align(center + horizon, 
alternatives[
  #grid(columns: 2,
  image("figures/wpc_sketch.png", height: 50%),
  image("figures/anatase_water_slab.png", height: 50%)
  )
][
  #image("figures/ba_fixed.png", height: 90%)
]
)

#v(-2em)
#blcite(<Stojkovic2024>)

])

== Optical spectra

Solve the BSE, using Koopmans eigenvalues in lieu of GW

#pause

#v(-1em)
#align(center + horizon,
grid(columns: 2,
image("figures/silicon_bse_spectra.png", height: 50%),
image("figures/silicon_bse_excitons.png", height: 50%)
))

#v(-1em)

#show table.cell: it => {
  set text(size: 0.8em)
  it
}
#table(align: center + horizon, columns: (auto, 1fr, 1fr, 1fr, 1fr), inset: 0.5em, stroke: none,
table.header([silicon], [indirect gap ], [direct gap ], [first excitonic peak ], [excitonic binding energy ]),
table.hline(),
[*qKI+BSE*], [1.12], [3.31], [3.42], [0.09], 
[G#sub[0]W#sub[0]+BSE], [1.17], [3.25], [3.34], [0.09],
)


// == 
// #image("figures/supercell_workflow.png", width: 100%)
// 
// #image("figures/primitive_workflow.png", width: 65.5%)

#focus-slide()[
#align(center, image(width: 80%, "media/logos/koopmans_white_on_transparent.svg"))
]

// #matrix-slide(alignment: horizon)[
//   #image("figures/website_cropped.png")
// ][
//   
//   - automated workflows
//   - `Quantum ESPRESSO` backend
//   - easy installation
//   - python API
//   
//   See `koopmans-functionals.org`
// ]

#matrix-slide(alignment: horizon, columns: (3fr, 2fr))[
  #image("figures/black_box_filled_square.png")
][
 
  Our goal:
  + accurate
  + robust
  + minimal input
  + fast

]


== 
#slide()[
#show table.cell: it => {
  if it.x == 5 {
    set text(fill: primary, weight: "semibold")
    it
  } else {
    it
  }
}
#grid(columns: 3, align: horizon + center,
  [
  #set text(size: 0.45em)
  #raw(read("scripts/gaas_auto.json"), block: true, lang: "json")
  ],
  [
    #set text(size: 3em)
    #sym.arrow.r
  ],
  [
    #image("figures/Unfold_And_Interpolate_bandstructure.png", height: 60%)
    #table(columns: (auto, 1fr, 1fr, 1fr, 1fr, 1fr, 1fr), inset: 0.5em, stroke: none,
    table.header([], [LDA], [HSE], [GW#sub[0]], [scGW̃ ], [KI], [exp]),
    table.hline(),
    [$E_"gap"$], [0.26], [1.28], [1.55], [1.62], [1.54], [1.55],
    [$chevron.l epsilon_d chevron.r$], [-14.9], [-15.6], [-17.3], [-17.6], [-17.9], [-18.9],
    [$Delta$], [12.8], [13.9], [], [], [12.7], [13.1]
    )
  ]
)
]

== `koopmans`
#grid(columns: (2fr, 3fr), align: horizon,
  image("figures/jctc_cover_annotated.jpg", height: 90%),
  [
    - used by a Fortune Global 500 company
    - two schools (online and then in Pavia, IT)

    See `koopmans-functionals.org`

    #blcite(<Linscott2023>)
  ]
)
= An alternative approach: DFT + _U_

#slide(repeat: 5, self => [

  #let (uncover, only, alternatives, delayedmark, delayedmarkhl) = methods-with-marks(self)
  #set text(size: 2em)

  #align(horizon,
    [$ E_("DFT"+U) = E_"DFT" + sum_(delayedmark(#2, I sigma, tag: #<site_index>, color: primary)) delayedmark(#5, U^I, tag: #<hubbard_u>, color: primary) / 2 "Tr"[delayedmark(#3, bold(n)^(I sigma), tag: #<occ_matrix>, color: primary) (1 - bold(n)^(I sigma))] $]
  )
  #set text(size: 0.666em)
  #uncover("4-")[
  $ n_(m m')^(I sigma) = chevron.l phi^I_m|hat(rho)^(sigma)|phi^I_(m')chevron.r = sum_i chevron.l phi^I_m|psi_i chevron.r f_i chevron.l psi_i|phi^I_(m')chevron.r $
  ]

  #pause
  #annot(<site_index>, pos: bottom)[#align(center, [site and spin indices])]
  #pause
  #annot(<occ_matrix>, pos: bottom)[#align(center, [local occupation \ matrix])]
  #pause
  #pause
  #annot(<hubbard_u>, pos: top)[#align(center, [Hubbard parameter])]

  #v(-3em)

  #blcite(<Linscott2018>)

])

== The historical derivation of DFT+_U_
#align(center + horizon,
  image("figures/tmos.svg", width: 50%)
)
#blcite(<Rodl2009>)

== The historical derivation of DFT+_U_

#slide(self => [
  #set text(size: 0.5em)
  #set math.equation(numbering: "(1)")
  #set highlight(fill: primary.lighten(50%))

  Let a correlated subspace be defined by a set of basis orbitals (known as _Hubbard projectors_). Within this subspace, the operator associated with electron-electron interactions is

  $
     hat(U) = sum_(m n m'n') sum_(sigma sigma')U_(m n m'n')hat(c)^dagger_(m sigma) hat(c)^dagger_(n sigma') hat(c)_(m'sigma') hat(c)_(n'sigma),
  $
  where $(m,n,m',n')$ are Hubbard projector labels and ${sigma}$ are spin indices, and $hat(c)^dagger_(m sigma)$ are the associated creation operators. One can show that
  $
     E_"Hub" = chevron.l hat(U) chevron.r = &  1 / 2 sum_(m n m' n' sigma \ m != n, m' != n')
     (U_(m n m' n')-U_(m n n' m'))chevron.l n',sigma;m',sigma|hat(rho)_2|n,sigma;m,sigma chevron.r \
                                             & + 1 / 2 sum_(m n m'n' sigma)
     U_(m n m'n')chevron.l n',sigma;m',-sigma|hat(rho)_2|n,-sigma;m,sigma chevron.r -U_(m n n' m')chevron.l n',-sigma;m',sigma|hat(rho)_2|n,-sigma;m,sigma chevron.r
  $

  where $hat(rho)_2$ is the two-body density matrix. Adopting the ansatz that the #box(place(bottom, text(size: 2em, fill: primary, weight: "bold")[Hartree-Fock approximation], dy: -1em))#highlight(fill: primary.lighten(50%))[many-body wavefunction is a Slater determinant of single-particle states], the two-body density matrices $hat(rho)_2$ can be decomposed as determinants of single-body density.@Parr1989a In this case

  $
    E_"Hub"
    = & 1/2 sum_(m n m'n' sigma \ m != n, m' != n')
    (U_(m n n'm')-U_(m n m' n'))n^(sigma)_(m m') n^(sigma)_(n n')
    +1/2 sum_(m n m'n' sigma)U_(m n m'n')n_(m n')^(sigma)n_(n m')^(-sigma),
  $ <equation_EU1>
  where $n_(m m')^(sigma) = chevron.l m | hat(rho)^(sigma)|m'chevron.r$. At this stage the only approximation that has been introduced is the assertion that the state corresponds to a Slater determinant. If $U_(m n m' n')$ is obtained using the unscreened Coulomb potential, then @equation_EU1 is equivalent to a Hartree-Fock treatment of the system.

  // #annot(<hf_approx>, pos: bottom)[#align(center, [Hartree-Fock approx.])]
  Now, #box(place(bottom, text(size: 2em, fill: primary, weight: "bold")[two-site terms only], dy: -1em))#highlight[all but two-site terms are ignored]. Due to the symmetries of $U_(m n m'n')$, this leaves only two types of terms: $U_(m n n m)$ and $U_(m n m n)$. These are then #box(place(top, text(size: 2em, fill: primary, weight: "bold")[average a rank-2 \ tensor], dy: 0.5em))#highlight[averaged] over the Hubbard projectors to yield two scalars:
  $
     U = 1 / ((2l+1)^2)sum_(m n)U_(m n n m);
     J = 1 / ((2l+1)^2)sum_(m n)U_(m n m n).
  $
  Using these average values in place of the tensorial terms simplifies @equation_EU1 to
  $
     E_"Hub"
     = & 1 /2 sum_(m n sigma)
     U(n^(sigma)_(m m) n^(sigma)_(n n)-n^(sigma)_(m n) n^(sigma)_(n m)+n^(sigma)_(m m) n^(-sigma)_(n n))
     +1/2 sum_(m n sigma)J(n_(m n)^(sigma)n_(n m)^(sigma)-n_(m m)^(sigma)n_(n n)^(sigma)+n_(m n)^(sigma)n_(n m)^(-sigma))
     = & sum_sigma U / 2 ((n^sigma)^2+n^sigma n^(-sigma)-"Tr"(bold(n)^sigma bold(n)^sigma))
     + J / 2 ("Tr"(bold(n)^sigma bold(n)^sigma+bold(n)^sigma bold(n)^(-sigma))-(n^sigma)^2)
  $ <equation_EU1b>
  where $n^sigma = "Tr"(bold(n)^sigma)$. If at this stage @equation_EU1b was to be incorporated directly into the DFT formalism, interactions associated with the subsystems that are already being handled by the conventional exchange-correlation functional would be double-counted. To avoid this, the #box(place(top, text(size: 2em, fill: primary, weight: "bold")[adopt some double-counting term], dy: 0.5em))#highlight[fully localised limit]@Petukhov2003a is considered, where all correlated subspaces have integer occupancy. In this approximation
  $
     "Tr"(bold(n)^sigma bold(n)^sigma) arrow.r n^sigma;  "Tr"(bold(n)^sigma bold(n)^(-sigma)) arrow.r n^(sigma_"min"),
  $
  where $sigma_"min"$ denotes the minority spin. Thus in the fully localised limit, the double counting term becomes
  $
     E_"DC" = U / 2 n(n-1) - J / 2 sum_sigma n^sigma(n^sigma-1)+J n^(sigma_"min")
  $
  where $n=sum_(sigma) n^sigma$. Hence
  $
     E_"Hub" - E_"DC" = sum_(I sigma) (U^I-J^I)/(2)"Tr"(bold(n)^(I sigma)(1-bold(n)^(I sigma))) + sum_(I sigma)(J^I)/(2)("Tr"(bold(n)^(I sigma)bold(n)^(I-sigma))-2delta_(sigma sigma_(min))n^(I sigma)).
  $
  Note that the entire expression has now been generalised to allow for the possibility of multiple sites (labelled with the index $I$), to each of which a correction term is applied.
  As a final approximation, #box(place(top, text(size: 2em, fill: primary, weight: "bold")[neglect _J_], dy: 0.5em))#highlight[terms arising from interaction between opposite spin (those contained in the second sum) are neglected]. This leaves
  $
     E_(U) = E_"Hub" - E_"DC" = sum_(I sigma) (U^I_"eff") / (2)"Tr"(bold(n)^(I sigma)(1-bold(n)^(I sigma))),
  $ <equation_EU3>
  where the on-site Coulomb repulsion parameter $U^I$ has been effectively reduced by $J^I$ to $U^I_"eff"$. The DFT+_U_ correction to the KS potential is given by
  $
     hat(V)_(U) = sum_(I sigma m n)U^(I) |m chevron.r (1/2-n_(m n)^(I sigma) ) chevron.l n|.
  $
  // With this, the derivation is complete: the Hubbard-model formalism is in a form which can be incorporated into the framework of (DFT).
  // onslide<2->(
  //    begin(tikzpicture)[overlay,remember picture, seaborn_bright_magenta, font=normalsize]
  //       fill[opacity=0.3] (neglectj.south west) rectangle (neglectj.north east);
  //       draw[decorate] (neglectj.south) to[out=-40,in=180] ++ (0.6,-0.3) node[right, text width=6cm](neglect ``$J$" / the second term));
  //    end(tikzpicture)
  // )

])

== The modern interpretation
#grid(align: horizon, columns: (1fr, 1fr),
  image(
  width: 100%,
  "figures/fig_en_curve_dftu_correction.svg"
),
[
  In a basis such that $n^(I sigma)_(i j) = lambda^(I sigma)_i delta_(i j)$,

  $
    E_U = sum_(I sigma) (U^I) / 2 sum_i lambda^(I sigma)_i (1 - lambda^(I sigma)_i)
  $

  #pause

  $U$ can be measured by linear response@Cococcioni2005a #pause -- *critical for predictive calculations*
]
)

==
#slide(repeat: 3, self => [

  #let (uncover, only, alternatives, delayedmark, delayedmarkhl) = methods-with-marks(self)

  #align(horizon,
    grid(columns: (1fr, 1fr), align: horizon + center, inset: 1em,
    [$ E_("DFT"+U) = E_"DFT" + sum_(delayedmark(#2, I sigma, color: primary)) U^delayedmark(#3, I, color: secondary) / 2 "Tr"[bold(n)^(delayedmark(#2, I sigma, color: primary)) (1 - bold(n)^(delayedmark(#2, I sigma, color: primary)))] $],
    [$ U^delayedmark(#3, I, color: secondary) = [chi_0^(-1) - chi^(-1)]_(I I) $],
    [#pause functional treats spin channels separately],
    [#pause LR treats them together],
    )
  )

  #blcite(<Linscott2018>)
])

==

#align(top,
  image("figures/prb2018_pdf_cropped.png", width: 100%, height: 100%)
)

== Spin-resolved linear response

#slide(repeat: 8, self => [

  #let (uncover, only, alternatives, delayedmark, delayedmarkhl) = methods-with-marks(self)

  #grid(columns: (1fr, auto, 1fr), align: horizon + center, inset: 1em,
  [conventional], [#sym.arrow.r], [spin-resolved],
[#pause $ chi_(I J) = (d n^I) / (d v^J)$],
[#pause #sym.arrow.r],
[$chi_(I J)^(sigma sigma') = (d n^(I sigma)) / (d v^(J sigma')) $],
[#pause $ mat(chi_(1 1), chi_(1 2), dots; chi_(2 1), chi_(2 2), dots; dots.v, dots.v, dots.down)$],
[#pause #sym.arrow.r],
[$
mat(delayedmarkhl(#6, mat(delim: #none, chi^(arrow.t arrow.t)_(1 1), chi^(arrow.t arrow.b)_(1 1); chi^(arrow.b arrow.t)_(1 1), chi^(arrow.b arrow.b)_(1 1)), color: primary), 
    mat(delim: #none, chi^(arrow.t arrow.t)_(1 2), chi^(arrow.t arrow.b)_(1 2), dots; chi^(arrow.b arrow.t)_(1 2), chi^(arrow.b arrow.b)_(1 2), dots);
    mat(delim: #none, chi^(arrow.t arrow.t)_(2 1), chi^(arrow.t arrow.b)_(2 1); chi^(arrow.b arrow.t)_(2 1), chi^(arrow.b arrow.b)_(2 1); dots.v, dots.v), 
    mat(delim: #none, chi^(arrow.t arrow.t)_(2 2), chi^(arrow.t arrow.b)_(2 2), dots; chi^(arrow.b arrow.t)_(2 2), chi^(arrow.b arrow.b)_(2 2), dots; dots.v, dots.v, dots.down))
   
 $],
 [#pause #pause $U^I = [chi_0^(-1) - chi^(-1)]_(I I)$],
 [#pause #sym.arrow.r],
 [$f^(sigma sigma')_(I I) arrow.r^? U^(I sigma)$],
  )


])


== Advantages of spin-resolved LR
#pause
- conceptual consistency (spin-resolved functional #sym.arrow.l.r spin-resolved linear response) #pause
- can recover the conventional linear response results #pause
- $J$ is "free" #pause
- easily implemented #pause
- can perform unconstrained constrained linear response #pause

  _e.g._ suppose we want to compute $ lr((d^2E_"Hxc") / (d (n^I)^2) |)_(mu^I)$ #pause

  This is easy with spin-resolved LR: $lr((d^2E_"Hxc") / (d n^2) |)_(mu) = & 1/4 (f^(arrow.t arrow.t) + f^(arrow.b arrow.b) + f^(arrow.t arrow.b) + f^(arrow.b arrow.t))$
  #pause
- ability to pursue more flexible, tailored corrections

// = Initial results
// ==
// #slide(repeat: 5, self => [
// 
//   #let (uncover, only, alternatives, delayedmark, delayedmarkhl) = methods-with-marks(self)
// 
//   #align(center, [
// 
//   #alternatives()[
//     #image("figures/fig_combined_U_J_and_occ.svg", height: 50%)
//   ][
//     #image("figures/fig_MnO_magmom.svg", height: 70%)
//   ][
//     #image("figures/fig_MnO_bandgap.svg", height: 70%)
//   ]
//   ])
//   #v(-1em)
//   #pause
//   #pause
//   #pause
// #align(right, [... for more details see Linscott et al. 2018. #pause Since then used by many authors@Orhan2020@Lambert2023@MacEnulty2023@Moore2024@MacEnulty2024 and opened the door to DFT+$U$-inspired approaches@Burgess2023@Burgess2024a])
// 
// ])

== BLOR
#align(center,
  image("figures/blor.svg", width: 80%)
)

a DFT+_U_ type functional that...
- is inspired by the intrinsic errors of approximate DFT
- relies on spin-resolved linear response
- includes a term to correct for static correlation error
- is double-counting-free

#v(-2em)
#blcite(<Burgess2023>)#blcite(<Burgess2024a>)

== BLOR

#grid(columns: (2fr, 3fr), align: center + horizon, gutter: 1em,
  image("figures/blor_on_h2.svg", width: 100%),
  image("figures/blor_on_h5+.svg", width: 100%),
  [stretched H#sub[2]],
  [stretched H#sub[5]#super[+]]
)
#blcite(<Burgess2023>)#blcite(<Burgess2024a>)

== Spin-resolved LR used in Materials Project!
#align(top,
  image("figures/ptable_HubbardU.svg", width: 100%, height: 80%)
)
#blcite(<Moore2024>)


== Summary

Understanding and correcting the failures of approximate DFT can yield simple but predictive functionals

#align(center,
grid(columns: 3,
  image("figures/fig_en_curve_gradients_zoom.svg", height: 40%),
  image("figures/scatter_plot.png", height: 40%),
  image("figures/blor_on_h5+.svg", height: 40%),
)
)

== Acknowledgements
#slide()[

#set text(size: 0.6em)
#align(center + horizon, 
grid(columns: (10%, 10%, 10%, 10%, 10%), column-gutter: 0.5em, align: center, row-gutter: 0.5em,
  image("media/mugshots/david_oregan.jpg", height: 30%),
  image("media/mugshots/andrew_burgess.jpeg", height: 30%),
  image("media/mugshots/nicola_colonna.png", height: 30%),
  image("media/mugshots/miki_bonacci.jpg", height: 30%),
  image("media/mugshots/aleksandr_poliukhin.jpg", height: 30%),
  [David O'Regan], [Andrew Burgess], [Nicola Colonna], [Miki Bonacci], [Aleksandr Poliukhin],
  image("media/mugshots/marija_stojkovic.jpg", height: 30%),
  image("media/mugshots/junfeng_qiao.jpeg", height: 30%),
  image("media/mugshots/yannick_schubert.jpg", height: 30%),
  image("media/mugshots/nicola_marzari.jpeg", height: 30%),
  [... and many others!],
  [Marija Stojkovic], [Junfeng Qiao], [Yannick Schubert], [Nicola Marzari]
)
)

#align(
  center,
  grid(
    columns: 4,
    align: horizon + center,
    gutter: 2em,
    image("media/logos/royal_society.jpg", height: 15%),
    image("media/logos/epsrc_logo.png", height: 15%),
    image("media/logos/SNF_logo_standard_web_color_pos_e.svg", height: 15%),
    image("media/logos/marvel_color_on_transparent.png", height: 15%),
  ),
)
]

#focus-slide()[#align(center, text(size: 2em, [Thank you!]) + linebreak() + text(size: 0.5em, style: "italic", [these slides are available at #h(0.2em) #box[#move(dy: 0.1em, image("media/logos/github-mark-white.svg", height: 1em))] `elinscott-talks/psi_quantum_presentation`]))]

#show: appendix

#focus-slide()[#align(center, text(size: 2em, [spare slides]))]

== Imposing generalised piecewise linearity
#align(horizon,
grid(align: horizon, columns: (1fr, auto), column-gutter: 1em,
  [Formally, every orbital $i$ should have an eigenenergy
  $
    epsilon_i^"Koopmans" = ⟨
      phi_i mid(|)hat(H)mid(|)phi_i
    ⟩ = frac(dif E, dif f_i)
  $
  that is
  - independent of $f_i$
  - equal to $Delta E$ of explicit electron addition/removal
],[
  #image(width: 20em, "figures/fig_en_curve_gradients_zoom.svg")
]
))

== Electronic screening via parameters
#slide(repeat: 3, self => [

  #let (uncover, only, alternatives, delayedmark, delayedmarkhl) = methods-with-marks(self)
  $
    E^"KI" [{rho_i}] = &
    E^"DFT" [rho]
    +
    sum_i (
      - integral_0^f_i lr(chevron.l phi_i mid(|) hat(h)^"DFT" (f) mid(|) phi_i chevron.r) dif f
      + f_i integral_0^1 lr(chevron.l phi_i mid(|) hat(h)^"DFT" (f) mid(|) phi_i chevron.r) dif f
    )
    #pause
    \ = & E^"DFT" [rho]
    + sum_i {
      - (E^"DFT" [rho] - delayedmark(#3, E^"DFT" [rho^(f_i arrow.r 0)], tag: #<ENm1_hard>, color: primary))
      + f_i (delayedmark(#3, E^"DFT" [rho^(f_i arrow.r 1)], tag: #<ENp1_hard>, color: primary) - delayedmark(#3, E^"DFT" [rho^(f_i arrow.r 0)], tag: #<ENm1b_hard>, color: primary))
    }
    // uncover("5-", 
    // \ arrow.r E^"uKI" [{rho_i}] approx & 
    // E^"DFT" [rho]
    // \ & +
    // sum_i {
    //   - (E^"DFT" [rho] - delayedmark(#3, E^"DFT" [rho - rho_i], tag: #<ENm1>, color: primary))
    //   + f_i (delayedmark(#3, E^"DFT" [rho - rho_i + n_i], tag: #<ENp1>, color: primary) - delayedmark(#3, E^"DFT" [rho - rho_i], tag: #<ENm1b>, color: primary))
    // }
    // )
  $

  #pause
  #annot(<ENm1_hard>, pos: bottom)[#align(center, [cannot evaluate \ directly])]
  #annot(<ENp1_hard>, pos: bottom)[#align(center, [cannot evaluate \ directly])]
  #annot(<ENm1b_hard>, pos: bottom)[#align(center, [cannot evaluate \ directly])]
  #pause

  // Instead use a frozen-orbital picture:
  
  // $
  //  rho^(f_i arrow.r f)(bold(r)) approx rho(bold(r)) + (f - f_i) |phi^N_i (bold(r))|^2
  // $
  // 
  // very easy to evaluate -- but not at all accurate! Correct this _post hoc_ via a screening parameter i.e.
  // 
  // $
  //   E[rho^(f_i arrow.r f)] approx alpha_i E[rho + (f - f_i) |phi^N_i (bold(r))|^2]
  // $
])

#slide[
#align(center + horizon, 
  image("figures/fig_pwl_DFT.svg", height: 100%)
)
]
#slide[
#align(center + horizon, 
  image("figures/fig_pwl_uKI.svg", height: 100%)
)
]
#slide[
#align(center + horizon, 
  image("figures/fig_pwl_alphaKI.svg", height: 100%)
)
]

#slide(repeat: 5, self => [

  #let (uncover, only, alternatives, delayedmark, delayedmarkhl) = methods-with-marks(self)

$
  E^"KI"_bold(alpha) [rho, {rho_i}] approx & 
  E^"DFT" [rho]
  \ & +
  sum_i delayedmark(#3, alpha_i, tag: #<alpha>, color: primary) {
    - (E^"DFT" [rho] - delayedmark(#2, E^"DFT" [rho - rho_i], tag: #<ENm1>, color: primary))
    + f_i (delayedmark(#2, E^"DFT" [rho - rho_i + n_i], tag: #<ENp1>, color: primary) - delayedmark(#2, E^"DFT" [rho - rho_i], tag: #<ENm1b>, color: primary))
  }
$

#pause
#annot(<ENm1>, pos: bottom)[uses frozen orbitals]
#annot(<ENp1>, pos: bottom)[uses frozen orbitals]
#annot(<ENm1b>, pos: bottom)[uses frozen orbitals]
#pause
#annot(<alpha>, pos: bottom)[screening parameter]
#pause
which is easy to evaluate _e.g._
$ H^"KI"_(i j) = chevron.l phi_j|hat(h)^"DFT" + alpha_i hat(v)_i^"KI"|phi_i chevron.r #h(2cm) hat(v)^"KI"_i = - E_"Hxc" [rho - n_i] + E_"Hxc" [rho] - integral v_"Hxc" (bold(r)', [rho]) n_i d bold(r)' $

#pause
Screening parameters _not_ a fitting parameter!

])

// == Screening
// 
// #grid(columns: (1fr, 2fr), 
// [
// #align(center + horizon, 
//   image("figures/fig_pwl.png", width: 100%)
// )
// ],
// [
// 
// #pause
// Construct $alpha_i$ from explicit $Delta$SCF calculations@Nguyen2018@DeGennaro2022a
// 
// $
//   alpha_i = alpha_i^0 (Delta E_i - lambda_(i i)(0)) / (lambda_(i i)(alpha^0) - lambda_(i i)(0)) "where" lambda_(i i)(alpha) = chevron.l phi_i|hat(h)^"DFT" + alpha hat(v)_i^"KI"|phi_i chevron.r $
// 
// #pause
// Recast via linear response@Colonna2018:
// 
// $
//   alpha_i = (chevron.l n_i mid(|) epsilon^(-1) f_"Hxc" mid(|) n_i chevron.r) / (chevron.l n_i mid(|) f_"Hxc" mid(|) n_i chevron.r)
// $
// 
// which can be efficiently computed via DFPT@Colonna2022
// ],
// )

== Orbital-density dependence
#slide()[
The potential is orbital-density-dependent!
#v(-0.5em)
  $ v^"KI"_(i in"occ") = - E_"Hxc" [rho - n_i] + E_"Hxc" [rho] - integral v_"Hxc" (bold(r)', [rho]) n_i d bold(r)' $

#pause

- loss of unitary invariance@Nguyen2018
#v(-1em)
#align(center,
  grid(columns: (auto, auto), column-gutter: 1em,
  image("figures/fig_nguyen_variational_orbital.png", width: 10em),
  image("figures/fig_nguyen_canonical_orbital.png", width: 10em),
  [two variational orbitals],
  [a canonical orbital],
  )
) #pause
- we can use MLWFs@Marzari2012 #pause
- we know $hat(H)|phi_i chevron.r$ but not $hat(H)$ #pause
- a natural generalisation of DFT towards spectral functional theory@Ferretti2014
]
== Frozen orbital approximation

  #v(-5em)
  #align(center + horizon, 
  grid(align: center + horizon, columns: 3, column-gutter: 2cm, row-gutter: 1cm,
  cetz.canvas({
    import cetz.draw: *
    content((1.25, 1.5), [$rho$])
    circle((0, 0), radius: 1, fill: primary, stroke: none)
    circle((2.5, 0), radius: 1, fill: primary, stroke: none)

  }),
  cetz.canvas({
    import cetz.draw: *

    content((9, 1.5), [$rho^(f_1 arrow.r 0)$])
    arc((10.75, 0), start: 0deg, stop: 360deg, radius: (1.5, 1), fill: primary, stroke: none)
    circle((8, 0), radius: 1, fill: none, stroke: (thickness: 2pt, paint: primary))
    circle((8, 0), radius: 1, fill: none, stroke: (dash: "dashed", thickness: 2pt, paint: white))
    // content((8, -1.5), [$f_1 = 0$])
  }),
  cetz.canvas({
    import cetz.draw: *

    content((17.25, 1.5), [$rho - |psi^N_1(r)|^2$])
    circle((16, 0), radius: 1, fill: none, stroke: (dash: "dashed", thickness: 2pt, paint: primary))
    circle((18.5, 0), radius: 1, fill: primary, stroke: none)
  }),
  [2-electron solution],
  [what we'd like to evaluate],
  [what we can quickly evaluate]

  ))

#matrix-slide(columns: (3fr, 2fr))[
#align(center + horizon,
  {only("1")[#image("figures/alpha_calc/fig_alpha_calc_step_0.png", height: 80%)]
  only("2")[#image("figures/alpha_calc/fig_alpha_calc_step_1.png", height: 80%)]
  only("3")[#image("figures/alpha_calc/fig_alpha_calc_step_2.png", height: 80%)]
  only("4-5")[#image("figures/alpha_calc/fig_alpha_calc_step_3.png", height: 80%)]
  only("6-7")[#image("figures/alpha_calc/fig_alpha_calc_step_4.png", height: 80%)]
  }
)
][
#only("7")[$ alpha_i = alpha_i^0 (Delta E_i - lambda_(i i)(0)) / (lambda_(i i)(alpha^0) - lambda_(i i)(0)) $
$ lambda_(i i)(alpha) = chevron.l phi_i|hat(h)^"DFT" + alpha hat(v)_i^"KI"|phi_i chevron.r $]
]

== Issues with extended systems

#align(center + horizon, 
  image("figures/fig_nguyen_scaling.png", width: 60%)
)

Two options: #pause _1._ use a more advanced functional#pause, or _2._ stay in the "safe" region
#blcite(<Nguyen2018>)

#slide[
=== ZnO @Colonna2022
#v(-1em)
#align(center + horizon,
grid(align: center + horizon, columns: 3, column-gutter: 1em,
image("figures/ZnO_lda_cropped.png", height: 50%),
image("figures/ZnO_hse_cropped_noaxis.png", height: 50%),
image("figures/ZnO_ki_cropped_noaxis.png", height: 50%),
))
#show table.cell: it => {
  set text(size: 0.8em)
  if it.x == 5 {
    set text(fill: primary, weight: "semibold")
    it
  } else {
    it
  }
}
#table(columns: (auto, 1fr, 1fr, 1fr, 1fr, 1fr, 1.5fr), align: center, inset: 0.5em, stroke: none,
table.header([], [LDA ], [HSE ], [GW#sub[0] ], [scGW̃ ], [KI ], [exp ]),
table.hline(),
[$E_"gap"$], [0.79], [2.79], [3.0], [3.2], [3.68], [3.60],
[$chevron.l epsilon_d chevron.r$], [-5.1], [-6.1], [-6.4], [-6.7], [-6.93], [-7.5 to -8.81 ],
[$Delta$], [4.15], [], [], [], [4.99], [5.3]
)
  
]

== Model systems
=== Hooke's atom@Schubert2023

#align(center + horizon, 
  image("figures/schubert_vxc_only.jpeg", height: 70%)
)
= Non-collinear spin
== Non-collinear spin

$ rho_i (bold(r)) pause arrow.r bold(rho)_i (bold(r)) = (rho_i (bold(r)), m_i^x (bold(r)), m_i^y (bold(r)), m_i^z (bold(r))) $

#pause

e.g. for the corrective potential

$ v_i^"qKI" = - 1 / 2 integral dif bold(r) dif bold(r)' rho_i (bold(r)) f_"Hxc" (bold(r), bold(r)') rho_i (bold(r)') + (1 - f_i) integral d bold(r)' f_"Hxc" (bold(r), bold(r)') rho_i (bold(r)') $

#pause

#align(center, sym.arrow.b)

$ v_i^"qKI" = - 1 / 2 integral dif bold(r) dif bold(r)' bold(rho)_i (bold(r)) bb(F)_"Hxc" (bold(r), bold(r)') bold(rho)_i (bold(r)') sigma_0 + (1 - f_i) sum_alpha integral d bold(r)' [bb(F)_"Hxc" (bold(r), bold(r)') bold(rho)_i (bold(r)')]_alpha sigma_alpha $


#blcite(<Marrazzo2024>)

#pagebreak()

CsPbBr#sub[3] #blcite(<Marrazzo2024>)
#v(-2em)
#align(center + horizon,
image("figures/marrazzo_CsPbBr3_bands.svg", height: 50%)
)
#table(align: center, columns: (auto, 1fr, 1fr, 1fr, 1fr, 1fr, 1.5fr), inset: 0.5em, stroke: none,
table.header([], [LDA ], [HSE ], [G#sub[0]W#sub[0] ], [scGW̃ ], [*KI*], [exp ]),
table.hline(),
[*with SOC*], [0.18], [0.78], [0.94], [1.53], [*1.78*], [1.85],
[without SOC], [1.40], [2.09], [2.56], [3.15], [3.12], [],
)

= Caveats

== Limitations

- only valid for systems with $E_"gap"$ > 0 #pause
- empty state localisation in the bulk limit #pause
- can break crystal point group symmetry

== Resonance with other efforts

- Wannier transition state method of Anisimov and Kozhevnikov@Anisimov2005
- Optimally-tuned range-separated hybrid functionals of Kronik, Pasquarello, and others@Kronik2012@Wing2021
- Ensemble DFT of Kraisler and Kronik@Kraisler2013
- Koopmans-Wannier method of Wang and co-workers@Ma2016
- Dielectric-dependent hybrid functionals of Galli and co-workers@Skone2016a
- Scaling corrections of Yang and co-workers@Li2018
= Computational cost and scaling
== Computational cost and scaling
#align(center + horizon,
image("figures/timings/benchmark.svg", width: 80%)
)

#pagebreak()

The vast majority of the computational cost: determining screening parameters

$
  alpha_i = (chevron.l n_i|epsilon^(-1) f_"Hxc"|n_i chevron.r) / (chevron.l n_i|f_"Hxc"|n_i chevron.r)
$

#pause

- a local measure of screening of electronic interactions #pause
- one screening parameter per orbital
- must be computed #emph[ab initio] via... #pause
  - $Delta$SCF@Nguyen2018@DeGennaro2022a: embarrassingly parallel steps which each cost $cal(O)(N_"SC"^3) tilde cal(O)(N_bold(k)^3 N^3)$ #pause
  - DFPT@Colonna2018@Colonna2022: $cal(O)(N_bold(k)^2 N^3)$

== Machine-learned electronic screening
#slide[
  #grid(
    columns: (1fr, 1fr),
    align: center + horizon,
    gutter: 1em,
    image(
      "figures/convergence_key.png",
      height: 5%,
    ) +  v(-1em) +
    image(
      "figures/convergence_fig.png",
      height: 55%,
    ),
    image("figures/speedup.png", height: 60%),

    [*accurate* to within $cal("O")$(10 meV) _cf._ typical band gap accuracy of $cal("O")$(100 meV)],
    [*speedup* of $cal("O")$(10) to $cal("O")$(100)],
  )

  #blcite(<Schubert2024>)
]
== Machine-learned electronic screening

#pagebreak()

#slide[
  #align(
    center,
    grid(
      columns: 5,
      align: horizon,
      gutter: 1em,
      image("figures/orbital.emp.00191_cropped.png", height: 30%),
      $stretch(->)^("power spectrum decomposition")$,
      $vec(delim: "[", x_0, x_1, x_2, dots.v)$,
      $stretch(->)^("ridge regression")$,
      $alpha_i$,
    ),
  )

  $
    c^i_(n l m, k) & = integral dif bold(r) g_(n l) (r) Y_(l m)(theta,phi) n^i (
      bold(r) - bold(R)^i
    )
  $


  $
    p^i_(n_1 n_2 l,k_1 k_2) = pi sqrt(8 / (2l+1)) sum_m c_(n_1 l m,k_1)^(i *) c_(n_2 l m,k_2)^i
  $

  #blcite(<Schubert2024>)
]

#pagebreak()

#slide[
  #align(
    center,
    grid(
      columns: 2,
      align: horizon + center,
      gutter: 1em,
      image("figures/water.png", height: 70%),
      image("figures/CsSnI3_disordered.png", height: 70%),

      "water", "CsSnI" + sub("3"),
    ),
  )
  #blcite(<Schubert2024>)
]

The use-case

   #grid(columns: 8, column-gutter: 0.3em, row-gutter: 0.3em,
        image("figures/CsSnI3_disordered.png", width: 100%),
        image("figures/CsSnI3_disordered.png", width: 100%),
        image("figures/CsSnI3_disordered.png", width: 100%),
        image("figures/CsSnI3_disordered.png", width: 100%),
        image("figures/CsSnI3_disordered.png", width: 100%),
        image("figures/CsSnI3_disordered.png", width: 100%),
        image("figures/CsSnI3_disordered.png", width: 100%),
        grid.cell(align: center + horizon, [...]),
        grid.cell(inset: 0.4em, align: center, fill: primary, colspan: 3, text(fill: white, "train", size: 1em, weight: "bold")),
        grid.cell(inset: 0.4em, align: center, fill: secondary, colspan: 5, text("predict", size: 1em, weight: "bold")),
  )

  #pause
  N.B. not a general model


#slide[
  #grid(
    columns: (1fr, 1fr),
    align: horizon + center,
    gutter: 1em,
    image(
      "figures/water_cls_calc_vs_pred_and_hist_bottom_panel_alphas.svg",
      height: 70%,
    ),
    image(
      "figures/CsSnI3_calc_vs_pred_and_hist_bottom_panel_alphas.svg",
      height: 70%,
    ),

    "water", "CsSnI" + sub("3"),
  )
  #blcite(<Schubert2024>)
]

#slide[
  #grid(
    columns: (1fr, 1fr),
    align: center + horizon,
    gutter: 1em,
    image(
      "figures/convergence_key.png",
      height: 5%,
    ) +  v(-1em) +
    image(
      "figures/convergence_fig.png",
      height: 55%,
    ),
    image("figures/speedup.png", height: 60%),

    [*accurate* to within $cal("O")$(10 meV) _cf._ typical band gap accuracy of $cal("O")$(100 meV)],
    [*speedup* of $cal("O")$(10) to $cal("O")$(100)],
  )

  #blcite(<Schubert2024>)
]


== Taking advantage of symmetries
To compute screening parameters via DFPT...
#algorithm(inset: 0.3em, indent: 1em, {
  import algorithmic: *
  Function("CalculateAlpha", ($n$,), {
    For($bold(q) in "BZ"$,
    {
        For($bold(k) in "BZ"$, {Comment[Linear system $A x = b$ to obtain $Delta psi_(bold(k)+bold(q),v)(bold(r))$]})
          Assign[$Delta rho^(0n)_(q)$][$sum_(bold(k)v)psi^*_(bold(k)v) (bold(r))Delta psi_(bold(k)+bold(q),v)(bold(r)) + c.c.$]
          Assign[$Pi^((r))_(0 n, bold(q))$][$chevron.l Delta rho^(0 n)_(bold(q))|f_"Hxc"|rho^(0 n)_(bold(q)) chevron.r$]
          Assign[$Pi^((u))_(0 n, bold(q))$][$chevron.l rho^(0 n)_bold(q)|f_"Hxc"|rho^(0 n)_bold(q) chevron.r$]
    })
    Return[$1 + sum_bold(q) Pi^((r))_(0 n, bold(q)) \/ sum_bold(q) Pi^((u))_(0 n, bold(q))$]
  })
})

#pagebreak()

#align(center,
  image("figures/bz-to-ibz-outer.svg", height: 80%)
)
$bold(q) in "BZ" $ $arrow.r$ $bold(q) in "IBZ"(n)$ (the symmetry of the perturbation; lower than that of the primitive cell)
#pagebreak()
#align(center,
  image("figures/bz-to-ibz-inner.svg", height: 80%)
)
$bold(k) in "BZ"$ $arrow.r$ $bold(k) in "IBZ"(bold(q))$ (can only use symmetries that leave $bold(q)$ invariant)

#align(horizon + center, image("figures/bz-to-ibz-speedup.svg", height: 100%))

= Automated Wannierisation
== 
#slide()[
#set text(size: 0.8em)
#raw(read("scripts/gaas.json"), block: true, lang: "json")
]
== Automated Wannierisation
#slide()[
  Koopmans functionals rely heavily on Wannier functions...
  - to initialise the minmising orbitals, _or_
  - in place of the minimising orbitals entirely

#pause

#grid(
  columns: (2fr, 2fr, 3fr),
  align: center + horizon,
  gutter: 1em,
  image("figures/proj_disentanglement_fig1a.png", height: 45%),
  image("figures/new_projs.png", height: 45%),
  image("figures/target_manifolds_fig1b.png", height: 45%),

  text("projectability-based disentanglement") + cite(<Qiao2023>),
  text("use PAOs found in pseudopotentials"),
  text("parallel transport to separate manifolds") + cite(<Qiao2023a>),
)
]

==
#align(center + horizon,
image("figures/supercell_workflow.png", width: 100%)
)


== 
#blcite(<Huber2020>)
#v(-2em)
#align(center,
  [
  #grid(columns: 3, align: horizon, column-gutter: 0.5em,
    image("media/logos/koopmans_grey_on_transparent.svg", height: 3em),
    image("figures/handshake.png", height: 2em, alt: "handshake"),
    image("media/logos/aiida.svg", height: 3em)
  )
  #pause `$ koopmans run tio2.json` #pause $arrow.r$ `$ koopmans run --engine=aiida tio2.json`
  ]
)

remote compute, parallel step execution, provenance-tracking, (requires configuration, WIP...)

#pause
#align(center, 
  image("figures/aiida-speed-up.svg", width: 70%)
)

== Connections with approx. self-energies

#blcite(<Ferretti2014>)#blcite(<Colonna2019>)

Orbital-density functional theory:

$ (h + alpha_i v^(K I)_i)|psi_i chevron.r = lambda_i|psi_i chevron.r $ $v_i^(K I)(bold(r))$ is real, local, and state-dependent #pause

cf. Green's function theory:

$ (h + Sigma_i)|psi_i chevron.r = z_i|psi_i chevron.r $ $Sigma_i (bold(r), bold(r)')$ is complex, non-local, and state-dependent

#slide[
Hartree-Fock self-energy in localized representation

$Sigma_x (bold(r), bold(r)') = - sum_(k sigma)^("occ") psi_(k sigma)(bold(r)) & f_H (bold(r), bold(r'))psi^*_(k sigma)(bold(r)') \
& arrow.r.double.long chevron.l phi_(i sigma)|Sigma_x|phi_(j sigma') chevron.r approx - chevron.l phi_(i sigma)|v_H [n_(i sigma)]|phi_(i sigma)chevron.r delta_(i j)delta_(sigma sigma')$

Unscreened KIPZ#sym.at Hartree ($v_"xc" arrow.r 0$; $f_"Hxc" arrow.r f_H$; $epsilon^(-1) arrow.r 1$)

$chevron.l phi_(i sigma)|v^"KIPZ"_(j sigma',"xc")|phi_(j sigma') chevron.r
approx {(1/2 - f_(i sigma)) chevron.l n_(i sigma)|f_H|n_(i sigma) chevron.r - E_H [n_(i sigma)]}
approx - chevron.l phi_(i sigma)|v_H [n_(i sigma)]|phi_(i sigma)chevron.r delta_(i j)delta_(sigma sigma')$

]

#slide[
Screened exchange plus Coulomb hole (COHSEX)

$ Sigma^"SEX"_"xc" (bold(s), bold(s)') = - sum_(k sigma)^"occ" psi_(k sigma)(bold(r)) psi_(k sigma)^*(bold(r)) W(bold(r), bold(r)') $

$ Sigma^"COH"_"xc" (bold(s), bold(s)') = 1/2 delta(bold(s), bold(s)'){W(bold(r), bold(r)') - f_H (bold(r), bold(r)')} $

$ arrow.r.double.long chevron.l phi_(i sigma)|Sigma^"COHSEX"_"xc"|phi_(j sigma')chevron.r approx {(1/2 - f_(i sigma)) chevron.l n_(i sigma)|W|n_(i sigma)chevron.r - E_H [n_(i sigma)]}delta_(i j) delta_(sigma sigma')$

KIPZ#sym.at Hartree with RPA screening ($v_"xc" arrow.r 0$; $f_"Hxc" arrow.r f_H$; $epsilon^(-1) arrow.r "RPA"$)

$ chevron.l phi_(i sigma)|v^"KIPZ"_(j sigma',"xc")|phi_(j sigma')chevron.r approx{(1/2 - f_(i sigma)) chevron.l n_(i sigma)|W|n_(i sigma)chevron.r - E_H [n_(i sigma)]}delta_(i j) delta_(sigma sigma')$
]

#slide[
  Static GWΓ#sub[xc] --- local (DFT-based) vertex corrections@Hybertsen1987@DelSole1994

  $ Sigma^(G W Gamma_"xc")_"xc"(1, 2) = i G(1, 2) W_(t-e) (1, 2) $
  
  $ W_(t-e) = (1 - f_"Hxc" chi_0)^(-1) f_H $

  $ arrow.r.double.long chevron.l phi_(i sigma)|Sigma^(G W Gamma_"xc")_"xc"|phi_(j sigma')chevron.r approx{(1/2 - f_(i sigma)) chevron.l n_(i sigma)|W_(t-e)|n_(i sigma)chevron.r - E_H [n_(i sigma)]}delta_(i j) delta_(sigma sigma')$

  KIPZ#sym.at DFT ($v_"xc" arrow.r$ DFT; $f_"Hxc" arrow.r$ DFT; $epsilon^(-1) arrow.r$ DFT)

  $ chevron.l phi_(i sigma)|v^"KIPZ"_(j sigma',"xc")|phi_(j sigma')chevron.r approx{chevron.l phi_(i sigma)|v^"DFT"_(sigma,"xc")|phi_(i sigma)chevron.r + (1/2 - f_(i sigma)) chevron.l n_(i sigma)|epsilon^(-1)_(t-e) f_"Hxc"|n_(i sigma)chevron.r - E_H [n_(i sigma)]}delta_(i j) delta_(sigma sigma')$
]

== Open questions
- why does correcting _local_ charged excitations correct the description of delocalized excitations?
- is there a good metric for selecting variational orbitals (_i.e._ the subspace with respect to which we enforce piecewise linearity)?
- are off-diagonal corrections appropriate? What form should they take?
- how to extend to metallic systems?
- can we provide a formal basis for the Koopmans correction?
  - GKS
  - spectral functional theory@Ferretti2014
  - ensemble DFT
  - RDMFT

== Deviation from the flat plane condition

#align(center + horizon,
  image("figures/he_deviation_from_flat_plane.svg", width: 50%)
)

== What is screening $U$?

#grid(columns: (1fr, 1fr), align: horizon + center, inset: 0.5em,
cetz.canvas({
  import cetz.draw: *

  let counter = 0
  let positions = ((0, 0), (2, 5), (3, 2), (4, -1), (6, 2), (5, 5), (-1, 3),)
  for pos in positions {
    counter = counter + 1
    for i in range(1, counter) {
      line(name: "line", positions.at(i - 1), pos, stroke: gray + 1pt)
      // content("line", text(fill: gray, [$chi_(#i #counter)$]), midpoint: 0.5)
    }
  }
  counter = 0
  for pos in positions {
    counter = counter + 1
    circle(pos, radius: 0.5, fill: primary, stroke: none, name: "circle")
    content("circle", text(fill: white, [#counter]))
  }
}),
cetz.canvas({
  import cetz.draw: *

  rect((-4, -3), (4, 3), stroke: none, fill: secondary, alpha: 0.5)
  content((4, 3), [bath], anchor: "north-east", padding: 0.5em)
  circle((0, 0), radius: 0.5, fill: primary, stroke: none)
  
}),
  [all sites included in response matrix],
  [only one site included in response matrix],
  [bare $U$],
  [fully-screened $U$],
)

#slide(self=> [

  #let (uncover, only, alternatives, delayedmark, delayedmarkhl) = methods-with-marks(self)
  #table(columns: (auto, auto, 1fr), align: horizon + center, inset: 1em,
    [
      fully-screened],
    [
      #set text(size: 0.6em)
      $
        mat(mat(delim: #none, chi^(arrow.t arrow.t)_(1 1), ; , chi^(arrow.b arrow.b)_(1 1)),,;,
            mat(delim: #none, chi^(arrow.t arrow.t)_(2 2), ; , chi^(arrow.b arrow.b)_(2 2)),;,,,dots.down)
      $ 
    ],
    [
      $
        U^(I sigma) = 1 / (chi_0)^(sigma sigma)_(I I) - 1 / chi^(sigma sigma)_(I I)
      $
    ],

    [
      not screened by \ opposite spin],
    [
      #set text(size: 0.7em)
      $
      mat(mat(delim: #none, chi^(arrow.t arrow.t)_(1 1), chi^(arrow.t arrow.b)_(1 1); chi^(arrow.b arrow.t)_(1 1), chi^(arrow.b arrow.b)_(1 1)),,;,
          mat(delim: #none, chi^(arrow.t arrow.t)_(2 2), chi^(arrow.t arrow.b)_(2 2); chi^(arrow.b arrow.t)_(2 2), chi^(arrow.b arrow.b)_(2 2)),;,,dots.down)
      $ 
    ],
    [
      $
      f^(sigma sigma')_I = [(chi_0)^(sigma sigma)_(I I)]^(-1) - [chi^(sigma sigma)_(I I)]^(-1) \ f^(sigma sigma')_(I) stretch(arrow.r)^(???) U^I "or" U^(I sigma)
      $
    ],
    [
      also not screened by \ other Hubbard sites],
    [
      #set text(size: 0.7em)
    $
    mat(mat(delim: #none, chi^(arrow.t arrow.t)_(1 1), chi^(arrow.t arrow.b)_(1 1); chi^(arrow.b arrow.t)_(1 1), chi^(arrow.b arrow.b)_(1 1)),
        mat(delim: #none, chi^(arrow.t arrow.t)_(1 2), chi^(arrow.t arrow.b)_(1 2), dots; chi^(arrow.b arrow.t)_(1 2), chi^(arrow.b arrow.b)_(1 2), dots);
        mat(delim: #none, chi^(arrow.t arrow.t)_(2 1), chi^(arrow.t arrow.b)_(2 1); chi^(arrow.b arrow.t)_(2 1), chi^(arrow.b arrow.b)_(2 1); dots.v, dots.v), 
        mat(delim: #none, chi^(arrow.t arrow.t)_(2 2), chi^(arrow.t arrow.b)_(2 2), dots; chi^(arrow.b arrow.t)_(2 2), chi^(arrow.b arrow.b)_(2 2), dots; dots.v, dots.v, dots.down))
    $ 
    ],
    [
      $ f^(sigma sigma')_(I J) = ... $ (left as an exercise to the reader)
    ]
  )
])

== Advantages of spin-resolved linear response
#focus-slide()[1. Conceptual consistency]
== 

#align(center + horizon, [spin-resolved linear response #sym.arrow.l.r spin-resolved DFT+_U_ functional])

#v(3em)
#pause
#align(right, [... we didn't explore DFT+$U^sigma$; instead see BLOR@Burgess2023)])

#focus-slide()[2. Unconstrained constrained linear response]
== Unconstrained constrained linear response

Suppose we want to compute $ lr((d^2E_"Hxc") / (d (n^I)^2) |)_(mu^I)$ #pause

This is easy with spin-resolved LR: #pause

$
  (d^2 E_"Hxc") / (d (n^I)^2) = & 1 / 2 (d v_"Hxc"^arrow.t + d v_"Hxc"^arrow.b) / d(n^arrow.t + n^arrow.b)
  = 1 / 2 (f^(arrow.t arrow.t ) d n^arrow.t + f^(arrow.t arrow.b) d n^arrow.b + f^(arrow.b arrow.t) d n^arrow.t + f^(arrow.b arrow.b) d n^arrow.b) / (d n^arrow.t + d n^arrow.b)
$

#pause
"Impose" the constraint by setting $d n^arrow.t = d n^arrow.b$ to get...

$
  lr((d^2E_"Hxc") / (d n^2) |)_(mu) = & 1/4 (f^(arrow.t arrow.t) + f^(arrow.b arrow.b) + f^(arrow.t arrow.b) + f^(arrow.b arrow.t))
$

#pause
This simple average is one choice (of many) for $M: f^(sigma sigma')_I arrow.r U^I$

#focus-slide()[3. We can recover conventional linear response]

== Conventional linear response

For conventional LR, $ d v^(I arrow.t) = d v^(I arrow.b) = d v$ #pause, in which case:

$
d n = sum_sigma d n^(sigma) = sum_(sigma sigma') chi^(sigma sigma') d v^(sigma') = sum_(sigma sigma') chi^(sigma sigma') d v
arrow.double.r.long chi_"conv" = (d n) / (d v) = sum_(sigma sigma') chi^(sigma sigma')
$
#pause
Likewise,
$
  (epsilon^(-1))_"conv" = ... = 1 / 2 sum_(sigma sigma') (f chi)^(sigma sigma')
$
#pause
And thus

$
  U = (epsilon^(-1) - 1) chi^(-1) = 1 / 2 (sum_(sigma sigma') (f chi)^(sigma sigma')) / (sum_(sigma sigma') chi^(sigma sigma'))
$
#focus-slide()[4. $J$ is free]

==
As defined by 

$
  J = - 1 / 2 (d v_"Hxc"^arrow.t - d v_"Hxc"^arrow.b) / (d (n^arrow.t - n^arrow.b)) = - 1 / 4 ((f^(arrow.t arrow.t) - f^(arrow.b arrow.t)) d n^arrow.t - (f^(arrow.b arrow.b) - f^(arrow.t arrow.b)) d n^arrow.b)/ (d (n^arrow.t - n^arrow.b))
$

Different ways to define $J$: #pause
+ while keeping $n = n^arrow.t + n^arrow.b$ fixed:
  $
    J = - 1 / 4 (f^(arrow.t arrow.t) - f^(arrow.b arrow.t) - f^(arrow.t arrow.b) + f^(arrow.b arrow.b))
  $
  #pause
+ for a perturbation where $d v^arrow.t = - d v^arrow.b$

#focus-slide()[5. Easy to implement]


= References
== References
#bibliography("references.bib", style: "nature-footnote.csl", title: none)