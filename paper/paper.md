---
title: 'eagar-tsai: A Python library for the Eagar-Tsai moving heat source model'
tags:
  - Python
  - additive manufacturing
  - laser powder bed fusion
  - melt pool
  - heat transfer
authors:
  - name: Doğuhan Sarıtürk
    orcid: 0000-0003-0008-1165
    corresponding: true
    affiliation: 1
  - name: Brent Vela
    orcid: 0000-0001-6465-3738
    affiliation: 1
  - name: Kyle Swartz
    orcid: 0009-0001-5618-4203
    affiliation: 1
  - name: Raymundo Arróyave
    orcid: 0000-0001-7548-8686
    affiliation: "1, 2, 3"
affiliations:
  - index: 1
    name: "Department of Materials Science and Engineering, Texas A&M University, USA"
    ror: 01f5ytq51
  - index: 2
    name: "J. Mike Walker '66 Department of Mechanical Engineering, Texas A&M University, USA"
    ror: 01f5ytq51
  - index: 3
    name: "Wm Michael Barnes '64 Department of Industrial and Systems Engineering, Texas A&M University, USA"
    ror: 01f5ytq51
date: 11 September 2026
bibliography: paper.bib
---

# Summary

`eagar-tsai` is a Python library for predicting melt pool dimensions and temperatures during laser processing of metals using the analytical moving heat source model of @eagar1983. The model takes laser power, scan velocity, beam diameter, and absorptivity as inputs, together with the material’s liquidus temperature, thermal conductivity, density, and specific heat. The outputs are the melt pool length, width, and depth, along with the peak temperature. With a compiled C extension, the library evaluates each condition in milliseconds and can process them individually or in batches across multiple CPU cores. It also generates printability maps across a grid of laser power and scan velocity, classifying each condition as keyhole, lack-of-fusion, balling, or defect-free. Built-in plotting functions visualize temperature fields as two-dimensional cross-sections or three-dimensional volume renders, and display printability map results as color-coded diagrams.

# Statement of need

Laser powder bed fusion, directed energy deposition, and laser welding all depend on melt pool geometry, because pool dimensions set the solidification conditions that govern grain structure and the mechanical properties of the deposited material. High-fidelity finite element and computational fluid dynamics simulations resolve pool geometry accurately, but a single condition takes minutes to hours of compute time. That cost puts full-physics methods out of reach for broad parameter sweeps, uncertainty quantification, and the generation of training data for machine learning models.

Analytical heat-source models are a practical alternative. The Rosenthal point source solution [@rosenthal1946] is the simplest of these, but its singularity at the source produces unrealistically high peak temperatures, and it underestimates pool width at finite beam diameters. The Eagar-Tsai model [@eagar1983] replaces the point source with a two-dimensional Gaussian distribution, which removes the singularity and improves agreement with measured pool dimensions. The model was developed for welding and is now applied widely in metal additive manufacturing [@honarmandi2021; @whalen2021; @menon2024; @rangaswamy2024], yet no maintained Python implementation with a batch-capable API has been publicly available. `eagar-tsai` is a maintained implementation of the model, with a batch API that evaluates thousands of conditions in a single call.

# State of the field

Several categories of tools exist for predicting thermal fields in laser processing, each occupying a different point in the accuracy-cost trade-off space. Finite element solvers and computational fluid dynamics codes such as FEniCS [@dolfinx2023], Firedrake [@firedrake2023], and OpenFOAM-based frameworks including AdditiveFOAM [@additiveFOAM] predict melt pool geometry well, because they model melt convection, surface depression, and phase change directly. Resolving that physics requires a transient solve on a fine mesh for every condition, which restricts them to a handful of conditions per study.

Analytical models occupy the other end of the trade-off space. The closest public implementation of the Eagar-Tsai model is an application rather than a library. The printability-map framework of @sheikh2024 maps alloy composition to thermophysical properties, evaluates analytical, scaled, and neural-network forms of the model, and applies twelve defect criteria across lack-of-fusion, keyholing, and balling. The model is not separable from the pipeline. A user configures the scripts by editing variables in the source, runs them from the command line, and collects results from CSV files [@printability_map_repo]. Applying the model to a different question means going back into the source. Other public implementations take the same form, released as artifacts attached to individual studies. `eagar-tsai` is instead a library, and its conditions are function arguments rather than edits to a script.

# Model description

The model assumes a semi-infinite solid with a planar top surface, constant material properties evaluated at the liquidus temperature, a Gaussian heat source with constant absorptivity, and a steady-state temperature field that moves rigidly with the beam. It omits latent heat of fusion, melt convection, and vaporization. Published assessments report melt pool widths within ten percent of measured values, with larger errors in depth under keyhole conditions [@honarmandi2021; @rangaswamy2024].

The resulting temperature at position $(x, y, z)$ in the co-moving frame follows from the solution of @eagar1983, written here in the non-dimensional form of @rubenchik2018 that the library implements:

$$T(x, y, z) = T_0 + T_s \int_0^\infty f(t;\, x, y, z, p)\, \mathrm{d}t$$

where $T_0 = 298\,\mathrm{K}$ is the ambient temperature, $t$ is a dimensionless integration variable, and $T_s$ is a dimensional prefactor defined by

$$T_s = \frac{A\,P}{\pi\,\rho\,c_p\,\sqrt{\pi\,\alpha\,v\,\sigma^3}}$$

with $A$ the absorptivity, $P$ the laser power, $v$ the scan velocity, $\rho$ the density, $c_p$ the specific heat, $k$ the thermal conductivity, and $\alpha = k/(\rho c_p)$ the thermal diffusivity. The beam width is $\sigma = \sqrt{2}\,d/2$ for a beam of diameter $d$. It is $\sqrt{2}$ times the standard deviation of the Gaussian source, which is the quantity @eagar1983 use as their distribution parameter. The integrand is:

$$f(t;\, x, y, z, p) = \frac{1}{(4pt + 1)\sqrt{t}} \exp\!\left(-\frac{z^2}{4t} - \frac{y^2 + (x - t)^2}{4pt + 1}\right)$$

The integrand takes non-dimensional coordinates, with $x$ and $y$ scaled by $\sigma$ and $z$ by $\sqrt{\alpha\sigma/v}$. The remaining parameter is $p = \alpha/(v\sigma)$. Heat transport is diffusion-dominated at large $p$ and advection-dominated at small $p$.

![A representative `eagar-tsai` workflow and output. Panels **a** and **b** show the top-surface and centerline temperature fields, respectively, for a steel-like benchmark condition ($P = 200\,\mathrm{W}$, $v = 0.5\,\mathrm{m\,s^{-1}}$, $d = 100\,\mu\mathrm{m}$, $A = 0.35$), with the liquidus contour overlaid. Panel **c** shows the three-dimensional temperature volume for the same benchmark condition. Panel **d** shows parameter maps for melt-pool length, width, and depth across a regular power-velocity sweep at fixed beam diameter and material properties. Panel **e** shows the printability map classifying each power-speed combination as keyhole, lack-of-fusion, balling, or defect-free. Panel **f** compares batch-runtime scaling with increasing numbers of conditions using the compiled C integrand and the pure-Python fallback.\label{fig:figure1}](figure_1.pdf)

The temperature field is evaluated on a discrete spatial grid covering two planes: the top surface ($z = 0$, \autoref{fig:figure1}**a**) for melt pool length and width, and the beam centerline ($y = 0$, \autoref{fig:figure1}**b**) for depth. The solver can also evaluate the field on a full three-dimensional grid (\autoref{fig:figure1}**c**). The melt pool boundary is identified as the liquidus isotherm, $T = T_\mathrm{liq}$. When the melt pool reaches the trailing $x$, outer $y$, or bottom $z$ boundary, the solver enlarges the domain in that direction and repeats the calculation, up to twenty times. The leading $x$ boundary stays fixed.

# Software design

The most computationally intensive step is the numerical evaluation of the integral at each grid point. SciPy's `quad` function [@virtanen2020] handles each evaluation using adaptive quadrature. To reduce overhead, the integrand is also implemented in C and called directly by the integration routine, bypassing the Python interpreter on each function evaluation. Benchmarks show the C path is approximately ten times faster than the pure-Python equivalent (\autoref{fig:figure1}**f**). The package retains a pure-Python fallback that runs whenever the compiled extension cannot be imported. Pre-built wheels ship the extension, and most installations never use the fallback.

The single-condition API takes laser and material parameters as typed, immutable dataclasses, which parallel runs can distribute to worker processes without risk of mutation. Values derived from them, such as the beam width and thermal diffusivity, are computed outside the integration loops. The batch API accepts conditions in a pandas DataFrame and can distribute them across CPU cores. A power-velocity sweep run through this API yields parameter maps for melt pool length, width, and depth (\autoref{fig:figure1}**d**).

The library also supports optional plotting. It draws temperature fields on the top surface and the beam centerline, with the liquidus isotherm and melt pool dimensions marked (\autoref{fig:figure1}**a** and \autoref{fig:figure1}**b**), or renders them as a three-dimensional volume (\autoref{fig:figure1}**c**). The volume can also be exported as a VTK ImageData file for external viewers. A dedicated entry point automates printability map generation and plotting, which labels each point on a power-velocity grid by defect regime with the criteria of @sheikh2024 (\autoref{fig:figure1}**e**).

# Research impact statement

The package is intended for research settings where analytical melt-pool estimates are needed at the scale of parameter sweeps (\autoref{fig:figure1}**d**), screening studies, and workflow integration rather than for one-off calculations. Its typed API and batch-processing interface make it practical to embed the Eagar-Tsai model in calibration loops, process-window studies, and data-generation pipelines without each group maintaining its own private implementation.

The Eagar-Tsai model underpins a continuing line of alloy and process design work in the authors' group, spanning calibration of the model against single-track experiments [@honarmandi2021], automated printability maps [@sheikh2024], Bayesian data augmentation of melt pool dimensions [@morcos2024], high-throughput alloy and process design [@sheikh2025npj], and printability assessment for high entropy alloys [@sheikh2025hea]. `eagar-tsai` packages the implementation that work relies on, so later studies start from a tested, installable library rather than a private copy of the solver.

# AI usage disclosure

Claude Code (Anthropic, Claude Opus 5) and Codex (OpenAI, GPT-5) were used in the development of `eagar-tsai` to draft docstrings and documentation, scaffold unit tests, and suggest refactorings of existing code, and to copy-edit this manuscript. The authors reviewed, edited, and validated all AI-assisted output before including it, accepting suggested code changes only after the full test suite passed, and are responsible for the correctness of the software, its documentation, and this paper.

# Acknowledgements

The C integrand reformulation implemented in this library was developed by Sasha Rubenchik at Lawrence Livermore National Laboratory (2015). R.A. and D.S. acknowledge support from the U.S. Army Research Office (ARO) through Grant No. W911NF-22-2-0117 and the U.S. Army Contract No. W911NF-25-1-0112. Portions of this research were conducted with the advanced computing resources provided by Texas A&M High Performance Research Computing.

# References
