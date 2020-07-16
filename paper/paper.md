---
title: 'NetworkChange: Analyzing Network Changes in R'
tags:
  - R
  - longitudinal networks
  - latent space 
  - hidden Markov models
  - Bayesian inference
authors:
- affiliation: 1
  name: Jong Hee Park
  orcid: 0000-0002-4598-3456
- affiliation: 2
  name: Yunkyu Sohn
  orcid: 0000-0002-9496-1907
affiliations:
- index: 1
  name: Department of Political Science and International Relations, Seoul National University
- index: 2
  name: Faculty of Political Science and Economics, Waseda University
date: 9 July 2020
bibliography: paper.bib
---
\newcommand{\A}{\mathcal{A}}
\newcommand{\B}{\mathcal{B}}
\newcommand{\Y}{\mathbf{Y}}
\newcommand{\Z}{\mathbf{Z}}
\renewcommand{\E}{\mathbf{E}}
\newcommand{\M}{\mathbf{M}}
\newcommand{\U}{\mathbf{U}}
\newcommand{\V}{\mathbf{V}}
\renewcommand{\L}{\mathbf{L}}
\renewcommand{\u}{\mathbf{u}}
\renewcommand{\v}{\mathbf{v}}
\renewcommand{\l}{\boldsymbol{\lambda}}
\newcommand{\z}{\mathbf{z}}
\newcommand{\y}{\mathbf{y}}
\newcommand{\x}{\mathbf{x}}
\newcommand{\X}{\mathbf{X}}
\newcommand{\W}{\mathbf{W}}
\newcommand{\F}{\mathbf{F}}
\newcommand{\G}{\mathbf{G}}
\newcommand{\R}{\mathbf{R}}
\def\real{\mathbb{R}}
\newcommand{\s}{\mathbf{s}}
\renewcommand{\S}{\mathbf{S}}
\renewcommand{\P}{\mathbf{P}}
\newcommand{\Sig}{\boldsymbol{\Sigma}}
\renewcommand{\a}{\alpha}
\renewcommand{\b}{\boldsymbol{\beta}}
\renewcommand{\t}{\boldsymbol{\theta}}
\newcommand{\T}{\boldsymbol{\Theta}}
\newcommand{\mb}{\mathbf}
\newcommand{\RD}{\mathbf{R}^{D}}
\newcommand{\e}{\boldsymbol{\varepsilon}}
\renewcommand{\r}{\rho}
\newcommand{\g}{\boldsymbol{\gamma}}
\renewcommand{\d}{\boldsymbol{\delta}}
\newcommand{\bs}{\boldsymbol}
\newcommand{\normdist}[2]{\ensuremath{\mathcal{N}(#1,#2)}}
\newcommand{\normdistk}[3]{\ensuremath{\mathcal{N}_{#3}(#1,#2)}}
\newcommand{\wish}[2]{\ensuremath{\mathcal{W}(#1,#2)}}
\newcommand{\invwish}[2]{\ensuremath{\mathcal{IW}(#1,#2)}}
\newcommand{\gamdist}[2]{\ensuremath{\mathcal{G}(#1,#2)}}
\newcommand{\invgam}[2]{\ensuremath{\mathcal{IG}(#1,#2)}}
\newcommand{\studt}[3]{\ensuremath{t_{#3}(#1,#2)}}
\newcommand{\binomial}[2]{\ensuremath{\mathcal{B}in(#1,#2)}}
\newcommand{\bern}[1]{\ensuremath{\mathcal{B}ernoulli(#1)}}
\newcommand{\diri}[1]{\ensuremath{\mathcal{D}irichlet(#1)}}
\newcommand{\unif}[2]{\ensuremath{\mathcal{U}(#1,#2)}}
\newcommand{\chisqr}[1]{\ensuremath{\chi_{#1}^{2}}}
\newcommand{\invchisqr}[1]{\ensuremath{\mathcal{I}nv}\textnormal{-}\ensuremath{\chi_{#1}^{2}}}
\newcommand{\betadist}[2]{\ensuremath{\mathcal{B}eta(#1,#2)}}
\newcommand{\poisson}[1]{\ensuremath{\mathcal{P}oisson(#1)}}
\newcommand{\expo}[1]{\ensuremath{\mathcal{E}xp(#1)}}
\newcommand{\Dir}{\mathrm{Dir}}
\newcommand{\thh}{^\mathrm{th}}
\newcommand{\modtwo}{\mathrm{[mod~2]}}
\newcommand{\thetaof}[2]{\theta \langle #1;#2\rangle}
\newcommand{\Mpa}{M_\mathrm{P,A}}
\newcommand{\Ma}{M_\mathrm{A}}
\newcommand{\rjaccept}{\mathcal{A}}

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
require(NetworkChange)
require(sna)
knitr::opts_chunk$set(
  dpi=300,fig.width=7,
  warning=FALSE,
  message=FALSE,
  collapse = TRUE,
  ## fig.asp = 1, 
  comment = "#>"## , 
  ## eval=FALSE
)
```

````{r img-setup, include=FALSE, cache=FALSE}
out.format <- knitr::opts_knit$get("out.format")
img_template <- switch( out.format,
                     word = list("img-params"=list(fig.width=6,
                                                   fig.height=6,
                                                   dpi=150)),
                     {
                       # default
                       list("img-params"=list( dpi=150,
                                               fig.width=6,
                                               fig.height=6,
                                               out.width="504px",
                                               out.height="504px"))
                     } )

knitr::opts_template$set( img_template )
````

# Summary

**NetworkChange** is an R package that detects multiple structural
changes in longitudinal network data using the latent space approach [@hoff2002latent].
Based on the Bayesian multi-array representation of longitudinal
networks [@Hoff2015], **NetworkChange** performs
Bayesian hidden Markov analysis to discover changes in structural
network features across temporal layers using the hidden Markov model formulation. **NetworkChange** can detect various forms of changes in network structure such as block-splitting, block-merging,
and core-periphery changes. **NetworkChange** also provides functions
for model diagnostics using WAIC, average loss, and log marginal
likelihoods as well as visualization tools for dynamic analysis results
of longitudinal networks. 

#  Statement of Need

The package is designed for R users who need to analyze longitudinal network data to discover latent node-level characteristics including cases when there are discrete changes of the underlying states governing the node-level characteristics. This is in contrast to an R package for latent space and cluster analysis of networks [@krivitsky2008fitting] which does not incorporate a state space model (e.g. hidden Markov model) and a Python code for longitudinal network analysis [@peel2015detecting] under a distinct formulation (hierarchical random graph model) with a changepoint detection function. In addition to functions for the statistical analysis, **NetworkChange** provides visualization functions for summary of the analysis results (\autoref{fig:list}). The complete guide for using core functions of the package is presented at https://github.com/jongheepark/NetworkChange as its vignette with an empirical data set analysis example. @ParkSohn2020 provide methodological details of the algorithms implemented in the package. 

#  Empirical Data Analysis Example

In this section, we analyze changes in the international military alliance network among major powers. The data set is originally from [@Gibler2009] and users can call this data set by `data(MajorAlly)`.

Our goal in this section is to detect structural changes in the longitudinal alliance network among major powers using HNC.  We follow the COW dataset's coding of "major powers" (the United Kingdom,  Germany, Austria-Hungary, France, Italy, Russia, the United States, Japan, and China) in the analysis. We aggregated every 2 year network from the original annual binary networks to increase the density of each layer.


```{r ally}
data(MajorAlly)
Y <- MajorAlly
time <- dim(Y)[3]
drop.state <- c(which(colnames(Y) == "USA"), which(colnames(Y) == "CHN"))
newY <- Y[-drop.state, -drop.state, 1:62]
```

First, we fit a pilot model to elicit reasonable inverse gamma prior values for $\v_t$ ($v_0$ and $v_1$).
```{r test}
G <- 100
set.seed(1990)
test.run <- NetworkStatic(newY, R=2, mcmc=G, burnin=G, verbose=0,
                          v0=10, v1=time*2)
V <- attr(test.run, "V")
sigma.mu = abs(mean(apply(V, 2, mean)))
sigma.var = 10*mean(apply(V, 2, var))
v0 <- 4 + 2 * (sigma.mu^2/sigma.var)
v1 <- 2 * sigma.mu * (v0/2 - 1)
```

Then, we diagnose the break number by comparing model-fits of several models with a varying number of breaks.
```{r, fig.asp = 0.25, out.width="100%"}
set.seed(11223);
detect2 <- BreakDiagnostic(newY, R=2, break.upper=2,
                           mcmc=G, burnin=G, verbose=0,
                           v0=v0, v1=v1)
detect2[[1]]
```

The test results from WAIC, log marginal likelihood, and average loss indicate that HNC with two breaks is most reasonable.

Based on the test result, we fit the HNC with two breaks to the major power alliance network and save the result in *R* object `fit`.
```{r hncally}
G <- 100
K <- dim(newY)
m <- 2
initial.s <- sort(rep(1:(m+1), length=K[[3]]))
set.seed(11223);
fit <- NetworkChange(newY, R=2, m=m, mcmc=G, initial.s = initial.s,
                     burnin=G, verbose=0, v0=v0, v1=v1)
```

First, we can examine transitions of hidden regimes by looking at posterior state probabilities ($p(\S | \mathcal{Y}, \Theta)$) over time. `plotState()` in `MCMCpack` pacakge provides a function to draw the posterior state probabilities from changepoint analysis results. Since our input data is an array, we need to change the input data as a vector.

```{r, fig.asp = 0.8, out.width="100%"}
attr(fit, "y") <- 1:K[[3]]
plotState(fit, start=1)
```

Next, we draw regime-specific latent node positions of major powers using `drawPostAnalysis`. Users can choose the number of clusters in each regime by `n.cluster}.
```{r, fig.asp = 0.33, out.width="100%"}
p.list <- drawPostAnalysis(fit, newY, n.cluster=c(4, 4, 3))
multiplot(plotlist = p.list, cols=3)
```

Then, using `drawRegimeRaw()`, we can visualize original network connections for each regime by collapsing network data within each regime.

```{r, fig.asp = 0.33, out.width="100%"}
drawRegimeRaw(fit, newY)
```

Identifying hidden regimes of the military alliance network makes it clear the central role of Austria-Hungary during the first two regimes in the military alliance network among major powers.

# Acknowledgements

This work was supported by the Japan Society for the Promotion of Science Early-Career Scientists Grant [JP19K13606] to Y.S.

![Summary of selected features and functions of the package.\label{fig:list}](list.png)

# References