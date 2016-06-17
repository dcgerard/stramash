[![Build Status](https://travis-ci.org/stephens999/ashr.svg)](https://travis-ci.org/stephens999/ashr)

This repository contains an R package for performing "Adaptive Shrinkage."

To install the ashr package first you need to install devtools
```
install.packages("devtools")
library(devtools)
install_github("stephens999/ashr")
```

## Running Adaptive Shrinkage


The main function in the ashr package is `ash`. To get minimal help:
```
> library(ashr)
> ?ash
```

## More background

The ashr (``Adaptive SHrinkage") package aims to provide simple, generic, and flexible methods to derive ``shrinkage-based" estimates and credible intervals for unknown quantities $\beta=(\beta_1,\dots,\beta_J)$, given only estimates of those quantities ($\hat\beta=(\hat\beta_1,\dots, \hat\beta_J)$) and their corresponding estimated standard errors ($s=(s_1,\dots,s_J)$). Although shrinkage-based estimation can be motivated in various ways, our key goal here is to combine information across the multiple measurements $j=1,\dots,J$ to improve inference for each individual $\beta_j$. By improved inference, we mean both
improved average accuracy of point estimates, which is
the traditional focus of shrinkage-based methods, \emph{and} improved assessments of uncertainty. 

The ``adaptive" nature of the shrinkage is two-fold. First, the appropriate amount of shrinkage is determined from the data, rather than being pre-specified. Second, the amount of shrinkage undergone by each $\hat\beta_j$ will depend on the standard error $s_j$: measurements with high standard error will undergo more shrinkage than measurements with low standard error.

As an important special case, these methods address the "multiple comparisons" setting, where interest focuses on which $\beta_j$ can be confidently inferred to be non-zero. Such problems are usually tackled by computing a $p$ value for each $j$, often by applying a $t$ test to $\hat\beta_j/s_j$,
and then applying a generic procedure, such as that of Benjamini 
and Hochberg (1995?) or Storey (2001?), designed to control or
estimate the false discovery rate (FDR) or the positive FDR (Storey, 2001?). 
The ashr package provides analagous
generic methods that work directly with two numbers for each 
measurement $(\hat\beta_j,s_j$), rather than a single number (e.g.~ the $p$ value, or $t$ statistic). Working with these two numbers has two important benefits: first, it permits estimation and not only testing; second, the 
uncertainty in each measurement $\hat\beta_j$ can be more fully accounted for, reducing the impact of ``high-noise" measurements (large $s_j$) that can reduce the effectiveness of a standard FDR analysis. 
For more on the potential for shrinkage-based estimation to
address multiple comparisons see Greenland and Robins (1991),
Efron (2008) and Gelman et al (2012). [Note, check also Louis, JASA, 1984] 

### Methods Outline

The methods are based on treating the vectors $\hat\beta$
and $s$ as ``observed data", and then performing inference for $\beta$ from these observed data, using a standard hierarchical modelling framework
to combine information across $j=1,\dots,J$.

Specifically, 
we assume that the true 
$\beta_j$ values are independent
and identically distributed from some distribution $g(\cdot;\pi)$
where $\pi$ are hyperparameters to be estimated.

The key assumption is that $g$ is *unimodal* about zero.
We implement this by assuming that $g$ is a mixture
of uniforms, or a mixture of normals, each centered at 0.

Then, given $\beta$, we assume that $(\hat\beta_j,s_j)$ are independent across $j$, and depend only on $\beta_j$. Putting these together, the joint model for the unobserved $\beta$ and the observed $\hat\beta, s$ is:
\begin{align}
p(\hat\beta, s, \beta | \pi) & = \prod_j g(\beta_j ; \pi) p(\hat\beta_j, s_j | \beta_j) \\
& = \prod_j g(\beta_j ; \pi) L(\beta_j; \hat\beta_j,s_j).
\end{align}
The specific choices of $g$ and $L$ are described below.

We fit this hierarchical model using the following "Empirical Bayes" approach. First we estimate the hyper-parameters $\pi$ by maximizing the likelihood
$$L(\pi; \hat\beta, s) := p(\hat\beta, s | \pi) = \int p(\hat\beta, s, \beta | \pi) d\beta.$$
Then, given this estimate $\hat\pi$, we compute the conditional distributions $$p(\beta_j | \hat\pi, \hat\beta, s) \propto g(\beta_j; \pi) L(\beta_j; \hat\beta_j, s_j).$$ 
In principle we would
prefer to take a full Bayes approach that accounts for uncertainty in $\pi$, but, at least for now, we compromise this principle for the simplicity of the EB approach.
[Note: a Variational Bayes version of this is
also implemented, and may become our preferred approach
after testing]


These conditional distributions can be conveniently summarized
in various ways, including point estimates (e.g. the posterior means or medians),
and credible intervals/regions.


The key components of this hierarchical model
are the distribution $g$ and the likelihood $L(\beta_j; \hat\beta_j, s_j)$. We make the following choices for these.

1. The likelihood for $\beta_j$ is t, with known degrees of freedom, centered on $\hat\beta_j$, with scale parameter $s_j$.

2. The distribution $g(\cdot; \pi)$ is a mixture of zero-centered normal distributions, or zero-centered uniform distributions. For example, for the normals we would have:
$$g(\cdot; \pi) = \sum_{k=1}^K pi_k N(\cdot; 0, \sigma^2_k).$$
In practice, we currently fix the number of components $K$ to be large, and take the variances $\sigma_1<\sigma_2<\dots<\sigma_K$ to be fixed, and vary from very small (possibly 0), to very large --  sufficiently large that typically $\hat\pi_K=0$.


