# AvrutinSearch
A Julia implementation via interval arithmetic of the algorithm described in *Avrutin V et al. Calculation of homoclinic and heteroclinic orbits in 1D maps. Commun Nonlinear Sci Numer Simulat (2014), http://dx.doi.org/10.1016/j.cnsns.2014.07.008*.

# Usage

## Finding homoclinic orbits
The algorithm from the paper is given below:

![Listing 1 from *Avrutin V et al.*](AvrutinSearch/docs/resources/Listing1_pseudocode.png)

You must supply $f$, $x^\*$, $\mathcal{T}\_0$, $\mathcal{I}$, $r\_{\max}$, and another parameter $n\_{\text{iter}}$ to the function `AvrutinSearch.homoclinic_to_equilibrium`.

- $f$ should be a function `f(x::Float64)::Float64`.
  We might later supply helper functions for converting discrete-sampled maps to such functions by various interpolation methods, but you can for now use analytic descriptions for maps or implement your own linear interpolation.
  Some restrictions on $f$ apply but are discussed elsewhere (e.g., in the description of the invariant interval $\mathcal{I}$).

- $x^\*$ should be a `Float64` value. Make sure that this is a repelling fixed point of the map!
  The easiest way to check this is by representing $x^\*$ as a variable `x` and running `f(x) == x`; this should evaluate to `true`.

- $\mathcal{T}\_0$ is the initial target subset of $\mathbb{W}^u(x^\*)$ being searched for.
  Practically, $\mathcal{T}\_0$ is supplied as an `Interval` value and should contain $x^\*$.
  Some guidelines for selecting $\mathcal{T}\_0$ are given in the paper by Avrutin et al.:

  > ![Determination of target sets from *Avrutin V et al.*](AvrutinSearch/docs/resources/target_set_determination.png)

  As $\{f\_j^{-1}, \mathcal{V}\_j\}\_{j=1}^k$ cannot currently be supplied, we only implement the calculation of forward iterates of $\mathcal{T}\_0$; please refer to $n\_{\text{iter}}$ for more on this.

- $\mathcal{I}$ is some closed bounded interval $[a, b]$ satisfying $x^\* \in \mathcal{I}$ and $f(\mathcal{I}) = \mathcal{I}$.
  In practice, this should be a `Tuple{Float64, Float64}` value `(a, b)`.
  Currently, you do need to determine this interval by manual inspection of $f$; in the future we may implement an algorithm for determination of this interval under certain assumptions about $f$ (e.g., $\mathcal{C^1}$ differentiability).
  
  Note that $\mathcal{I}$ can and should be chosen to be as small as possible; as long as $\mathcal{I}$ satisfies the two aforementioned requirements, any points lying outside of $\mathcal{I}$ cannot be involved in homoclinic orbits.
  Moreover, there should be some positive integer $k$ such that for each $x \in \mathcal{I}$ the pre-image set $f|\_\mathcal{I}^{-1}(x)$ has at most $k$ many elements; otherwise the pre-image search tree may become intractibly broad.
  In particular, $f$ should have $k$ many interval ``branches'' $\mathcal{V}\_1, \dots, \mathcal{V}\_k$ in $\mathcal{I}$ on which its restriction is invertible.

- $\mathcal{r\_{\max}}$ should be a `UInt64` value.
  This is a maximal depth in the pre-image search tree.
  If $\mathcal{r\_{\max}}$ is zero, there will be no maximal rank of pre-image imposed (which is probably not the best idea).

- $n\_{\text{iter}}$ should be a `UInt64` value.
  This controls how many forward iterates of the initially supplied target interval $\mathcal{T}\_0$ are taken as a true target set $\mathcal{T}$.
  If zero, then no forward iterates are used.
  More iterates mean more memory consumption for storing the multiple interval (connected) components of $\mathcal{T}$, but less pre-image computation steps during search; this can drastically improve search time.

$\{f\_j^{-1}, \mathcal{V}\_j\}\_{j=1}^k$ is computed automatically, but in the future, we may allow these data to be manually supplied.