# The mathematics behind CppCrystal

This document sets out the mathematics the library rests on, in general terms.
Nothing here is about a particular space group; everything applies to all of
them at once, and to the layer, rod and point groups as well. The emphasis is on
the group theory, and in particular on the orbit-stabiliser theorem, which is
the single fact that ties multiplicities, site symmetries, Wyckoff positions,
irreducible k-points and structure generation together. Where an algorithm is
mentioned, it is described as a piece of mathematics, not as code.

Notation. $\mathbb{R}^3$ is Euclidean space, $E(3)$ its isometry group,
$O(3)$ the orthogonal group. Column vectors throughout. For a group $G$ acting
on a set $X$, $G\cdot x$ is the orbit of $x$ and $G_x$ its stabiliser. $[G:H]$
is the index of a subgroup, $|G|$ the order.

---

## 1. Affine maps and the Seitz calculus

Every isometry of $\mathbb{R}^3$ is an affine map

$$
x \;\mapsto\; R\,x + t ,\qquad R\in O(3),\ t\in\mathbb{R}^3 ,
$$

written in Seitz form as $(R\,|\,t)$. Composition and inversion follow from
substituting one map into the other:

$$
(R_1|t_1)(R_2|t_2) = (R_1R_2 \,|\, R_1t_2 + t_1),\qquad
(R|t)^{-1} = (R^{-1}\,|\,-R^{-1}t).
$$

The linear part $R$ is a group homomorphism $E(3)\to O(3)$; its kernel is the
group of pure translations $(I|t)$, which is therefore normal in $E(3)$.

### Fractional coordinates

A lattice basis is a matrix $L$ whose columns are three linearly independent
vectors $a,b,c$. Cartesian $r$ and fractional $x$ coordinates are related by
$r = Lx$. An isometry $(R|t)$ in Cartesian form becomes, in fractional form,

$$
(W\,|\,w) = (L^{-1}RL \;|\; L^{-1}t).
$$

Everything in the library is done in fractional coordinates, so an operation is
a pair $(W|w)$ with $W$ a $3\times3$ matrix and $w$ a vector. The composition
rule above is unchanged since it is purely algebraic.

### The metric tensor

The Euclidean inner product in fractional coordinates is
$\langle x, y\rangle = x^{\mathsf T} G\, y$ with $G = L^{\mathsf T}L$, the
metric tensor. A fractional matrix $W$ represents an isometry exactly when it
preserves $G$:

$$
W^{\mathsf T} G\, W = G .
$$

This is the condition used to search for symmetries without ever leaving the
lattice basis.

---

## 2. Lattices

A lattice is a discrete subgroup $T \subset \mathbb{R}^3$ of full rank,
$T = L\,\mathbb{Z}^3$. It is isomorphic to $\mathbb{Z}^3$, and the choice of
isomorphism is the choice of basis.

### Change of basis

Two bases $L$ and $L'$ span the same lattice exactly when $L' = LM$ with
$M \in GL(3,\mathbb{Z})$, i.e. $M$ integral with $\det M = \pm1$ (unimodular).
The set of all bases of one lattice is a single $GL(3,\mathbb{Z})$-orbit. A
*reduced* basis is a canonical representative of that orbit, chosen by
inequalities on the metric tensor:

* **Delaunay (Selling) reduction.** Extend the basis by $b_4 = -(b_1+b_2+b_3)$
  and require every pairwise inner product $b_i\cdot b_j \le 0$. The reduction
  is an iterative folding; the resulting *obtuse superbase* has the property
  that the shortest lattice vectors all lie among the seven vectors
  $b_1,b_2,b_3,b_4,b_1+b_2,b_2+b_3,b_3+b_1$. That finiteness is what makes the
  symmetry search of section 5 a bounded enumeration.
* **Niggli reduction.** The Krivý–Gruber conditions on the six parameters
  $(a\cdot a,\ b\cdot b,\ c\cdot c,\ 2\,b\cdot c,\ 2\,a\cdot c,\ 2\,a\cdot b)$.
  The Niggli cell is unique, which is what a canonical *standardized* cell
  needs.

Both reductions are compositions of unimodular steps, so the product of the
steps is the change of basis from the input to the reduced lattice.

### Sublattices, centering, index

If $T' \subseteq T$ is a sublattice with basis $L' = LM$, $M$ integral, then

$$
[T : T'] = |\det M| .
$$

A *conventional* cell is a sublattice $T_{\mathrm{conv}} \subseteq T$ chosen
for its symmetry-adapted shape; the *centering* vectors are the coset
representatives of $T_{\mathrm{conv}}$ in $T$, and there are exactly
$[T:T_{\mathrm{conv}}]$ of them (1, 2, 3 or 4 for the crystallographic
centerings). The primitive cell is a basis of $T$ itself.

### Reduced periodicity

A layer is periodic in a rank-2 lattice, a rod in a rank-1 lattice, a finite
cluster in the trivial lattice. The library keeps the three-vector formalism
throughout and records which coordinate axes are periodic. Every construction
below that says "modulo $T$" or "fold into the cell" applies only to the
periodic axes; on an aperiodic axis coordinates are compared exactly.

---

## 3. The symmetry group of a crystal

A crystal is a set of decorated points $\{(r_i, \tau_i)\}$ (position, species)
invariant under a lattice $T$. Its symmetry group is the stabiliser of the
crystal in $E(3)$:

$$
\mathcal{G} \;=\; \{\, g \in E(3) \;:\; g \text{ maps the crystal onto itself, species preserved} \,\}.
$$

This is a group because it is a stabiliser. Its pure translations form the
normal subgroup $T = \mathcal{G} \cap \{(I|t)\}$, and the linear parts form the
*point group*

$$
P \;=\; \{ W : (W|w) \in \mathcal{G} \text{ for some } w \} \;\cong\; \mathcal{G}/T .
$$

So a space group is an extension

$$
1 \longrightarrow T \longrightarrow \mathcal{G} \longrightarrow P \longrightarrow 1 ,
$$

with $T \cong \mathbb{Z}^3$ and $P$ finite. The extension splits (there is an
origin at which every operation has $w \in T$) exactly for the *symmorphic*
groups; otherwise some operation carries a translation not in $T$ for any
choice of origin, i.e. a screw axis or glide plane.

### Why the point group is finite and crystallographic

Since $P$ preserves the lattice, in fractional coordinates each $W \in P$ is an
integer matrix with $\det W = \pm1$. An integer matrix has an integer trace.
For a proper rotation by angle $\theta$, $\operatorname{tr} W = 1 + 2\cos\theta$,
so $2\cos\theta \in \mathbb{Z}$ and $\theta \in \{0, 60°, 90°, 120°, 180°\}$.
This is the crystallographic restriction: rotation orders $1, 2, 3, 4, 6$ only.
Each element of $P$ is classified by the pair $(\det W, \operatorname{tr} W)$,
which takes exactly ten values, the ten crystallographic rotation types
$1, 2, 3, 4, 6, \bar1, m, \bar3, \bar4, \bar6$. A subgroup of $GL(3,\mathbb{Z})$
preserving a positive-definite metric is finite, and the finite ones are the
32 crystallographic point groups; the largest has order 48.

### The finite quotient that does the work

$\mathcal{G}$ is infinite. For counting we pass to a finite quotient. Given a
conventional sublattice $T_{\mathrm{conv}} \subseteq T$, define

$$
G \;=\; \mathcal{G}/T_{\mathrm{conv}}, \qquad
|G| \;=\; |P|\cdot [T : T_{\mathrm{conv}}] .
$$

$G$ is finite. Its elements are the operations $(W|w)$ with $w$ taken modulo
$T_{\mathrm{conv}}$, i.e. the point-group operations *times* the centering
translations. This is what a crystallographic table calls the "general
position" list of a space group, and what the library calls the conventional
operation set. $G$ acts on the torus $\mathbb{R}^3 / T_{\mathrm{conv}}$, which
is the conventional cell with opposite faces identified. When
$T_{\mathrm{conv}} = T$ (primitive setting) $G = P$ as an abstract group, but
the action still remembers the translations $w$.

All of the orbit counting below takes place in $G$ acting on the torus. For a
layer or rod the torus has fewer periodic directions, and for a point group it
is $\mathbb{R}^3$ itself, but the statements are identical.

---

## 4. Group actions, orbits, stabilisers

Let a finite group $G$ act on a set $X$. For $x \in X$:

$$
G\cdot x = \{ g x : g \in G \}, \qquad
G_x = \{ g \in G : g x = x \}.
$$

$G_x$ is a subgroup (it is closed under composition and inverses). Orbits
partition $X$.

### The orbit-stabiliser theorem

> For every $x \in X$, $\;|G\cdot x| = [G : G_x]$, hence
> $\;|G\cdot x|\,\cdot\,|G_x| = |G|.$

*Proof.* Consider the map $\phi: G \to G\cdot x$, $g \mapsto gx$. It is onto by
definition of the orbit. Two elements have the same image exactly when
$g x = h x \iff h^{-1}g \in G_x \iff g \in h\,G_x$. So the fibres of $\phi$ are
precisely the left cosets of $G_x$, all of size $|G_x|$, and there are
$[G:G_x]$ of them. Lagrange's theorem $|G| = [G:G_x]\,|G_x|$ finishes it.
$\square$

Two immediate corollaries:

1. **Coset representatives generate the orbit.** If $g_1, \dots, g_m$ is one
   representative per left coset of $G_x$, then
   $G\cdot x = \{ g_1 x, \dots, g_m x \}$ with no repetitions. So an orbit can be
   built either by applying every element and discarding duplicates, or by
   applying only the coset representatives; the theorem says both give
   $m = |G|/|G_x|$ points.
2. **Stabilisers along an orbit are conjugate.**
   $G_{gx} = g\,G_x\,g^{-1}$. Every point of an orbit has an isomorphic
   stabiliser, and the stabilisers of one orbit form a single conjugacy class
   of subgroups.

### Burnside's counting lemma

The number of orbits of $G$ on a finite $X$ is

$$
|X/G| = \frac{1}{|G|} \sum_{g \in G} |\operatorname{Fix}(g)| ,
$$

where $\operatorname{Fix}(g) = \{x : gx = x\}$. This follows by counting the
pairs $(g, x)$ with $gx = x$ in two ways and applying orbit-stabiliser. It gives
the number of irreducible points of a symmetric mesh (section 7) without
enumerating them.

---

## 5. Applying the theorem to crystal positions

Take $G = \mathcal{G}/T_{\mathrm{conv}}$ acting on the torus
$X = \mathbb{R}^3/T_{\mathrm{conv}}$ as in section 3. For a point $x$ of the
cell:

* the orbit $G\cdot x$ is the set of symmetry-equivalent positions inside the
  conventional cell; its size is the **multiplicity** of $x$;
* the stabiliser $G_x$ is the **site-symmetry group** of $x$.

### The site-symmetry group is a point group

The restriction of the quotient map $\mathcal{G} \to P$, $(W|w) \mapsto W$, to
$G_x$ is injective: if $(W|w)$ and $(W|w')$ both fix $x$ then
$(I\,|\,w - w')$ fixes $x$, so $w - w' \in T_{\mathrm{conv}}$ and the two
operations coincide in $G$. Hence $G_x$ is isomorphic to a subgroup of $P$: a
site-symmetry group is always one of the 32 point groups (or a subgroup),
which is why it is written with a point-group symbol.

### Multiplicity times site-symmetry order

Orbit-stabiliser now reads

$$
\boxed{\;\text{multiplicity}(x)\;\times\;|G_x| \;=\; |G| \;=\; |P|\cdot[T:T_{\mathrm{conv}}]\;}
$$

This is the identity that the tables satisfy line by line, that the library's
group objects are checked against for every group in every family, and that
the refinement stage uses to decide which Wyckoff position an atom sits on.

### Fixed sets are affine subspaces

For a single operation $g = (W|w)$ acting on the torus,

$$
\operatorname{Fix}(g) = \{ x : Wx + w = x + n,\ n \in T_{\mathrm{conv}} \}
= \bigcup_{n} \{ x : (W - I)\,x = n - w \}.
$$

For each lattice vector $n$ the solution set is either empty or an affine
subspace parallel to $\ker(W - I)$, the $+1$ eigenspace of $W$. Its dimension
is 3 for the identity, 2 for a mirror, 1 for a rotation, 0 for an inversion or
rotoinversion. If $n - w$ never lies in the image of $W - I$ the operation has
no fixed point: it is a screw or a glide. Distinct lattice shifts $n$ can give
distinct, non-parallel-translated components inside one cell (a two-fold axis
at $x = 0$ and another at $x = \tfrac12$); on the torus these are different
pieces of $\operatorname{Fix}(g)$ and need not be related by any element of
$G$.

For a subgroup $H \le G$ the pointwise-fixed set is an intersection:

$$
\operatorname{Fix}(H) = \bigcap_{h \in H} \operatorname{Fix}(h),
$$

again a finite union of affine pieces. Conversely the pointwise stabiliser of
any set $A$, $\{g : gx = x\ \forall x \in A\}$, is a subgroup.

### Wyckoff positions as the stabiliser stratification

Two points are in the same **Wyckoff position** when their stabilisers are
conjugate in $G$. By corollary 2 above, a whole orbit lies in one Wyckoff
position, so a Wyckoff position is a union of orbits, all of the same
multiplicity. A Wyckoff position is therefore a *type* of orbit, and one orbit
is one instance of it.

The set of all points fixed by at least $H$ is $\operatorname{Fix}(H)$. Consider
the collection of all such sets, for all subgroups $H$: the *arrangement of
stabiliser loci*. It contains the whole space (for $H$ trivial), every
$\operatorname{Fix}(g)$, and is closed under intersection, since
$\operatorname{Fix}(H_1) \cap \operatorname{Fix}(H_2) = \operatorname{Fix}(\langle H_1, H_2\rangle)$.
Order it by inclusion. A point's stabiliser is read off the *smallest* locus
containing it: a point on a mirror plane but not on any axis inside that plane
has stabiliser $\{1, m\}$; a point at the intersection of the plane with an
axis has the larger group generated by both. So

* a **locus** $A$ of the arrangement, minus its proper sub-loci, is the set of
  points whose stabiliser is exactly the pointwise stabiliser $H_A$ of $A$;
* $G$ permutes the loci, because
  $h\cdot\operatorname{Fix}(g) = \operatorname{Fix}(hgh^{-1})$: the image of
  $g$'s fixed set under $h$ is the fixed set of the conjugate;
* the **Wyckoff positions are the $G$-orbits of loci**, and a locus and its
  image have conjugate stabilisers, $H_{hA} = h H_A h^{-1}$.

The **degrees of freedom** of a Wyckoff position are the dimension of its locus,
$\dim \ker(W - I)$ intersected across its stabiliser: the number of coordinates
that can vary while the stabiliser stays the same. The **general position** is
the locus "whole space", with trivial stabiliser, multiplicity $|G|$, and full
dimension. It always exists, and is the only Wyckoff position with trivial
stabiliser.

A *generic* point of a locus is any point of it not on a proper sub-locus. Its
stabiliser is exactly $H_A$; a non-generic point has a strictly larger one. The
orbit of a generic point has exactly $|G|/|H_A|$ points, and the points of the
locus that are generic form an open dense subset, so a random point of the
locus is generic with probability one.

### Projection onto a locus

An affine locus $A = a + \operatorname{span}(B)$ (with $B$ a $3\times d$ basis
of its directions, $d$ the degrees of freedom) comes with the affine projection

$$
\pi_A(x) = \Pi\, x + (a - \Pi a),\qquad \Pi = B\,(B^{\mathsf T}B)^{-1}B^{\mathsf T},
$$

which is idempotent and fixes $A$ pointwise. The tabulated *coordinate
operator* of a Wyckoff position (the triple such as $(x, x, 0)$ read as an
affine map) is exactly such a projection: its rank is the number of degrees of
freedom and its image is the locus. Projecting a point that is only
approximately on the locus onto it, then expanding the orbit, gives an exact
orbit of the correct multiplicity.

### Symmetrising a point onto its locus

If $H$ is a finite group of affine maps and $x$ any point, then

$$
\bar{x} = \frac{1}{|H|} \sum_{h \in H} h\,x
$$

is fixed by $H$: for any $h' \in H$,
$h'\bar x = \frac{1}{|H|}\sum_h h'h\,x = \bar x$ because $h \mapsto h'h$ is a
bijection of $H$. So averaging a point over its (approximately determined)
stabiliser projects it onto the exact fixed locus. On the torus the images
$h\,x$ must first be brought to the representatives nearest $x$, otherwise the
mean is taken across a lattice vector. This is the refinement step that turns
a numerically noisy special position into an exact one.

---

## 6. Finding the symmetry of a given crystal

The determination problem is the inverse of the description above: given
$(L, \{x_i\}, \{\tau_i\})$, find $\mathcal{G}$. Everything is up to a tolerance
$\varepsilon$, so "equal" means "within $\varepsilon$ modulo the lattice".

### The lattice point group

The candidates for the linear parts are the elements of the stabiliser of the
metric tensor in $GL(3,\mathbb{Z})$,

$$
P_L = \{\, M \in GL(3,\mathbb{Z}) : M^{\mathsf T} G M = G \,\},
$$

the *holohedry* of the lattice. An element of $P_L$ sends each basis vector to
a lattice vector of the same length, and in a Delaunay-reduced basis the
vectors of any given short length lie in a small explicit set (the 26 vectors
with coordinates in $\{-1,0,1\}$ suffice). So $P_L$ is found by trying every
triple of candidate images, keeping those that are unimodular and preserve $G$
within tolerance. $|P_L| \le 48$; more than that means the tolerance is too
loose to separate the lattice from a higher-symmetry one. The result is
transported to the input basis by the similarity transform of section 1.

$P \subseteq P_L$: the crystal's point group is a subgroup of its lattice's.

### Translation parts

For each $W \in P_L$, the translations $w$ such that $(W|w)$ permutes the
decorated points are wanted. A translation is fixed by the image of a single
point: if $(W|w)$ maps atom $i_0$ onto atom $j$ then $w = x_j - W x_{i_0}$
modulo $T$. Choosing $i_0$ among the rarest species, there are at most as many
candidates as there are atoms of that species; each candidate is verified
against every atom. The set of $(W|w)$ that pass is the symmetry set of the
cell in $G = \mathcal{G}/T_{\mathrm{input}}$.

In exact arithmetic this set is automatically a group. Under tolerance it is
only approximately closed; the refinement stage (symmetrise positions, re-run
with the exact positions) restores exact closure.

### The true translation lattice and the primitive cell

The operations with $W = I$ found above are the pure translations
$T/T_{\mathrm{input}}$: the input cell may be a supercell. These finitely many
vectors, together with the input basis, generate $T$; a basis of $T$ is chosen
from them (any three that span the volume $|\det L| / [T:T_{\mathrm{input}}]$),
Delaunay-reduced, and the atoms folded and de-duplicated into it. The number of
atoms drops by the index, and every operation's translation is re-expressed.

### Identifying the point group

The multiset of rotation types $(\det W, \operatorname{tr} W)$ over $W \in P$
is a class function, invariant under any change of basis. It is a complete
invariant for the 32 crystallographic point groups: no two have the same
histogram of the ten types. So the point group is identified by a histogram
lookup, with no need to find an explicit isomorphism. Choosing conventional
axes, the *arithmetic* part of the identification, means choosing a basis in
which the distinguished rotation axis and the mirrors line up with the
tabulated setting; the axes are found as eigenvectors ($\ker(W - I)$ for a
proper rotation, a primitive integer vector) and ordered by the conventions of
the Laue class.

### Identifying the space-group type

Two space groups are of the same *type* when they are conjugate by an affine
map $(M|s)$ with $M \in GL(3,\mathbb{Z})$ (bases of the same lattice) and $s$
an origin shift; up to that equivalence there are 230 types in three dimensions
(219 if $M$ is restricted to $\det M = +1$). Given $G$ in some conventional
setting, identification means finding $(M|s)$ that carries $G$ onto one of the
tabulated settings.

The linear part $M$ ranges over a finite list of candidate axis and setting
choices for the crystal system, so it is enumerated. For each candidate the
origin shift is a linear problem. Under conjugation by a translation,

$$
(I|s)\,(W|w)\,(I|s)^{-1} = (W \,|\, w - (W - I)\,s),
$$

so if the found operation is $(W_i|w_i)$ and the tabulated one with the same
linear part is $(W_i|w_i^0)$, the shift must satisfy

$$
(W_i - I)\, s \;\equiv\; w_i - w_i^0 \pmod{T_{\mathrm{conv}}}
$$

for every generator $i$ of the tabulated setting. Stacking three generators
gives a $9\times3$ integer system $A s = d$. With a Smith-type decomposition
$A = V\,S\,U$ ($V, U$ unimodular, $S$ diagonal) the system is solvable modulo
the lattice exactly when the components of $V^{-1}d$ along zero rows of $S$ are
integral, and then $s = U^{-1} S^{+} V^{-1} d$. The decomposition depends only
on the tabulated generators and is precomputed once per setting. Because two
found operations with the same $W$ differ by a lattice vector, the
right-hand side is well defined modulo the lattice; every representative is
tried so the answer does not depend on the order the operations were found in.
A candidate $(M|s)$ is accepted when it maps *every* operation of $G$ onto a
distinct tabulated operation.

The stabiliser of the tabulated $G$ under this action, the set of $(M|s)$ that
map the group onto itself, is its *affine normaliser*; it is why the same group
can be matched at several origin shifts, and why the first successful choice is
taken by a fixed preference order rather than being unique.

### Assigning Wyckoff letters

Once the setting and origin are fixed, each independent atom is matched to a
tabulated Wyckoff position. A candidate position with multiplicity $m$ and
site-symmetry generator set $S$ is consistent with the atom's orbit exactly when
some point of the orbit is fixed by $S$ and the count of orbit points
coinciding with it, times $m$, equals $|G|$: this is orbit-stabiliser applied
in reverse, using the theorem to reject positions whose stabiliser is too
small or too large for the observed orbit.

---

## 7. Reciprocal space and irreducible points

The reciprocal lattice is the dual $T^* = \{ k : k\cdot t \in \mathbb{Z}\ \forall t \in T\}$,
with basis $L^{-\mathsf T}$. If real-space fractional coordinates transform by
$W$, reciprocal fractional coordinates transform by $W^{-\mathsf T}$, the
contragredient. Translations act trivially on $k$, so the acting group is the
point group $P$ (or its contragredient image, which is isomorphic). Time
reversal sends $k \mapsto -k$; including it means acting by
$P \cup (-I)P$, the Laue group of $P$ when $-I \notin P$, and $P$ itself
otherwise.

A mesh $\{ (i/N_1, j/N_2, l/N_3) + \text{shift} \}$ is a finite set. When the
group maps the mesh onto itself the mesh is a finite $P$-set, and

* the **irreducible k-points** are the orbit representatives,
* the **weight** of each is its orbit size, which by orbit-stabiliser divides
  $|P|$ and equals $|P|/|P_k|$; points on symmetry elements have smaller orbits,
* the number of irreducible points is given by Burnside's lemma,
* a canonical representative is the smallest grid index in the orbit, reachable
  by applying every element since the group is closed.

If the mesh is not preserved by the group (a rotation sends a grid point off the
grid) the action is computed in a finer integer lattice and only images that
land on grid points are accepted; that is the "distortion" path. Brillouin-zone
folding maps each grid point to its representative of minimal length among the
translates $k + K$, $K \in T^*$; ties on the zone boundary are kept as the set
of equally short translates.

---

## 8. Magnetic and spin groups

A magnetic crystal decorates each site with a moment $m_i$, an axial vector.
Under an isometry with linear part $W$, an axial vector transforms by
$\det(W)\,W$, and time reversal $1'$ negates it. The magnetic symmetry group is
the stabiliser of the decorated crystal in $E(3) \times \{1, 1'\}$.

Let $\mathcal{M}$ be that group, $F$ its image in $E(3)$ (time reversal
forgotten, the *family* group), and $D \le F$ the operations that occur
unprimed (the *maximal ordinary* subgroup). Then $[F : D] \in \{1, 2\}$, because
the primed operations, if any, form a single coset of $D$ in $F$: the product of
two primed operations is unprimed. The four Opechowski–Guccione types follow
from that index and from whether pure time reversal is present:

* $[F:D] = 1$ and $1' \notin \mathcal{M}$: type I (an ordinary space group);
* $[F:D] = 1$ and $1' \in \mathcal{M}$: type II (grey);
* $[F:D] = 2$ and the coset representative is a primed point operation: type III;
* $[F:D] = 2$ and the coset representative is a primed pure translation
  (an anti-translation): type IV.

Identification then proceeds as in section 6 on $F$ (types I–III) or on $D$
(type IV), and the coset structure is transported by the same conjugation.

For a *spin* space group the spatial and spin actions are decoupled: an
operation is a pair (spatial isometry, spin rotation). For collinear moments the
spin rotations reduce to a sign, and the search asks, for each spatial operation
that permutes the sites, whether a single sign $\pm1$ carries every moment onto
the moment of its image. Sites with zero moment do not constrain the sign, and
an operation constrained by none of them is admitted with both.

---

## 9. Subgroups

A subgroup $H \le \mathcal{G}$ is *translationengleiche* (a t-subgroup) when it
keeps the whole translation lattice, $H \cap T = T$, and *klassengleiche* (a
k-subgroup) when it keeps the point group and loses translations. Hermann's
theorem says a maximal subgroup is one or the other.

The t-subgroups of $\mathcal{G}$ correspond one-to-one with the subgroups of
$P$: for $Q \le P$, $H_Q = \{ (W|w) \in \mathcal{G} : W \in Q \}$, and
$[\mathcal{G} : H_Q] = [P : Q]$. So the maximal t-subgroups are read off the
maximal subgroups of the point group. The subgroup lattice of a finite $P$ is
generated by its cyclic subgroups and closed under join
($\langle Q_1 \cup Q_2 \rangle$); iterating joins to a fixed point enumerates
every subgroup, including those that are not generated by two elements. The
maximal ones are the proper subgroups not strictly contained in another proper
subgroup. Each surviving operation set $H_Q$ is a space group in its own right,
identified by the procedure of section 6, which yields an edge
$\mathcal{G} \to H_Q$ with its index. The transitive closure of the maximal
relation is the full t-subgroup relation.

---

## 10. Generating a crystal with prescribed symmetry

Generation is the inverse problem of section 5: given $G$ and a composition
(species $\tau$ with $n_\tau$ atoms each), produce a crystal whose symmetry
group is $G$.

### Seating a composition on Wyckoff positions

Each atom must lie on some orbit, and an orbit is either a whole Wyckoff
position with zero degrees of freedom (a fixed, finite set of points, usable at
most once in the entire structure since a second use would coincide with the
first) or one instance of a position with free coordinates (usable any number
of times with different parameters). So for each species one must choose a
multiset of Wyckoff positions $\{A_j\}$ with

$$
\sum_j \text{multiplicity}(A_j) = n_\tau ,
$$

subject to the once-only rule across all species. This is a bounded partition
problem; the search is a depth-first walk that prunes a branch as soon as the
remaining positions cannot supply the remaining count. A composition is
*compatible* with $G$ exactly when at least one assignment exists; for
instance a count smaller than the minimal multiplicity is never compatible.

### Sampling and expansion

For each chosen instance of a position with $d$ degrees of freedom, $d$
parameters are drawn, the point $a + B\,\theta$ is formed, and its orbit is
expanded through the coset representatives of its stabiliser, folding only the
periodic axes. A generic point (probability one under a continuous
distribution) has exactly the stabiliser of its locus, so the orbit has the
tabulated multiplicity and the structure's symmetry group is $G$ exactly, no
larger.

Under a finite tolerance, though, the set of parameters that a determination
would *read* as having extra symmetry has positive measure, and it grows with
the number of atoms on special sites. Using only the general position avoids
this: a structure whose atoms all lie on trivial-stabiliser orbits has no
mechanism to acquire a mirror or an inversion by accident.

### Distances on the torus

Two atoms at $x, y$ are at distance
$\min_{n \in T} \| L\,(x - y - n) \|$, the *minimal image*; for a reduced basis
the minimum is attained with $n$ in $\{-1,0,1\}^3$ on the periodic axes and
$n = 0$ on the aperiodic ones. A generated structure is accepted only if every
minimal-image distance exceeds a species-dependent threshold; lattice and
parameters are resampled otherwise.

---

## 11. What is exact and what is approximate

The group theory of sections 3–5 is exact and applies to the tabulated groups
as given: integer matrices, rational translations, exact fixed loci computed
by integer linear algebra. The only tolerance there is the one that compares
a projected generic seed to its images, and for rational data with an
irrational-looking seed it never matters.

The determination of sections 6–7 is approximate in two places: comparing
metric tensors, and comparing positions modulo the lattice. Consequences:

* the found operation set is only approximately a group; closure is restored
  by refinement, and a *count* of found operations that exceeds the maximum
  (48 rotations, or a multiplicity product that disagrees with $|G|$) is
  a signal that the tolerance is too loose;
* the same crystal may be reported with different symmetry at different
  tolerances, and both answers are correct for their $\varepsilon$;
* an origin shift or standardisation is unique only up to the affine
  normaliser, so a canonical choice among equivalent answers is a convention,
  not a theorem.

---

## 12. Invariants that validate the implementation

Several of the identities above are cheap to check and independent of any
external reference; the test suite asserts them across every group of every
family.

* **Orbit-stabiliser, per Wyckoff position.** For every group $G$ in the
  catalogue and every position $A$ of it,
  $\text{multiplicity}(A) \cdot |H_A| = |G|$. The two factors come from
  different places: for space and layer groups the multiplicity is tabulated
  while the stabiliser is computed by applying the tabulated operations to a
  generic point of the locus; for point and rod groups both are derived from
  the operations alone. Section 4 says the identity holds for *any* action of
  *any* finite group at *any* point, so what the check certifies is that the
  operations form a group modulo the lattice, that the folding and equality
  used to compare images are consistent with that group structure, and (where
  the multiplicity is tabulated) that the locus, the operations and the table
  agree with one another.
* **General position.** The last position has multiplicity $|G|$, full
  dimension, and a trivial stabiliser. A non-trivial stabiliser at a generic
  point would mean a repeated operation in the table.
* **Origin.** For a point group the origin is fixed by everything, so the first
  position has multiplicity one.
* **Orders.** The point-group orders are the textbook values, and the number of
  operations of a space group in its conventional setting is
  $|P|\cdot[T:T_{\mathrm{conv}}]$.
* **Self-consistency of generation.** A generated structure is invariant under
  every operation of its target group, and the determination pipeline recovers
  the target group from it. The second of these exercises sections 5, 6 and 10
  against one another without any external oracle.
* **Oracle comparison.** Independently of the invariants, the determination,
  standardisation, reciprocal-mesh and magnetic outputs are compared with a
  reference implementation over a corpus of structures.

The invariants are necessary conditions. What they do not establish is
completeness (that every point of the cell lies in *some* listed position) or
the identification of a derived position with a particular tabulated letter;
those rest on the tabulated data for space and layer groups and on the
arrangement construction of section 5 for point and rod groups.
