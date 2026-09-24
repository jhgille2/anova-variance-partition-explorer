# ANOVA Variance Partition Explorer

An interactive, dependency-free web application that visualizes how analysis
of variance partitions phenotypic variation among sources in a balanced
multi-environment agricultural trial — and how that partitioning determines
heritability.

**Live app:** https://jhgille2.github.io/anova-variance-partition-explorer/

## What it shows

1. **ANOVA table** — df, SS, MS, F, p, and method-of-moments variance-component
   estimates for genotype, environment, genotype × environment,
   rep(environment), and residual.
2. **Variance partitioning** — grouped bars comparing the hidden true variance
   components (used to simulate the trait) against the ANOVA estimates from the
   trial, showing how well ANOVA recovers each source.
3. **Heritability** — true vs. estimated entry-mean h², plus theoretical curves
   of h² against replications and against environments, marking the current
   design.
4. **Estimator behavior** — histogram of the ĥ² estimate over Monte Carlo
   repetitions of the whole trial, showing how design choices affect the
   precision (and truncation frequency) of variance partitioning.

## Controls

- **Trial design:** number of lines (genotypes), replications, environments.
- **Hidden trait parameters:** mean μ and the five variance components
  σ²_G, σ²_E, σ²_GE, σ²_R, σ²_ε. A "target h²" solver back-computes the σ²_G
  needed for a desired entry-mean heritability under the current design.
- **Simulation:** Monte Carlo trial count and random seed (fully reproducible).

## Statistical model

Balanced g × e × r trial, RCBD within environments (reps nested in
environments), all effects random and independent:

y_ijk = μ + G_i + E_j + GE_ij + R(E)_jk + ε_ijk

- Variance components by the ANOVA method of moments from the expected mean
  squares, e.g. σ̂²_G = (MS_G − MS_GE)/(r·e).
- Entry-mean heritability: h² = σ²_G / (σ²_G + σ²_GE/e + σ²_ε/(r·e))
  (Holland, Nyquist & Cervantes-Martinez, 2003).
- With r = 1, GE and residual are confounded; MS_GE serves as the residual,
  as is conventional for unreplicated trials.

See the "Statistical methodology" section in the app for the full EMS table,
estimators, assumptions, and references.

## Run locally

Any static file server works, e.g.:

```sh
python3 -m http.server 8000
# open http://localhost:8000/
```

No build step, no dependencies.

## Files

- `index.html` — page structure and controls
- `css/styles.css` — styling
- `js/anova-core.js` — pure statistical core (simulation, ANOVA, variance
  components, heritability, F-distribution, Monte Carlo); also loadable in
  Node for testing
- `js/app.js` — UI wiring and canvas rendering

## References

- Holland JB, Nyquist WE, Cervantes-Martinez CT (2003) Estimating and
  interpreting heritability for plant breeding: an update.
  *Plant Breeding Reviews* 22:9–112. https://doi.org/10.1002/9780470650202.ch2
- Fehr WR (1987) *Principles of Cultivar Development, Vol. 1: Theory and
  Technique.* Macmillan.
- Hallauer AR, Miranda JB (1988) *Quantitative Genetics in Maize Breeding*,
  2nd ed. Iowa State University Press.

## License

MIT.
