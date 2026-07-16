# GADMA2 Parameter Files

Three GADMA2 optimization sets were run to infer the joint demographic history
of *S. incompletum*, *S. sandwicense* (Oahu), and *S. kavaiensis* (Kauai).
All runs used the moments engine with a three-population [1,1,1] model structure
and asymmetric migration as free parameters. Full methodological details are
provided in the manuscript Methods and Supplementary Methods.

## Run Order

**Run A** (`runA_Si_outgroup_primary.txt`) — Primary optimization under the
*S. incompletum*-outgroup topology (Si(Ss,Sk)). 10 replicates from random
starting points. The lowest-AIC replicate from this run is the primary best-fit
model reported in the manuscript.

**Run B** (`runB_Sk_outgroup_topology_test.txt`) — Topology test using an
alternative *S. kavaiensis*-outgroup topology (Sk(Si,Ss)), with a separately
generated SFS reflecting the alternative population label order. 10 replicates.
Compare best AIC from Run B against Run A to confirm topology; the Si-outgroup
topology was decisively supported (ΔAIC = 269.51).

**Run C** (`runC_Si_outgroup_convergence_mixed.txt`) — Convergence confirmation
under the Si-outgroup topology. 5 replicates warm-started from the Run A
best-fit model and 5 replicates from random starting points. Agreement between
warm and cold starts on key demographic parameters (Ne, divergence times) was
taken as evidence of convergence on the global optimum.

## Execution

```bash
gadma --params runA_Si_outgroup_primary.txt
gadma --params runB_Sk_outgroup_topology_test.txt
gadma --params runC_Si_outgroup_convergence_mixed.txt
```

Run A and B can be run independently. Run C should be run after Run A, as it
requires the best-fit model file from Run A for warm-starting (specified via
`Custom filename` in the params file). Update the `Custom filename` and
`Output directory` paths in each params file to match your local directory
structure before running.
