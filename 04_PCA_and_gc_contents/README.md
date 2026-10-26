# PCA and sequence/signature features

The current portable implementation is [analyze_robustness_results.py](../07_robustness_analysis/analyze_robustness_results.py). Follow the [manuscript-analysis commands](../07_robustness_analysis/README.md) to provide your audited human, bat and experimental exports and generate PCA, GC3 and signature-feature tables/figures. The calculation uses the pinned Python environment and does not require access to the authors' cluster.

`reduce_signatures.py` stays at its original path because the frozen [signature-group manifest](../02_run_sinabro/manuscript_signature_groups.json) records its checksum. The current exporter already emits the 60-signature and 41-signature tables when given that manifest; rerunning the old reduction script is unnecessary.

The authors computed these panels together with downstream statistics on SeaWulf. The actual submitted script and execution record are retained in [07_robustness_analysis/author_hpc/](../07_robustness_analysis/author_hpc/README.md); this stage did not have a separate Slurm job. Earlier PCA/R plotting code is in [historical/04_PCA_and_gc_contents/](../historical/04_PCA_and_gc_contents/).
