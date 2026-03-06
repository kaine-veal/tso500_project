## Cross-system validation

A comparative analysis was performed to ensure that the tool behaves consistently across different systems, not only on the developer’s environment.

For this validation, the pipeline was installed and executed on a separate server using the same dummy VCF file. The resulting outputs were then compared using the comparison script provided in this repository.

Overall, **99.2% of the compared values were identical**. The small number of differences observed is mainly attributable to the use of different versions of external annotation databases (e.g. ClinVar).
