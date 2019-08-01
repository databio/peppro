# Browse the output directory

This is an interactive display of exactly what results you'll get as output from PEPPRO:

* ??? danger ":fa-folder-open-o: reports/"
    * :fa-file-code-o: [fastp_report.html](../files/examples/tutorial/reports/fastp_report.html)
    * :fa-file-code-o: [fragment_distribution.html](../files/examples/tutorial/reports/fragment_distribution.html)
    * :fa-file-code-o: [library_complexity.html](../files/examples/tutorial/reports/library_complexity.html)
    * :fa-file-code-o: [minus_frif.html](../files/examples/tutorial/reports/minus_frif.html)
    * :fa-file-code-o: [plus_frif.html](../files/examples/tutorial/reports/plus_frif.html)
    * :fa-file-code-o: [mrna_contamination.html](../files/examples/tutorial/reports/mrna_contamination.html)
    * :fa-file-code-o: [objects.html](../files/examples/tutorial/reports/objects.html)
    * :fa-file-code-o: [pause_index.html](../files/examples/tutorial/reports/pause_index.html)
    * :fa-file-code-o: [samples.html](../files/examples/tutorial/reports/samples.html)
    * :fa-file-code-o: [status.html](../files/examples/tutorial/reports/status.html)
    * :fa-file-code-o: [tss_enrichment.html](../files/examples/tutorial/reports/tss_enrichment.html)
    * :fa-file-code-o: [tutorial.html](../files/examples/tutorial/reports/tutorial.html)
* ??? danger ":fa-folder-open-o: results_pipeline/"
    * ??? danger ":fa-folder-open-o: tutorial/"
        * ??? danger ":fa-folder-open-o: aligned_hg38/"
            * :fa-file-o: tutorial_fail_qc.bam
            * :fa-file-o: tutorial_plus.bam
            * :fa-file-o: tutorial_minus.bam
            * :fa-file-o: tutorial_sort.bam
            * :fa-file-o: tutorial_unmap.bam
        * ???+ danger ":fa-folder-open-o: fastq/"
        * ??? danger ":fa-folder-open-o: fastqc/"
            * [:fa-file-code-o: tutorial_rmAdapter.html](../files/examples/tutorial/results_pipeline/tutorial/fastqc/tutorial_rmAdapter.html)
            * [:fa-file-code-o: tutorial_rmAdapter.json](../files/examples/tutorial/results_pipeline/tutorial/fastqc/tutorial_rmAdapter.json)
            * [:fa-file-text-o: tutorial_rmAdapter.txt](../files/examples/tutorial/results_pipeline/tutorial/fastqc/tutorial_rmAdapter.txt)
        * ??? danger ":fa-folder-open-o: prealignments/"
            * :fa-file-text-o: tutorial_rCRSd_3k_bt_aln_summary.log
            * :fa-file-archive-o: tutorial_rCRSd_3k_unmap_R1.fq.gz
            * :fa-file-archive-o: tutorial_rCRSd_3k_unmap_R2.fq.gz
            * :fa-file-text-o: tutorial_human_rDNA_bt_aln_summary.log
            * :fa-file-archive-o: tutorial_human_rDNA_unmap_R1.fq.gz
            * :fa-file-archive-o: tutorial_human_rDNA_unmap_R2.fq.gz
         * ??? danger ":fa-folder-open-o: QC_hg38/"
            * :fa-file-o: tutorial_bamQC.tsv
            * :fa-file-text-o: tutorial_fragCount.txt
            * :fa-file-text-o: tutorial_fragLen.txt
            * [:fa-file-pdf-o: tutorial_fragLenDistribution.pdf](../files/examples/tutorial/results_pipeline/tutorial/QC_hg38/tutorial_fragLenDistribution.pdf)
            * [:fa-file-image-o: tutorial_fragLenDistribution.png](../files/examples/tutorial/results_pipeline/tutorial/QC_hg38/tutorial_fragLenDistribution.png)
            * :fa-file-text-o: tutorial_fragLenDistribution.txt
            * [:fa-file-pdf-o: tutorial_minus_frif.pdf](../files/examples/tutorial/results_pipeline/tutorial/QC_hg38/tutorial_minus_frif.pdf)
            * [:fa-file-pdf-o: tutorial_minus_frif.png](../files/examples/tutorial/results_pipeline/tutorial/QC_hg38/tutorial_minus_frif.png)
            * :fa-file-text-o: tutorial_minus_TssEnrichment.txt
            * [:fa-file-pdf-o: tutorial_mRNA_contamination.pdf](../files/examples/tutorial/results_pipeline/tutorial/QC_hg38/tutorial_mRNA_contamination.pdf)
            * [:fa-file-pdf-o: tutorial_mRNA_contamination.png](../files/examples/tutorial/results_pipeline/tutorial/QC_hg38/tutorial_mRNA_contamination.png)
            * [:fa-file-pdf-o: tutorial_plus_frif.pdf](../files/examples/tutorial/results_pipeline/tutorial/QC_hg38/tutorial_plus_frif.pdf)
            * [:fa-file-pdf-o: tutorial_plus_frif.png](../files/examples/tutorial/results_pipeline/tutorial/QC_hg38/tutorial_plus_frif.png)
            * :fa-file-text-o: tutorial_plus_TssEnrichment.txt           
            * :fa-file-text-o: tutorial_preseq_counts.txt
            * :fa-file-text-o: tutorial_preseq_coverage.txt
            * :fa-file-text-o: tutorial_preseq_out.txt
            * [:fa-file-pdf-o: tutorial_preseq_plot.pdf](../files/examples/tutorial/results_pipeline/tutorial/QC_hg38/tutorial_preseq_plot.pdf)
            * [:fa-file-pdf-o: tutorial_preseq_plot.png](../files/examples/tutorial/results_pipeline/tutorial/QC_hg38/tutorial_preseq_plot.png)
            * :fa-file-text-o: tutorial_preseq_yield.txt
            * [:fa-file-pdf-o: tutorial_TSSenrichment.pdf](../files/examples/tutorial/results_pipeline/tutorial/QC_hg38/tutorial_TSSenrichment.pdf)
            * [:fa-file-pdf-o: tutorial_TSSenrichment.png](../files/examples/tutorial/results_pipeline/tutorial/QC_hg38/tutorial_TSSenrichment.png)
        * ??? danger ":fa-folder-open-o: raw/"
            * :fa-file-archive-o: tutorial_r1.fastq.gz
            * :fa-file-archive-o: tutorial_r2.fastq.gz
            * :fa-file-archive-o: hg38_annotations.bed.gz 
        * ??? danger ":fa-folder-open-o: signal_hg38/"
            * :fa-file-o: tutorial_minus_body_0-mer.bw
            * :fa-file-o: tutorial_plus_body_0.bw      
        * [:fa-file-text-o: objects.tsv](../files/examples/tutorial/results_pipeline/tutorial/objects.tsv) 
        * :fa-file-code-o: PEPPRO_cleanup.sh
        * [:fa-file-code-o: PEPPRO_commands.sh](../files/examples/tutorial/results_pipeline/tutorial/PEPPRO_commands.sh)
        * :fa-file-o: PEPPRO_completed.flag
        * [:fa-file-o: PEPPRO_log.md](../files/examples/tutorial/results_pipeline/tutorial/PEPPRO_log.md)
        * [:fa-file-text-o: PEPPRO_profile.tsv](../files/examples/tutorial/results_pipeline/tutorial/PEPPRO_profile.tsv)
        * [:fa-file-text-o: stats.tsv](../files/examples/tutorial/results_pipeline/tutorial/stats.tsv)
* ??? danger ":fa-folder-open-o: submission/"
    * :fa-file-code-o: peppro.py_tutorial.sub
    * :fa-file-code-o: tutorial.yaml
    * :fa-file-text-o: peppro.py_tutorial.log
* ??? danger ":fa-folder-open-o: summary/"
    * :fa-file-pdf-o: [tutorial_libComplexity.pdf](../files/examples/tutorial/summary/tutorial_libComplexity.pdf)
    * :fa-file-image-o: [tutorial_liComplexity.png](../files/examples/tutorial/summary/tutorial_libComplexity.png)
* :fa-file-text-o: [tutorial_stats_summary.tsv](../files/examples/tutorial/tutorial_stats_summary.tsv)
* :fa-file-code-o: [tutorial_summary.html](../files/examples/tutorial/tutorial_summary.html)
