# Preprocessing of reads by mapping them to the reference genome and collapsing them into UMI count matrices (10x's cellranger)

This folder has been created for internal purposes only, mainly for the review of this section to be performed by Kevin.
Major things to check:
- Correctness of the CITE-seq sample sheets. I've left two examples, for sample-id-1 and sample-id-2 within sample project <code>sample_project_1</code>
- Run of at least one gene expression sample per sequencing run.

### sequencing_data/

* <code>date_1/</code>
  * <code>aggr/</code>
  * <code>count/</code>
    * <code>data/</code>
        * <code>experiment_layout/</code>
            * <code>sample_project_1/</code>
                * <code>**CITE-seq_sample_sheet_for_feat_barcoding_NV0XX_sample-id-1_CITE.csv**</code>
                * <code>**CITE-seq_sample_sheet_for_feat_barcoding_NV0XX_sample-id-2_CITE.csv**</code>
            * sample_project_2/
                * ...
    * sample_project_1/
    * sample_project_2/
  * mkfastq/
  * jobs_scripts/
    * aggr/
    * count/
        * sample_project_1/
            * sample-id-1_Gex
                * **sample-id-1_Gex_cellranger_count.job.sh**
            * sample-id-1_CITE
            * sample-id-2_Gex
            * sample-id-2_CITE
        * sample_project_2/
            * ...
    * mkfastq/
* date_2/
  * aggr/
  * count/
  * mkfastq/
* date_n/
    * ...