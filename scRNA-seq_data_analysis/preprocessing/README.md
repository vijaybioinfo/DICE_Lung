# Preprocessing of reads by mapping them to the reference genome and collapsing them into UMI count matrices (10x's cellranger)

This folder has been created for internal purposes only, mainly for the review of this section to be performed by Kevin.
Major things to check:
- Correctness of the CITE-seq sample sheets. I've left two examples, for <code>sample-id-1</code> and <code>sample-id-2</code> within sample project <code>sample_project_1</code>.
- Run of at least one gene expression sample per sequencing run. I've left the path to an example for sample <code>sample-id-1_Gex</code>
- Recreate aggr tables. You can find my aggr tables in the following path: <code>/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/sequencing_data/05-14-2022/aggr/data</code>

### sequencing_data/

* <code>date_1/</code>
  * <code>aggr/</code>
  * <code>count/</code>
    * <code>data/</code>
        * <code>experiment_layout/</code>
            * <code>sample_project_1/</code>
                * <code>**CITE-seq_sample_sheet_for_feat_barcoding_NV0XX_sample-id-1_CITE.csv**</code>
                * <code>**CITE-seq_sample_sheet_for_feat_barcoding_NV0XX_sample-id-2_CITE.csv**</code>
            * <code>sample_project_2/</code>
                * </code>...</code>
    * <code>sample_project_1/</code>
    * <code>sample_project_2/</code>
  * <code>mkfastq/</code>
  * <code>jobs_scripts/</code>
    * <code>aggr/</code>
    * <code>count/</code>
        * <code>sample_project_1/</code>
            * <code>sample-id-1_Gex</code>
                * <code>**sample-id-1_Gex_cellranger_count.job.sh**</code>
            * <code>sample-id-1_CITE</code>
            * <code>sample-id-2_Gex</code>
            * <code>sample-id-2_CITE</code>
        * <code>sample_project_2/</code>
            * <code>...</code>
    * <code>mkfastq/</code>
* <code>date_2/</code>
  * <code>aggr/</code>
  * <code>count/</code>
  * <code>mkfastq/</code>
* <code>date_n/</code>
    * <code>...</code>