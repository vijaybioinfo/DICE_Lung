rule merge_cells:
    input:
        get_overlaps
    output:
        CFG.prefix + "/overlap/{project}/{disease}/{dataset}.tsv"
    resources:
        mem_gb = 20,
        walltime = 60
    run:
        from pandas import concat
        from R24.Overlap.module.utils import (read_file, 
                                               get_filename, 
                                               trim_file_extension,
                                               Columns)

        combined = []
        for file in input:
            data = read_file(file)
            data[Columns.CELL] = trim_file_extension(get_filename(file))
            combined.append(data)

        combined = concat(combined)
        combined.to_csv(output[0], index=False, sep="\t")

rule merge_datasets:
    input:
        merge_by_category
    output:
        merge = CFG.prefix + "/summary/{project}/{disease}/{disease}_{category}.tsv",
        matrix = CFG.prefix + "/summary/{project}/{disease}/{disease}_{category}.mtx"
    params:
        metadata = QTL.get_metadata,
        prefix = CFG.prefix + "/summary/{project}/{disease}/{disease}_{category}"
    resources:
        mem_gb = 20,
        walltime = 30
    script:
        "../scripts/summary.py"