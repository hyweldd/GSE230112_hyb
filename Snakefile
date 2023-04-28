# -*- coding: utf-8; mode: python; tab-width: 4; -*-
# vim: ft=python

# ______________________________________________________________________________
#
#     Snakefile
#     Copyright (C) 2023  Hywel.Dunn-Davies@ed.ac.uk
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Affero General Public License as
#     published by the Free Software Foundation, either version 3 of the
#     License, or (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Affero General Public License for more details.
#
#     You should have received a copy of the GNU Affero General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ______________________________________________________________________________

# Snakemake pipeline to generate the hyb files in GEO accession GSE230112 from 
# the corresponding source data.

import os
import re


# Constants ----

HYB_COMMAND = "hyb align=bowtie2 word=11 format=fasta type=all pref=none"
HYB_DB = "meg3_hyb_db"


# Setup ----
r1_prefixes = glob_wildcards("00_source_fastq_gz/{prefix}_1.fq.gz").prefix
r2_prefixes = glob_wildcards("00_source_fastq_gz/{prefix}_2.fq.gz").prefix

singularity: 'ezh2_flash_hyb_analysis_pipeline.sif'

# Snakemake rules ----

# Target rule ---
rule target:
    input:
        expand("01_cutadapt/{prefix}_trimmed_{read}.fq", prefix=r1_prefixes, read = [1,2]),
        expand("02_collapsed_r1_fasta/{prefix}_trimmed_1_comp.fasta", prefix=r1_prefixes),
        expand("03_hyb_alignment_results/{prefix}_trimmed_1_comp_meg3_hyb_db_mtophits.blast", prefix=r1_prefixes),
        "04_combined_ref_file/combined_meg3.ref",
        expand("05_hyb_analysis_results/{prefix}_trimmed_1_comp_meg3_hyb_db_hybrids_ua_dg.hyb", prefix=r1_prefixes)


# Pipeline rules ---
rule trim_adapters:
    input:
        r1 = "00_source_fastq_gz/{r1_prefix}_1.fq.gz",
        r2 = "00_source_fastq_gz/{r1_prefix}_2.fq.gz"
    output:
        r1 = "01_cutadapt/{r1_prefix}_trimmed_1.fq",
        r2 = "01_cutadapt/{r1_prefix}_trimmed_2.fq"
    log:
        out = "01_cutadapt/{r1_prefix}_trimmed.fq.out",
        err = "01_cutadapt/{r1_prefix}_trimmed.fq.err"
    threads:
        5
    shell:
        r"""
        cutadapt -j {threads} -e 0.1 -q 16 -O 3 --trim-n --minimum-length 15 -a TGGAATTCTCGGGTGCCAAGGC -A GCCTTGGCACCCGAGAATTCCA -o {output.r1} -p {output.r2} {input.r1} {input.r2} > {log.out} 2> {log.err}
        """

rule collapse_duplicates:
    input:
        "01_cutadapt/{r1_prefix}_trimmed_1.fq"
    output:
        "02_collapsed_r1_fasta/{r1_prefix}_trimmed_1_comp.fasta"
    log:
        out = "02_collapsed_r1_fasta/{r1_prefix}_trimmed_1_comp.fasta.out",
        err = "02_collapsed_r1_fasta/{r1_prefix}_trimmed_1_comp.fasta.err"
    shell:
        r"""
        pyFastqDuplicateRemover.py -f {input} -o {output} > {log.out} 2> {log.err}
        """

rule align:
    input:
        "02_collapsed_r1_fasta/{prefix}_trimmed_1_comp.fasta"
    output:
        tab = "03_hyb_alignment_results/{prefix}_trimmed_1_clipped_qf.tab",
        fasta = "03_hyb_alignment_results/{prefix}_trimmed_1_comp.fasta",
        blast = "03_hyb_alignment_results/{prefix}_trimmed_1_comp_meg3_hyb_db.blast",
        mtophits_blast = "03_hyb_alignment_results/{prefix}_trimmed_1_comp_meg3_hyb_db_mtophits.blast"
    log:
        out = "03_hyb_alignment_results/{prefix}.out",
        err = "03_hyb_alignment_results/{prefix}.err"
    threads:
        28
    params:
        script = HYB_COMMAND,
        db = HYB_DB
    shell:
        r"""
        cd 03_hyb_alignment_results
        {params.script} id={wildcards.prefix} in=../{input} db={params.db} `basename {output.mtophits_blast}` `basename {output.tab}` `basename {output.fasta}` `basename {output.blast}` > ../{log.out} 2> ../{log.err}
        """

rule create_combined_ref_file:
    input:
        expand("03_hyb_alignment_results/{prefix}_trimmed_1_comp_meg3_hyb_db_mtophits.blast", prefix=r1_prefixes)
    output:
        ref = "04_combined_ref_file/combined_meg3.ref",
        blast = "04_combined_ref_file/combined_meg3.blast"
    log:
        err = "04_combined_ref_file/combined_meg3.ref.err"
    params:
        script = "LANG=C create_reference_file.pl"
    shell:
        r"""
        cat {input} > {output.blast}
        {params.script} {output.blast} > {output.ref} 2> {log.err}
        """

rule hyb_analyse:
    input:
        fasta = "02_collapsed_r1_fasta/{prefix}_trimmed_1_comp.fasta",
        blast = "03_hyb_alignment_results/{prefix}_trimmed_1_comp_meg3_hyb_db.blast",
        mtophits_blast = "03_hyb_alignment_results/{prefix}_trimmed_1_comp_meg3_hyb_db_mtophits.blast",
        ref = "04_combined_ref_file/combined_meg3.ref"
    output:
        ua_dg_hyb = "05_hyb_analysis_results/{prefix}_trimmed_1_comp_meg3_hyb_db_hybrids_ua_dg.hyb",
        merged_hyb = "05_hyb_analysis_results/{prefix}_trimmed_1_comp_meg3_hyb_db_hybrids_ua_merged.hyb",
        hyb_stats = "05_hyb_analysis_results/{prefix}_trimmed_1_comp_meg3_hyb_db_hybrids.hyb_stats.txt"
    log:
        out = "05_hyb_analysis_results/{prefix}.out",
        err = "05_hyb_analysis_results/{prefix}.err"
    params:
        script = HYB_COMMAND,
        db = HYB_DB
    shell:
        r"""
        cd 05_hyb_analysis_results

        ln -s ../{input.fasta}
        ln -s ../{input.blast}
        ln -s ../{input.mtophits_blast}

        {params.script} id={wildcards.prefix} in=../{input.fasta} db={params.db} ref=../{input.ref} analyse > ../{log.out} 2> ../{log.err}
        """

