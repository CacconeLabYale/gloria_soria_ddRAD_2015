import csv

from collections import defaultdict

import pandas as pd
pd.set_option('display.max_columns', 60)

import munch

from IPython.display import display, HTML

import pybedtools as pbt

import sh

import tabulate
tbl = tabulate.tabulate

from spartan.utils.annotations.ensembl.gff3 import parse_gff3
from spartan.utils.annotations.ensembl.gff3 import parse_gff3_attributes
from spartan.utils.files import tableFile2namedTuple

from spartan.utils.genome_specific.GfusI1 import GfusI1_0


def make_BED_line(in_line, name_map):
    chrom = name_map[in_line[1]]
    chromstart = str(int(in_line[2]) - 1)
    chromend = str(int(in_line[2]))
    
    return chrom, chromstart, chromend

def get_in_lines(path, skip_first_row=True):
    with open(path, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        if skip_first_row:
            reader.next()
        for line in reader:
            yield line
            
def write_line(out_file, line):
    line = "%s\n" % ('\t'.join(line))
    out_file.write(line)

def convert_a_file(in_path, out_path, name_map):
    with open(out_path, 'wb') as out_file:
        
        in_lines = get_in_lines(in_path)

        for line in in_lines:
            bed_data = make_BED_line(line, name_map)
            write_line(out_file,bed_data)

def snp_vs_gff_to_DF(bedtools_out):
    headers = ["bed3_seq",
               "bed3_start",
               "bed3_end",
               "gff3_seq",
               "gff3_source",
               "gff3_type",
               "gff3_start",
               "gff3_end",
               "gff3_score",
               "gff3_strand",
               "gff3_phase",
               "gff3_attributes",]
    df = pd.read_csv(bedtools_out, sep='\t', names=headers)
    
    gene_id = lambda x: parse_gff3_attributes(x)['ID']
    
    df['gff3_rec'] = df.gff3_attributes.apply(gene_id)
    
    return df

def answer_set_1(df, title="None"):
    p = df.query("Aspect == 'P'")
    f = df.query("Aspect == 'F'")
    c = df.query("Aspect == 'C'")

    num_genes = len(df.gene_id.unique())
    num_contigs = len(df.bed3_seq.unique())
    num_func_annos = len(df.Name.unique())
    num_p = len(p.Name.unique())
    num_f = len(f.Name.unique())
    num_c = len(c.Name.unique())

    anno_stats = tbl(pd.DataFrame(df['Total Score'].describe()), headers=('Metric','Value'))
    anno_stats_p = tbl(pd.DataFrame(p['Total Score'].describe()), headers=('Metric','Value'))
    anno_stats_f = tbl(pd.DataFrame(f['Total Score'].describe()), headers=('Metric','Value'))
    anno_stats_c = tbl(pd.DataFrame(c['Total Score'].describe()), headers=('Metric','Value'))
    
    top10_annos = tbl(pd.DataFrame(df.Name.value_counts().head(10)), headers=('Annotation','Genes'))
    top10_annos_p= tbl(pd.DataFrame(p.Name.value_counts().head(10)), headers=('Annotation','Genes'))
    top10_annos_f= tbl(pd.DataFrame(f.Name.value_counts().head(10)), headers=('Annotation','Genes'))
    top10_annos_c= tbl(pd.DataFrame(c.Name.value_counts().head(10)), headers=('Annotation','Genes'))

    top10_annos_TS_800 = tbl(pd.DataFrame(df[df['Total Score'] >= 800].Name.value_counts().head(10)), headers=('Annotation','Genes'))
    top10_annos_TS_800_p = tbl(pd.DataFrame(p[p['Total Score'] >= 800].Name.value_counts().head(10)), headers=('Annotation','Genes'))
    top10_annos_TS_800_f = tbl(pd.DataFrame(f[f['Total Score'] >= 800].Name.value_counts().head(10)), headers=('Annotation','Genes'))
    top10_annos_TS_800_c = tbl(pd.DataFrame(c[c['Total Score'] >= 800].Name.value_counts().head(10)), headers=('Annotation','Genes'))


    return """
## {title} ##


- __Genes:__ {num_genes}
- __Contigs:__ {num_contigs}
- __Unique Terms (all):__ {num_func_annos}




### Annotation Scores ###

Table: All domains

{anno_stats}

Table: Process

{anno_stats_p}

Table: Function

{anno_stats_f}

Table: Cellular

{anno_stats_c}

----

### Top annotations ###

Table: All domains

{top10_annos}

Table: Process

{top10_annos_p}

Table: Function

{top10_annos_f}

Table: Cellular

{top10_annos_c}


----

### Top annotations (TS >= 800) ###

Table: All domains

{top10_annos_TS_800}

Table: Process

{top10_annos_TS_800_p}

Table: Function

{top10_annos_TS_800_f}

Table: Cellular

{top10_annos_TS_800_c}


\\newpage

""".format(
    title=title,
    num_genes=num_genes,
    num_contigs=num_contigs,
    num_func_annos=num_func_annos,
    anno_stats=anno_stats,
    anno_stats_p=anno_stats_p,
    anno_stats_f=anno_stats_f,
    anno_stats_c=anno_stats_c,
    top10_annos=top10_annos,
    top10_annos_p=top10_annos_p,
    top10_annos_f=top10_annos_f,
    top10_annos_c=top10_annos_c,
    top10_annos_TS_800=top10_annos_TS_800,
    top10_annos_TS_800_p=top10_annos_TS_800_p,
    top10_annos_TS_800_f=top10_annos_TS_800_f,
    top10_annos_TS_800_c=top10_annos_TS_800_c,
    )