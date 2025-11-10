import argparse
import csv
import gzip
import re

from pybedtools import BedTool


def merge_exome(gtf_fp: str, out_fp: str):
    """Merges overlapping exon features within genes.

    Args:
        gtf_fp: The `.gtf` file path.
        out_fp: The output `.csv` file path.

    Writes:
        A `.csv` with merged exon, i.e., meta-exon, information.
    """
    exome = read_exome(gtf_fp)
    exome_merged = exome.merge(c=[4, 5, 6], o='distinct')

    with open(out_fp, 'w') as f:
        writer = csv.writer(f)
        writer.writerow([
            "gene_id",
            "gene_biotype",
            "chrom",
            "meta_beg",
            "meta_end",
            "strand",
            "meta_type",
            "make_up"
        ])

        for row in exome_merged:
            gene_func, meta_beg, meta_end, make_up, chrom, strand = row
            gene, func = gene_func.split("=")

            meta_type = "_".join(
                sorted(set(re.findall(r'\w{2}(?==)', make_up)))
            )

            writer.writerow([
                gene,
                func,
                chrom,
                meta_beg,
                meta_end,
                strand,
                meta_type,
                make_up
            ])


def read_exome(gtf_gz: str):
    """Reads exome from a `.gtf` file.

    Args:
        gtf_gz: The `.gtf.gz` file path.

    Returns:
        A `BedTool` instance, where field 1 stores the gene and field 5 the
        chromosome name — the idea being to prevent between-gene merging
        without having run a separate merge operation on each gene, which is
        much slower.
    """
    exome = str()

    with gzip.open(gtf_gz, 'rt') as f:
        tran_beg = None
        tran_end = None

        for row in f:
            if row.startswith("#"):
                continue

            row = row.strip("\n").split("\t")
            chrom, _, feat, *pos, _, strand, _, meta_data = row

            beg, end = map(int, pos)

            if feat == "transcript":
                tran_beg = beg
                tran_end = end

                continue

            if feat == "exon":
                gene = get_tag('gene_id', meta_data)
                func = get_tag('gene_biotype', meta_data)

                boundary_match = int(beg == tran_beg), int(end == tran_end)

                if boundary_match == (1, 1):
                    exon_type = "SE"

                if boundary_match == (1, 0):
                    exon_type = "FE" if strand == "+" else "LE"

                if boundary_match == (0, 1):
                    exon_type = "LE" if strand == "+" else "FE"

                if boundary_match == (0, 0):
                    exon_type = "IE"

                exome += (
                    f"{gene}={func}\t"
                    f"{beg}\t"
                    f"{end}\t"
                    f"{exon_type}={beg}-{end}\t"
                    f"{chrom}\t"
                    f"{strand}\n"
                )

    return BedTool(exome, from_string=True).sort()


def get_tag(tag: str, meta_data: str):
    """Parses the `attribute` field of a `.gtf` file.

    Args:
        tag: The name of an attribute tag in a `.gtf` file.
        meta_data: The meta-data, i.e., attribute string, to parse.

    Returns:
        The value for `tag`, or `None` when `tag` cannot be found.
    """
    match = re.search(rf'(?<={tag} ")\w+', meta_data)

    if match:
        return match.group(0)
    else:
        return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-g', help="The <.gtf> file path.")
    parser.add_argument('-o', help="The <.csv> path to write output to.")

    args = parser.parse_args()
    merge_exome(args.g, args.o)
