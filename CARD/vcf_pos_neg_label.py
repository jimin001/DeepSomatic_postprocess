import pysam
import argparse

"""
running command:
python3 vcf_pos_neg_label.py \
    --positive TP.vcf.gz \
    --negative FP.vcf.gz \
    --outfile output_labeled.vcf
"""

parser = argparse.ArgumentParser()
parser.add_argument('--positive', '-p', type=str, required=True, help='vcf file path of positive variants (TP)')
parser.add_argument('--positivevcfindex', '-q', type=str, help='vcf index file path of positive variants (TP)')
parser.add_argument('--negative', '-n', type=str, required=True, help='vcf file path of negative variants (FP)')
parser.add_argument('--negativevcfindex', '-m', type=str, help='vcf index file path of negative variants (FP)')
parser.add_argument('--outfile', '-o', type=str, required=True, help='output file path')

args = parser.parse_args()


class VCFHandler:
    """
    VCF Handler class.
    """
    def __init__(self, input_vcf_path, input_vcf_index_path=None):
        self.vcf = pysam.VariantFile(input_vcf_path, 'r', index_filename=input_vcf_index_path)
        self.header = self.vcf.header

    def get_records(self):
        return self.vcf.fetch()

    def reformat_record(self, record, label):
        """
        Reformat a single record into a tab-separated VCF string.
        Returns None for records that should be skipped (e.g. RefCall).
        :param record: pysam VariantRecord object
        :param label: 'positive' (TP) or 'negative' (FP)
        """
        filter_key = list(record.filter.keys())[0] if record.filter.keys() else '.'

        if filter_key == 'RefCall':
            return None

        alt_str = record.alts[0] if record.alts else '.'
        qual_str = str(round(record.qual, 2)) if record.qual is not None else '.'
        filter_str = filter_key if filter_key else '.'

        return f'{record.chrom}\t{record.pos}\t.\t{record.ref}\t{alt_str}\t{qual_str}\t{filter_str}\ttype={label}\n'


TYPE_INFO_LINE = '##INFO=<ID=type,Number=1,Type=String,Description="Variant label: positive (TP) or negative (FP)">'


def build_header(header):
    """
    Return a filtered header string containing only lines relevant to the
    output VCF (no ##FORMAT or ##INFO lines from the input, since the output
    records contain only CHROM/POS/ID/REF/ALT/QUAL/FILTER/INFO=type).
    Inserts TYPE_INFO_LINE before #CHROM.
    """
    lines = []
    for line in str(header).splitlines():
        if line.startswith('##FORMAT=') or line.startswith('##INFO='):
            continue
        lines.append(line)

    chrom_idx = next(i for i, l in enumerate(lines) if l.startswith('#CHROM'))
    lines[chrom_idx] = '\t'.join(lines[chrom_idx].split('\t')[:8])
    lines.insert(chrom_idx, TYPE_INFO_LINE)

    return '\n'.join(lines) + '\n'


def chrom_sort_key(chrom):
    """
    Sort key that orders chromosomes numerically (chr1 < chr2 < ... < chr10)
    followed by non-numeric contigs (chrX, chrY, chrM, etc.) lexicographically.
    """
    c = chrom.lower().lstrip('chr')
    try:
        return (0, int(c), '')
    except ValueError:
        return (1, 0, c)


if __name__ == '__main__':
    tp_vcf = VCFHandler(args.positive, args.positivevcfindex)
    fp_vcf = VCFHandler(args.negative, args.negativevcfindex)

    all_records = (
        [(record, 'positive') for record in tp_vcf.get_records()] +
        [(record, 'negative') for record in fp_vcf.get_records()]
    )
    all_records.sort(key=lambda x: (chrom_sort_key(x[0].chrom), x[0].pos))

    with open(args.outfile, 'w') as out:
        out.write(build_header(tp_vcf.header))

        for record, label in all_records:
            line = tp_vcf.reformat_record(record, label)
            if line is not None:
                out.write(line)
