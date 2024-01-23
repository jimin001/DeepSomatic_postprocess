import pysam
import argparse

"""
running command
VCF=/Users/jimin/PycharmProjects/PEPPER/VCF/1108_1395_Ilmn_DeepSomatic_v1.6.0.vcf.gz
INDEX_FILE=/Users/jimin/PycharmProjects/PEPPER/VCF/1108_1395_Ilmn_DeepSomatic_v1.6.0.vcf.gz.tbi
OUTPUT=/Users/jimin/PycharmProjects/PEPPER/VCF/tester_strings.vcf
HEAD=/Users/jimin/PycharmProjects/PEPPER/VCF/tester_header.vcf
python3 somatic_vcf_training_format.py -v $VCF -i $INDEX_FILE -o $OUTPUT -x $HEAD

# post-process steps
cd $VCF_OUTPUT_PATH
head -n -1 $HEAD > temp_delete_last_line_header.vcf
cat temp_delete_last_line_header.vcf $OUTPUT > final_training.vcf
"""



parser = argparse.ArgumentParser()
parser.add_argument('--vcf', '-v', type=str, required=True, help='vcf file path')
parser.add_argument('--vcfindex', '-i', type=str, help='vcf index file path')
parser.add_argument('--outfile', '-o', type=str, help='output file path')
parser.add_argument('--outfile_header', '-x', type=str, help='output file path')

args = parser.parse_args()

class DeepSomaticSampleData:
    """
    Store information for GT:AP:GQ:DP:AD:VAF:REP.
    """
    def __init__(self, samples):
        self.samples = samples

        for sample in samples:
            self.sample = sample
            self.GT = self.sample[0]
            self.GQ = self.sample[1]
            self.DP = self.sample[2]
            self.AD = self.sample[3]
            self.VAF = self.sample[4]
            self.PL = self.sample[5]

class VCFHandler:
    """
    VCF Handler class.
    """
    def __init__(self, input_vcf_path, input_vcf_index_path):
        self.vcf = pysam.VariantFile(input_vcf_path, 'r', index_filename=input_vcf_index_path)
        self.header = self.vcf.header
        self.contig = None

    def get_header(self):
        """
        Get header of VCF.
        """
        return self.header

    def get_records(self):
        """
        Read in vcf line by line in specified contig region
        """
        # returns iterator object
        return self.vcf.fetch(self.contig)

    def get_record_attributes(self, records, attribute):
        """
        Get specified feature of vcf record.
        """
        # returns generator object of specified attribute
        for record in records:
            if attribute == 'chrom':
                yield record.chrom
            elif attribute == 'pos':
                yield record.pos
            elif attribute == 'id':
                yield record.id
            elif attribute == 'ref':
                yield record.ref
            elif attribute == 'alt':
                yield record.alts
            elif attribute == 'qual':
                yield record.qual
            elif attribute == 'filter':
                yield record.filter.keys()
            elif attribute == 'info':
                yield record.info.keys()
            # elif attribute == 'format':
                # yield record.format.keys()
            # elif attribute == 'sample':
                # yield [y.values() for y in record.samples.itervalues()]

    def get_record_samples(self, samples_list, type):
        """
        Return generator object of specified data type of sample.
        :param samples_list: samples attribute generator object or list
        :param type: GT:AP:GQ:DP:AD:VAF:REP
        :return:
        """
        for s in samples_list:
            samples = SampleData(s)
            if type == 'GT':
                yield samples.GT
            elif type == 'GQ':
                yield samples.GQ
            elif type == 'DP':
                yield samples.DP
            elif type == 'AD':
                yield samples.AD
            elif type == 'VAF':
                yield samples.VAF
            elif type == 'PL':
                yield samples.PL
            else:
                return 'Incorrect sample type'

    def filter_min_vaf(self, records, min_vaf):
        for record in records:
            sample_info = [y.values() for y in record.samples.itervalues()]
            #print(sample_info)
            sample_obj = SampleData(sample_info)
            vaf = sample_obj.VAF[0]
            #print(vaf)
            if vaf >= min_vaf:
                yield record

    def reformat_records(self, records):
        for record in records:

            chrom = record.chrom
            pos = record.pos
            id = record.id
            ref = record.ref
            alt = record.alts
            qual = record.qual

            reformat_string = chrom + '\t' + str(pos) + '\t' + '.' + '\t' + ref + '\t' + alt[0] + '\t' + str(round(qual, 2))

            if record.filter.keys()[0] == 'RefCall':
                continue
            elif record.filter.keys()[0] == 'PASS':
                reformat_string += '\t' + 'PASS' + '\t' + 'type=somatic' + '\n'
                
                yield reformat_string
            elif record.filter.keys()[0] == 'GERMLINE':
                reformat_string += '\t' + 'PASS' + '\t' + 'type=germline' + '\n'
               
                yield reformat_string

if __name__ == '__main__':
    vcf_file = args.vcf
    vcf_index = args.vcfindex
    out_file_name = args.outfile_header
    out_file_strings = args.outfile
    vcf_h = VCFHandler(vcf_file, vcf_index)

    out_file = pysam.VariantFile(out_file_name, 'w', header=vcf_h.header)

    records = vcf_h.get_records()
    reformatted_records = vcf_h.reformat_records(records)

    out = open(out_file_strings, 'w')
    header_strings = "#CHROM" + '\t' + "POS" + '\t' + "ID" + '\t' + "REF" + '\t' + "ALT" + '\t' + "QUAL" + '\t' \
                     + "FILTER" + '\t' + "INFO" + '\n'
    out.write(header_strings)
    for record in reformatted_records:
        out.write(record)

