#!/usr/bin/env python

import argparse
import sys

def smartopen(fname):
    if fname == '-' : return sys.stdin
    else: return open(fname)

def smartclose(infile):
    if infile == sys.stdin : return None
    else: return infile.close()

def main():
    parser = argparse.ArgumentParser("Evaluate delins variant calls.")
    parser.add_argument('--base-vcf', type=str, required=True,
            help="Baseline VCF file containing true positive variants. ")
    parser.add_argument('--call-vcf', type=str, required=True,
            help="Call VCF file containing variant candidates predicted to be positive. ")
    parser.add_argument('--field', type=str,
            help="The VCF tag field used to rank variant candidates.", default = 'QUAL')
    args = parser.parse_args()
    
    #print(args)
    #print(args.base_vcf)
    #print(args.call_vcf)
    base_variants = set([])
    base_samples = []
    basefile = smartopen(args.base_vcf) #(open(args.base_vcf) if (args.base_vcf != '-') else sys.stdin)
    for line in basefile:
        if line.startswith('##'): continue
        if line.startswith('#CHROM'):
            for tokens in line.rstrip().split('\t'):
                base_samples = tokens[9:]
            continue
        tokens = line.rstrip().split('\t')
        assert len(tokens) > 4, 'The tokens ( {} ) are invalid'.format(tokens)
        chrom, pos, ref, alt = (tokens[0], int(tokens[1]), tokens[3], tokens[4])
        for sample, format in zip(base_samples, tokens[9:]):
            GT = format.split(':')[0]
            if '1' in GT: base_variants.add((sample, chrom, pos, ref, alt))
        if len(base_samples) == 0: base_variants.add(('-', chrom, pos, ref, alt))
        #sys.stderr.write('base_var = {}\n'.format(tokens))
    smartclose(basefile)
    n_base_snv = 0
    n_base_indel = 0
    n_base_delins = 0
    n_base_error = 0
    for base_variant in base_variants:
        ref = base_variant[3]
        alt = base_variant[4]
        if len(ref) == 1 and len(alt) == 1: n_base_snv += 1
        else:
            if (ref[0] == alt[0]) and (len(ref) == 1 or len(alt) ==1):
                n_base_indel += 1
            elif len(ref) >= 1 and len(alt) >= 1:
                n_base_delins += 1
            else: 
                n_base_error += 1
    call_variants = []
    call_samples = []
    with open(args.call_vcf) as callfile:
        for line in callfile:
            if line.startswith('##'): continue
            if line.startswith('#CHROM'):
                for tokens in line.rstrip().split('\t'):
                    call_samples = tokens[9:]
                continue
            tokens = line.rstrip().split('\t')
            chrom, pos, ref, alt = (tokens[0], int(tokens[1]), tokens[3], tokens[4])
            subfields = args.field.split('/')
            if args.field == 'QUAL':
                field_idx = -2
                varqual = float(tokens[5])
            elif len(subfields) == 2 and subfields[0] == 'INFO':
                field_idx = -1
                for infotoken in tokens[7].split(';'):
                    k = infotoken.split('=')[0]
                    if k == subfields[1]:
                        varqual = float(infotoken.split('=')[1])
            elif len(subfields) == 2 and subfields[0] == 'FORMAT':
                field_idx = tokens[8].split(':').index(subfields[1])
            for format in tokens[9:]:
                GT = format.split(':')[0]
                if field_idx >= 0: varqual = float(format.split(':')[field_idx].split(',')[-1])
                sample = ('-' if len(base_samples) == 0 else format)
                isTP = (1 if (sample, chrom, pos, ref, alt) in base_variants else 0)
                call_variants.append((sample, chrom, pos, ref, alt, varqual, isTP))
            if not tokens[9:]:
                isTP = (1 if ('-', chrom, pos, ref, alt) in base_variants else 0)
                call_variants.append(('-', chrom, pos, ref, alt, varqual, isTP))
    call_variants2 = sorted(call_variants, key = lambda x:-x[5])
    tp_num = 0
    tot_num = 0
    tot_true_tp_num = len(base_variants)
    max_fscore = -1.0;
    lines = []
    lines.append('##num_base_vars={},num_call_vars={},n_base_snv={},n_base_indel={},n_base_delins={},n_base_error={}'.format(
            len(base_variants), len(call_variants),
            n_base_snv, n_base_indel, n_base_delins, n_base_error))
    lines.append('#' + '\t'.join(['sample', 'chrom', 'pos', 'ref', 'alt', 'qual', 'isTP', 'tpNum', 'fpNum', 'fnNum', 'precision', 'sensitivity', 'fscore']))
    for call_variant in call_variants2:
        tp_num += call_variant[-1]
        tot_num += 1
        assert (tot_num > 0)
        prec = tp_num / float(tot_num)
        sens = (tp_num / float(tot_true_tp_num) if (tot_true_tp_num > 0) else float('nan'))
        fscore = ((2.0 / (prec**(-1) + sens**(-1))) if (prec > 0 and sens > 0) else 0)
        if fscore > max_fscore: 
            max_fscore = fscore
        fp_num = tot_num - tp_num
        fn_num = tot_true_tp_num - tp_num
        lines.append('\t'.join([str(x) for x in call_variant]) + '\t' + '\t'.join([str(x) for x in [tp_num, fp_num, fn_num, prec, sens, fscore]]))
    if tot_true_tp_num == 0:
        max_fscore_str = 'is_not_available'
    elif max_fscore > 1-1e-6: 
        max_fscore_str = 'is_perfect'
    elif max_fscore < 1e-6:
        max_fscore_str = 'is_false_negative'
    else:
        max_fscore_str = 'is_between_0_and_1'
    print('##max_fscore={} max_fscore_type={}'.format(max_fscore, max_fscore_str))
    print('\n'.join(lines))
if __name__ == '__main__':
    main()

