#!/usr/bin/env python3
"""Sequence-only EcoRI--HpaII/MspI potential-fragment simulator (unmethylated model).

Assumptions: complete digest, linear contigs, EcoRI G^AATTC and C^CGG,
no methylation model, no PCR efficiencies. Selective extensions are read
inward from each recognition site in the adapter-primer's 5'->3' direction.
"""
import argparse
import csv
import gzip
import itertools
from collections import defaultdict
from pathlib import Path

SITES = {"E": "GAATTC", "H": "CCGG"}
RC_TRANS = str.maketrans("ACGT", "TGCA")

def rc(s): return s.translate(RC_TRANS)[::-1]

def fasta_records(path):
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rt') as f:
        name, chunks = None, []
        for line in f:
            if line.startswith('>'):
                if name is not None: yield name, ''.join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else: chunks.append(line.strip())
        if name is not None: yield name, ''.join(chunks).upper()

def cut_sites(seq):
    sites=[]
    for enzyme, motif in SITES.items():
        at=seq.find(motif)
        while at>=0:
            sites.append((at+1, at, enzyme))  # recognition cut offset +1
            at=seq.find(motif, at+1)
    return sorted(sites)

def candidates(name, seq, minimum=100, maximum=500, max_n=3):
    """Generate EcoRI-CCGG fragments bounded by consecutive digest sites.
    CCGG sites act identically to HpaII/MspI only in unmethylated model.
    All cuts, including same-enzyme cuts, break molecules.
    """
    sites=cut_sites(seq)
    for (cut_l,start_l,kind_l),(cut_r,start_r,kind_r) in zip(sites, sites[1:]):
        if kind_l==kind_r: continue
        fragment_size=cut_r-cut_l  # biological length (without adapters)
        if not minimum <= fragment_size <= maximum: continue
        if kind_l=='E':
            e_ext=seq[start_l+6:start_l+6+max_n]
            h_ext=rc(seq[start_r-max_n:start_r])
        else:
            h_ext=seq[start_l+4:start_l+4+max_n]
            e_ext=rc(seq[start_r-max_n:start_r])
        if len(e_ext)<max_n or len(h_ext)<max_n: continue
        if any(x not in 'ACGT' for x in e_ext+h_ext): continue
        yield {'contig':name,'start_1based':cut_l+1,'end_1based':cut_r,
               'bp':fragment_size,'E_extension':e_ext,'HM_extension':h_ext,
               'orientation':kind_l+'-'+kind_r}

def analyze(fasta, minimum, maximum, e_len, h_len, bin_bp, output, max_n=3):
    if max(e_len,h_len)>max_n: raise ValueError('extension lengths exceed max_n')
    output=Path(output); output.mkdir(parents=True, exist_ok=True)
    pairs=defaultdict(list)
    total=0
    fragment_file=output/'predicted_fragments.csv'
    with fragment_file.open('w', newline='') as out:
        writer=csv.DictWriter(out, fieldnames=['contig','start_1based','end_1based','bp','E_extension','HM_extension','orientation'])
        writer.writeheader()
        for contig,seq in fasta_records(fasta):
            for fragment in candidates(contig,seq,minimum,maximum,max_n):
                writer.writerow(fragment)
                pairs[(fragment['E_extension'][:e_len], fragment['HM_extension'][:h_len])].append(fragment)
                total+=1
    # explicit binning is an approximate measure of gel peak crowding, not actual resolution
    summary=[]
    for e,h in itertools.product(map(''.join,itertools.product('ACGT',repeat=e_len)),map(''.join,itertools.product('ACGT',repeat=h_len))):
        frags=pairs.get((e,h),[])
        bins={int((f['bp']-minimum)//bin_bp) for f in frags}
        summary.append({'E_extension':e,'HM_extension':h,'fragments':len(frags),
            'occupied_size_bins':len(bins),'co_migrating_excess':len(frags)-len(bins),
            'unique_contigs':len({f['contig'] for f in frags}),
            'min_bp':min([f['bp'] for f in frags],default=''),
            'max_bp':max([f['bp'] for f in frags],default='')})
    # rank by occupied peaks, penalize crowding (heuristic, NOT empirically calibrated)
    summary.sort(key=lambda x:(x['occupied_size_bins']-2*x['co_migrating_excess'],x['occupied_size_bins']),reverse=True)
    with (output/'primer_ranking.csv').open('w',newline='') as out:
        w=csv.DictWriter(out,fieldnames=summary[0].keys());w.writeheader();w.writerows(summary)
    print(f'Sequence-only potential fragments in length range: {total}')
    print(f'Primer pairs evaluated: {len(summary)}')
    print(f'Wrote: {fragment_file}, {output / "primer_ranking.csv"}')
    print('Top 5 combinations (unvalidated heuristic):')
    for row in summary[:5]: print(row)
    return summary

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fasta',required=True,help='Plain or .gz FASTA')
    p.add_argument('--out',default='results')
    p.add_argument('--min-bp',type=int,default=100)
    p.add_argument('--max-bp',type=int,default=500)
    p.add_argument('--e-bases',type=int,default=2,choices=[1,2,3])
    p.add_argument('--hm-bases',type=int,default=2,choices=[1,2,3])
    p.add_argument('--bin-bp',type=int,default=5,help='Approximate gel resolution for ranking')
    args=p.parse_args()
    if args.min_bp<=0 or args.min_bp>args.max_bp or args.bin_bp<=0: p.error('Invalid size range or bin')
    analyze(args.fasta,args.min_bp,args.max_bp,args.e_bases,args.hm_bases,args.bin_bp,args.out)
if __name__=='__main__':main()
