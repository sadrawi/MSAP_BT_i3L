from pathlib import Path
import random
rng=random.Random(20261016)
out=Path(__file__).parent/'synthetic_demo.fasta'
seq=[]
# Create chain of alternating restriction sites with different nearby selective bases.
for i in range(120):
    motif='GAATTC' if i%2==0 else 'CCGG'
    seq.append(motif)
    seq.append(''.join(rng.choices('ACGT',k=rng.randint(130,310))))
s=''.join(seq)
with out.open('w') as w:
    w.write('>synthetic_demo_not_rice\n')
    for i in range(0,len(s),80):w.write(s[i:i+80]+'\n')
print('Wrote',out,len(s),'bp')
