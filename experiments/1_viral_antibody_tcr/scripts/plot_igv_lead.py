#!/usr/bin/env python3
"""Plot measured genome context and exact local alignment; no inferred allele call."""
import hashlib
import json
import os
import sys
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/bio-fishing-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Patch, Rectangle

from acquire import ROOT, sha256, utc


def main():
    output = ROOT / 'results/metavr/lead_followup'
    source = output / 'figure_data.json'
    data = json.loads(source.read_text())
    figures = ROOT / 'figures'
    figures.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'svg.hashsalt': 'igv-lead-20261116', 'pdf.fonttype': 42})
    ink, viral, candidate, gray = '#20364a', '#246b96', '#c45d38', '#a5b0bb'
    identity, cysteine, cdr = '#d8efe4', '#ffe1a0', '#e8b35e'
    fig = plt.figure(figsize=(15, 11.4), facecolor='white')
    grid = fig.add_gridspec(4, 1, height_ratios=[0.75, 1.65, 1.1, 2.8],
                           top=.88, bottom=.09, left=.07, right=.97, hspace=.75)
    fig.text(.07, .96, 'An IgV-like protein in an iridovirus-like assembly', fontsize=23, weight='bold', color=ink)
    fig.text(.07, .93, 'Harvard Forest soil  ·  Public MetaVR assembly IMGVR_UViG_3300079393_000044', fontsize=11, color=ink)

    def genes(ax, low, high, height):
        for gene in data['neighbors']:
            lo, hi = map(int, gene['protein'].rsplit('_', 2)[-2:])
            if hi < low or lo > high:
                continue
            strand = json.loads(gene['locations'])[0]['strand']
            color = gray
            if 'virus' in gene['description'].lower() and gene['evalue'] and float(gene['evalue']) <= 1e-8:
                color = viral
            if gene['protein'] == data['candidate_id']:
                color = candidate
            start = lo if strand == 1 else hi
            ax.add_patch(FancyArrow(start, 0, strand*(hi-lo+1), 0, width=height,
                                   head_width=height, head_length=min((high-low)*.012, (hi-lo+1)*.4),
                                   length_includes_head=True, facecolor=color, edgecolor='white', lw=.65))
        ax.set_xlim(low, high)
        ax.set_yticks([])
        for name in ['left','right','top']:
            ax.spines[name].set_visible(False)
        ax.spines['bottom'].set_color('#bcc6ce')
        ax.tick_params(axis='x', color='#bcc6ce', labelcolor=ink, labelsize=9)

    a = fig.add_subplot(grid[0])
    genes(a, 0, 17620, .35)
    a.add_patch(Rectangle((14000,-.34),3620,.68,facecolor='none',edgecolor=candidate,lw=1.2,linestyle='--'))
    a.set_ylim(-.52,.52)
    a.set_xticks([0,5000,10000,15000,17620], ['0','5','10','15','17.62 kb'])
    a.set_title('A   Entire assembly: 17,620 bp and 27 deposited proteins', loc='left', weight='bold', color=ink, pad=12)
    fig.legend(handles=[Patch(color=viral,label='Reviewed viral homolog (E ≤ 10⁻⁸)'),
                        Patch(color=candidate,label='IgV-like candidate'),Patch(color=gray,label='Other / unassigned')],
               loc='upper left',bbox_to_anchor=(.065,.79),ncol=3,frameon=False,fontsize=9,borderaxespad=0)

    b = fig.add_subplot(grid[1])
    genes(b, 14000, 17620, .25)
    b.set_ylim(-1.0,1.15)
    b.set_xticks([14000,15000,16000,17000,17620], ['14,000','15,000','16,000','17,000','17,620 bp'])
    b.set_title('B   Immediate genomic neighborhood — arrows are deposited CDSs, including stop codons', loc='left', weight='bold', color=ink, pad=13)
    for x,y,label,color in [
        (14375,.35,'TFIIS-like\nQ197A5 · IIV3\nE = 3.73 × 10⁻²⁰',viral),
        (15145,.36,'IgV-like candidate\n360 aa\n14,604–15,686',candidate),
        (15953,-.35,'045L-like\nQ6GZT1 · frog virus 3\nE = 8.72 × 10⁻⁹',viral),
        (16538,.43,'Unassigned',gray),
        (17375,.35,'Bro-N-like\nQ91FN5 · IIV6\nE = 2.89 × 10⁻⁷⁰',viral)]:
        b.text(x,y,label,ha='center',va='bottom' if y>0 else 'top',fontsize=10,color=color,
               weight='bold' if color==candidate else 'normal')

    c = fig.add_subplot(grid[2])
    c.set_xlim(1,360); c.set_ylim(-.6,1.0)
    c.add_patch(Rectangle((1,.18),359,.28,facecolor='#e1e6eb',edgecolor='none'))
    c.add_patch(Rectangle((25,.18),85,.28,facecolor=candidate,edgecolor='none'))
    c.text(67,.65,'IGK_V profile: aa 25–109',ha='center',color=candidate,fontsize=10,weight='bold')
    c.text(234,.32,'Remaining sequence: ancestry / function unassigned',ha='center',va='center',fontsize=10,color=ink)
    c.plot([21,107],[-.12,-.12],color=viral,lw=3)
    c.plot([21,21],[-.2,-.04],color=viral,lw=1)
    c.plot([107,107],[-.2,-.04],color=viral,lw=1)
    c.text(117,-.13,'BLASTP aligns aa 21–107: 87 of 360 residues (24.2%)',va='center',color=viral,fontsize=10)
    c.set_xticks([1,100,200,300,360], ['1','100','200','300','360 aa']); c.set_yticks([])
    for side in ['left','right','top']: c.spines[side].set_visible(False)
    c.spines['bottom'].set_color('#bcc6ce')
    c.set_title('C   The IgV match covers only the N-terminal portion of this protein', loc='left', weight='bold',color=ink,pad=14)

    d = fig.add_subplot(grid[3]); d.axis('off'); d.set_xlim(-9.2,48); d.set_ylim(-.65,4.75)
    d.set_title('D   Local alignment to mouse κ variable region T1 (P01637)',loc='left',weight='bold',color=ink,pad=32)
    d.text(-9.2,4.9,'33 / 89 identical positions (37.1%)  ·  55.1 bits  ·  E = 3.18 × 10⁻⁸ against reviewed UniProt',fontsize=10,color=ink)
    columns = data['alignment_columns']
    for chunk, start in enumerate(range(0,len(columns),45)):
        rows=columns[start:start+45]; y=3.75-chunk*2.3
        for index,row in enumerate(rows):
            exact=row['identity']; conserved_cys=exact and row['query_residue']=='C'
            if exact:
                d.add_patch(Rectangle((index,y-.76),1,.99,facecolor=cysteine if conserved_cys else identity,edgecolor='none'))
            for residue,dy in [(row['query_residue'],0),(row['reference_residue'],-.55)]:
                d.text(index+.5,y+dy,residue,fontfamily='DejaVu Sans Mono',fontsize=12.5,
                       ha='center',va='center',color=ink,weight='bold' if exact else 'normal')
        # Labels are the mouse reference's annotation, not predicted viral CDRs.
        regions=[]
        for i,row in enumerate(rows):
            name=row['reference_region']
            if not name: continue
            if regions and regions[-1][2]==name:
                regions[-1][1]=i+1
            else: regions.append([i,i+1,name])
        for lo,hi,name in regions:
            is_cdr=name.startswith('Complementarity')
            d.add_patch(Rectangle((lo,y+.38),hi-lo,.17,facecolor=cdr if is_cdr else '#c5ced7',edgecolor='white',lw=.4))
            if hi-lo>=3:
                label=('CDR' if is_cdr else 'FR')+name.rsplit('-',1)[-1]
                d.text((lo+hi)/2,y+.68,label,ha='center',va='center',fontsize=9,color=ink)
        for field,label,dy in [('query_position','Candidate',0),('reference_position','Mouse T1',-.55)]:
            coords=[r[field] for r in rows if r[field] is not None]
            d.text(-.65,y+dy,f'{label} {coords[0]:>3}',ha='right',va='center',fontsize=10,fontfamily='DejaVu Sans Mono',color=ink)
            d.text(len(rows)+.65,y+dy,str(coords[-1]),ha='left',va='center',fontsize=10,fontfamily='DejaVu Sans Mono',color=ink)
    d.text(-9.2,-.42,'Green: identity. Gold: conserved cysteines. Bars: mouse framework (FR) and CDR annotations. Gaps: “–”.',fontsize=10,color=ink)
    fig.text(.07,.035,'T1 CDR3 (aa 109–117) is outside the alignment. Viral context is supported; antibody ancestry and a donor allele are unproven.',
             fontsize=11,color=ink,weight='bold')
    products=[]
    for suffix in ['png','svg','pdf']:
        path=figures/f'igv_candidate_context.{suffix}'
        metadata={'Date':None} if suffix=='svg' else ({'CreationDate':None,'ModDate':None} if suffix=='pdf' else {})
        fig.savefig(path,dpi=180,metadata=metadata)
        products.append(path)
    plt.close(fig)
    report={'status':'complete','completed_utc':utc(),'command':[sys.executable,*sys.argv],
            'python':sys.version,'matplotlib':matplotlib.__version__,'script_sha256':sha256(Path(__file__)),
            'input_sha256':{str(source):sha256(source)},
            'outputs':{str(p.relative_to(ROOT)):sha256(p) for p in products},
            'interpretation':'Measured alignment and deposited assembly context; no allele, function, or capture inferred by the figure'}
    (output/'figure_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Generated',*[str(p.relative_to(ROOT)) for p in products])


if __name__ == '__main__':
    main()
