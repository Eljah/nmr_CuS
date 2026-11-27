import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from powder_fit import ROOT,NQR_SETS
from validate_and_plot import exact_result

joint=json.loads((ROOT/'joint_results.json').read_text())
best=min(joint,key=lambda r:r['SSE_joint'])
alternative=[r for r in joint if r['seed']==67][0]
fig,axs=plt.subplots(2,1,figsize=(9,6.4),sharex=True)
converged=[]
for record,color,label in [(best,'#087b94','Joint fit: eta1 = 0.150'),(alternative,'#cb733b','Alternative: eta1 = 0.925')]:
    sse=0
    for ax,which in zip(axs,['long','short']):
        rec=record[which]
        fit,pred,meta=exact_result(rec,160,320)
        sse+=meta['normalized_SSE']
        out={'seed':record['seed'],'series_ms':rec['series_ms'],**meta}
        converged.append(out)
        name=f'joint_converged_seed{record["seed"]}_{rec["series_ms"]}ms.csv'
        np.savetxt(ROOT/name,np.column_stack([fit.x,fit.y,pred,fit.mask]),delimiter=',',header='frequency_MHz,observed,predicted,used_in_fit',comments='')
        ax.plot(fit.x,pred,color=color,lw=1.4,label=label)
        if record is best:
            ax.plot(fit.x,np.where(fit.mask,fit.y,np.nan),'k.',ms=3,label='Observed outside coil mask')
        ax.set_title(f'Repetition time {rec["series_ms"]} ms',loc='left',fontsize=11)
    print({'seed':record['seed'],'converged_joint_SSE':sse},flush=True)
for ax in axs:
    for f in [75.941,81.346]:ax.axvspan(f-.12,f+.12,color='0.91',zorder=-1)
    ax.set_ylabel('Stored echo integral (a.u.)');ax.grid(alpha=.2);ax.set_ylim(bottom=0)
axs[0].legend(fontsize=8,loc='upper right');axs[-1].set_xlabel('Frequency (MHz)');axs[-1].set_xlim(72,84.8)
fig.tight_layout();fig.savefig(ROOT/'joint_powder_fit.png',dpi=240,bbox_inches='tight');fig.savefig(ROOT/'joint_powder_fit.pdf',bbox_inches='tight')
(ROOT/'joint_numerical_validation.json').write_text(json.dumps(converged,indent=2))

rows=[]
for name in ['thesis2009','publication2013']:
    rows+=json.loads((ROOT/f'summary_{name}.json').read_text())
for record in [best,alternative]:
    rows += [record['long'],record['short']]
with (ROOT/'model_comparison.csv').open('w',newline='') as f:
    fields=['NQR_set','model','series_ms','coil_mask_halfwidth_MHz','eta_s1','eta_s2','eta_s3','eta_s4','nuQ_s1_MHz','nuQ_s2_MHz','nuQ_s3_MHz','nuQ_s4_MHz','nuZ_s1_MHz','nuZ_s2_MHz','nuZ_s3_MHz','nuZ_s4_MHz','sigma_MHz','rho65','normalized_RMS','normalized_SSE','boundary_note']
    w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
    for r in rows:
        notes=[]
        if abs(r['sigma_MHz']-.25)<1e-5:notes.append('sigma upper bound')
        if abs(r['sigma_MHz']-.025)<1e-5:notes.append('sigma lower bound')
        if any(min(abs(v-75.2),abs(v-76.3))<1e-5 for v in r['effective_63_zeeman_MHz']):notes.append('frequency bound in separate-series exploratory model')
        if any(min(abs(v),abs(v-1))<1e-5 for v in r['eta']):notes.append('eta boundary')
        row={'NQR_set':r['dataset'],'model':r['tag'],'series_ms':r['series_ms'],'coil_mask_halfwidth_MHz':r['mask_half_width_MHz'],'sigma_MHz':r['sigma_MHz'],'rho65':r['isotope65_scale'],'normalized_RMS':r['normalized_RMS'],'normalized_SSE':r['normalized_SSE'],'boundary_note':'; '.join(notes)}
        for k in range(4):
            row[f'eta_s{k+1}']=r['eta'][k];row[f'nuQ_s{k+1}_MHz']=r['nuQ63_MHz'][k];row[f'nuZ_s{k+1}_MHz']=r['effective_63_zeeman_MHz'][k]
        w.writerow(row)

# A compact machine-readable table with the two competing joint solutions.
tab=[]
for record in [best,alternative]:
    r=record['long']
    for i,pair in enumerate(NQR_SETS['thesis2009']):
        tab.append({'solution_seed':record['seed'],'group':f's{i+1}','nuNQR63_MHz':pair[0],'nuNQR65_MHz':pair[1],'eta':r['eta'][i],'nuQ63_MHz':r['nuQ63_MHz'][i],'effective63_zeeman_MHz':r['effective_63_zeeman_MHz'][i],'amplitude500':r['amplitudes'][i],'amplitude20':record['short']['amplitudes'][i]})
with (ROOT/'joint_candidate_parameters.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(tab[0]));w.writeheader();w.writerows(tab)
