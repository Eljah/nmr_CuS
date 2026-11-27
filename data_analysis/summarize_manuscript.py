"""Reproduce descriptive numbers used in the anilite manuscript.

Run analyze_measurements.py first. Uses only exported stored measurements;
does not fit spectra, infer T1/T2, subtract baseline, or normalize amplitudes.
Literature NQR frequencies are explicitly transcribed constants, not new data.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import numpy as np

BASE = Path(__file__).resolve().parent
OUT = BASE / 'export'
WINDOWS = [
    ('broad_lower', 75.20, 75.70, 11),
    ('broad_upper', 76.10, 76.40, 7),
    ('narrow_lower', 75.85, 75.95, 3),
    ('narrow_upper', 81.25, 81.35, 3),
]
# Gainov et al., European Journal of Mineralogy 25 (2013), p. 572.
# DOI 10.1127/0935-1221/2013/0025-2325. Labels A-D are descriptive.
NQR_PAIRS = [('A',11.90,11.02),('B',15.21,14.10),('C',16.26,15.05),('D',16.81,15.55)]

def write_csv(name, fields, rows):
    with (OUT/name).open('w', newline='', encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields)
        w.writeheader(); w.writerows(rows)

def read_csv(name):
    p=OUT/name
    if not p.exists():
        raise SystemExit(f'Missing {name}. First run: python data_analysis/analyze_measurements.py')
    return np.genfromtxt(p, delimiter=',',names=True),hashlib.sha256(p.read_bytes()).hexdigest()

def window_summary():
    a,digest=read_csv('opj_Spec05122016.csv')
    # Frequency values encode a nominal 0.05 MHz grid. Convert to integer kHz
    # only for interval membership to avoid a floating-point boundary exclusion.
    kHz=np.rint(a['Freq']*1000).astype(np.int64)
    assert np.all(np.abs(a['Freq']*1000-kHz)<1e-6)
    assert np.all(kHz%50==0)
    valid=np.ones(len(a),dtype=bool)
    conflicts=[]
    for f in np.unique(kHz):
        ix=np.flatnonzero(kHz==f)
        if len(ix)>1:
            valid[ix]=False
            conflicts.extend(int(a['origin_row'][i]) for i in ix)
    assert conflicts==[188,189],conflicts
    rows=[]; points=[]
    for name,lo,hi,count in WINDOWS:
        select=valid&(kHz>=round(lo*1000))&(kHz<=round(hi*1000))
        s=a[select]
        assert len(s)==count,(name,len(s),count)
        assert np.all(s['IntQ']>0) and np.all(s['IntQ1']>0)
        short=math.fsum(float(x) for x in s['IntQ'])
        long=math.fsum(float(x) for x in s['IntQ1'])
        ratio=long/short
        row={'window':name,'lower_MHz_inclusive':lo,'upper_MHz_inclusive':hi,
             'number_of_points':len(s),'sum_I20':short,'sum_I500':long,
             'ratio_sum_I500_over_sum_I20':ratio}
        rows.append(row)
        for r in s:
            points.append({'window':name,'origin_row':int(r['origin_row']),
                'frequency_MHz':float(r['Freq']),'File_20ms':int(r['File']),
                'File_500ms':int(r['File1']),'I20':float(r['IntQ']),
                'I500':float(r['IntQ1']),'pointwise_I500_over_I20':float(r['IntQ1']/r['IntQ'])})
        # Independent summation check; this does not define a second statistic.
        assert math.isclose(ratio,float(np.sum(s['IntQ1'])/np.sum(s['IntQ'])),rel_tol=1e-14)
    write_csv('manuscript_window_ratios.csv',list(rows[0]),rows)
    write_csv('manuscript_window_points.csv',list(points[0]),points)
    return {'source':'opj_Spec05122016.csv','source_sha256':digest,
        'definition':'sum(I500 over included grid points) / sum(I20 over the same points)',
        'bounds':'inclusive, membership checked on nominal integer-kHz grid',
        'excluded_conflict_origin_rows':conflicts,'baseline_subtracted':False,'amplitude_normalization':False,
        'windows':rows}

def peak_summary():
    rows=[];hashes={}
    for name,reptime in [('spec20ms',20),('spec500ms',500)]:
        a,digest=read_csv(f'opj_{name}.csv');hashes[f'opj_{name}.csv']=digest
        peaks=[]
        for lo,hi in [(74.0,78.0),(80.0,83.0)]:
            ix=np.flatnonzero((a['FreqMHz']>=lo)&(a['FreqMHz']<=hi))
            m=ix[np.argmax(a['Real'][ix])]
            peaks.append((int(a['origin_row'][m]),float(a['FreqMHz'][m]),float(a['Real'][m])))
        row={'repetition_ms':reptime,'source':f'opj_{name}.csv',
             'lower_group_search_from_MHz':74.0,'lower_group_search_to_MHz':78.0,
             'higher_group_search_from_MHz':80.0,'higher_group_search_to_MHz':83.0,
             'lower_peak_origin_row':peaks[0][0],'lower_peak_MHz':peaks[0][1],
             'lower_peak_real_amplitude':peaks[0][2],
             'higher_peak_origin_row':peaks[1][0],'higher_peak_MHz':peaks[1][1],
             'higher_peak_real_amplitude':peaks[1][2],
             'higher_over_lower_frequency_ratio':peaks[1][1]/peaks[0][1]}
        rows.append(row)
    write_csv('manuscript_peak_summary.csv',list(rows[0]),rows)
    return {'method':'maximum stored Real value in each stated search interval; no fit',
        'source_sha256':hashes,'results':rows,
        'warning':'Grid coordinates are not peak-fit uncertainty estimates; isotope consistency does not identify the material producing a peak.'}

def nqr_bounds():
    rows=[]
    lower_factor=math.sqrt(3)/2
    for pair,f63,f65 in NQR_PAIRS:
        lo=f63*lower_factor;hi=f63
        assert math.isclose(lo*math.sqrt(1+1/3),f63,rel_tol=1e-14)
        assert math.isclose(hi*math.sqrt(1+0/3),f63,rel_tol=1e-14)
        rows.append({'pair_label':pair,'nu_NQR_63_MHz':f63,'nu_NQR_65_MHz':f65,
            'temperature_K':4.2,'abs_nuQ_63_min_MHz_eta_1':lo,'abs_nuQ_63_max_MHz_eta_0':hi})
    write_csv('manuscript_nqr_bounds.csv',list(rows[0]),rows)
    return {'source':'Gainov et al. (2013), European Journal of Mineralogy 25, 569-578, p. 572',
        'source_doi':'10.1127/0935-1221/2013/0025-2325','spin_I':1.5,
        'nuQ_definition':'3 e Q Vzz / [2 I (2 I - 1) h]',
        'relation':'nu_NQR = abs(nuQ) * sqrt(1 + eta^2/3)',
        'eta_range_inclusive':[0,1],'lower_factor':lower_factor,'results':rows,
        'warning':'Bounds are due to unknown EFG asymmetry, not confidence intervals. No temperature correction from 4.2 K to 5 K has been applied.'}

def main():
    OUT.mkdir(exist_ok=True)
    result={'window_ratios':window_summary(),'processed_peak_coordinates':peak_summary(),'nqr_coupling_bounds':nqr_bounds()}
    (OUT/'manuscript_numeric_summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
