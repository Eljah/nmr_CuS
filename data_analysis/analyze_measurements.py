"""Read-only extraction and auditable plotting of supplied Anilite measurement files.

The OPJ data-list parser follows the public Origin 7.0552 format description:
https://github.com/jgonera/openopj/blob/master/docs/opj_format.markdown
It does not execute formulas, macros, or embedded files. It parses only the data
list and records graph labels as unassociated strings in a separate file.
"""
from pathlib import Path
import csv, hashlib, json, re, struct, zipfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
INPUT = BASE.parent / 'upload'
OUT = BASE / 'export'
OUT.mkdir(exist_ok=True)

def extract_sources():
    """Materialize only the three user-supplied measurement archives."""
    manifest=[]
    for archive,dirname in [
        ('Anilite14072016.zi','anilite_raw'),('CuS.zi','cus_raw'),
        ('Attachments_irek.mukhamedshin@gmail.com_2016-07-04_09-39-32.zip','email_attachments')]:
        source=INPUT/archive;destination=(BASE/dirname).resolve()
        destination.mkdir(exist_ok=True)
        with zipfile.ZipFile(source) as z:
            infos=z.infolist()
            assert sum(i.file_size for i in infos)<100_000_000
            for i in infos:
                target=(destination/i.filename).resolve()
                assert target.is_relative_to(destination),i.filename
                assert not (i.external_attr>>16)&0o170000==0o120000,'Symlink member'
            z.extractall(destination)
            manifest.append({'archive':archive,'sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
                'destination':dirname,'members':[{'filename':i.filename,'size':i.file_size,
                    'archive_timestamp':list(i.date_time)} for i in infos]})
    (OUT/'source_archive_manifest.json').write_text(json.dumps(manifest,indent=2))

def write_csv(path, headers, rows):
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f); w.writerow(headers); w.writerows(rows)

def read_opj(path):
    data = path.read_bytes()
    signature, _ = data.split(b'\n', 1)
    assert signature == b'CPYA 4.2673 552#', signature
    cursor = len(signature) + 1
    def block():
        nonlocal cursor
        start = cursor
        n = struct.unpack_from('<I', data, cursor)[0]; cursor += 4
        assert data[cursor] == 10, (cursor, n)
        cursor += 1
        payload = data[cursor:cursor+n]; cursor += n
        assert len(payload) == n
        if n:
            assert data[cursor] == 10, (cursor, n)
            cursor += 1
        return start, payload
    _, global_header = block()
    assert len(global_header) == 39
    _, end = block(); assert not end
    result = {}; metadata = []
    while True:
        start, h = block()
        if not h: break
        assert len(h) == 123, (start, len(h))
        name = h[0x58:0x71].split(b'\0', 1)[0].decode('ascii')
        value_size = h[0x3d]
        data_type = struct.unpack_from('<H', h, 0x16)[0]
        total, first, last = struct.unpack_from('<III', h, 0x19)
        content_start, d = block()
        _, mask = block()
        assert not mask, 'Unexpected mask: needs additional parser'
        assert value_size == 10 and data_type & 0x100, (name, value_size, data_type)
        assert len(d) == total*value_size
        flags = np.frombuffer(d, dtype=np.dtype({'names':['flag','value'],
                          'formats':['<u2','<f8'],'offsets':[0,2],'itemsize':10}))
        assert np.all(flags['flag'] == 0), 'Text cell: needs explicit export'
        arr = flags['value'].copy()
        sentinel = arr == -1.23456789e-300
        arr[sentinel] = np.nan
        assert np.all(np.isfinite(arr[~sentinel]))
        result[name] = arr
        metadata.append({'name':name,'header_offset':start,'content_offset':content_start,
            'dtype_hex':hex(data_type),'value_size':value_size,'allocated_rows':total,
            'header_first_row':first,'header_last_row':last,'finite_rows':int(np.isfinite(arr).sum()),
            'missing_rows':int(sentinel.sum()),'first_values':arr[np.isfinite(arr)][:3].tolist(),
            'last_values':arr[np.isfinite(arr)][-3:].tolist()})
    return result, metadata, cursor

def save_opj_tables(columns):
    sheets = {}
    for key,val in columns.items():
        sheet, col = key.rsplit('_',1)
        sheets.setdefault(sheet,{})[col] = val
    for sheet,cols in sheets.items():
        names = list(cols)
        n = max(len(a) for a in cols.values())
        a = np.full((n,len(cols)),np.nan)
        for i,k in enumerate(names): a[:len(cols[k]),i] = cols[k]
        keep = np.any(np.isfinite(a),axis=1)
        rows = [[i+1]+[format(v,'.17g') if np.isfinite(v) else '' for v in row]
                for i,row in enumerate(a) if keep[i]]
        write_csv(OUT/f'opj_{sheet}.csv',['origin_row']+names,rows)
    return sheets

def graph_metadata(path, dataset_names):
    """Parse target window's first layer, annotations and dataset associations.

    Offsets for annotations and graph curves are documented by liborigin:
    https://sources.debian.org/src/labplot/2.5.0-2/liborigin/OriginAnyParser.cpp
    This intentionally does not reconstruct graph drawing or all axis options.
    """
    b=path.read_bytes();result={}
    def block(p):
        n=struct.unpack_from('<I',b,p)[0]
        assert b[p+4]==10
        d=b[p+5:p+5+n]
        assert len(d)==n
        if n:assert b[p+5+n]==10
        return p+5+n+(1 if n else 0),d
    for name in ['GraphT2','GrSpectra','T2a','Spec05122016']:
        candidates=[]
        for m in re.finditer(b'\x00\x00'+re.escape(name.encode())+b'\x00',b):
            p=m.start()-5
            if p<0:continue
            n=struct.unpack_from('<I',b,p)[0]
            if 195<=n<2048 and b[p+4]==10 and b[p+5+n]==10:candidates.append(p)
        assert len(candidates)==1,(name,candidates)
        p,h=block(candidates[0]);p,layer=block(p)
        entry={'window_header_offset':candidates[0],'annotations':[],'curves':[]}
        if name.startswith('Gr'):
            entry['x_axis_limits_raw']=list(struct.unpack_from('<ddd',layer,0x0f))
        while True:
            ah_offset=p;p,ah=block(p)
            if not ah:break
            p,d1=block(p);p,d2=block(p);p,d3=block(p)
            sec=ah[0x46:0x46+41].split(b'\0',1)[0].decode('cp1251')
            def txt(d):
                z=d.split(b'\0',1)[0]
                return z.decode('ascii') if all(32<=c<127 or c in [9,10,13] for c in z) else None
            entry['annotations'].append({'name':sec,'offset':ah_offset,'data1_text':txt(d1),'data2_text':txt(d2),'data3_text':txt(d3)})
        while True:
            curve_offset=p;p,c=block(p)
            if not c:break
            p,d=block(p)
            rec={'offset':curve_offset,'data_text':d.decode('cp1251',errors='replace')}
            if name.startswith('Gr'):
                yi=struct.unpack_from('<H',c,4)[0];xi=struct.unpack_from('<H',c,0x23)[0]
                rec.update({'x_dataset_index_1_based':xi,'y_dataset_index_1_based':yi,
                    'x_dataset':dataset_names[xi-1] if xi else None,'y_dataset':dataset_names[yi-1] if yi else None})
            entry['curves'].append(rec)
        result[name]=entry
    return result

def july_export():
    source = BASE/'anilite_raw/14072016_CuS_45_67.dat'
    a = np.genfromtxt(source,names=True)
    b = np.genfromtxt(BASE/'anilite_raw/14072016_CuS_45_58.dat',names=True)
    assert np.array_equal(a[:len(b)],b)
    rows=[]; masked=[]
    for i in range(len(a)):
        for j in [1,2,3]:
            file_no=int(a[f'File{j}'][i]); param=int(a[f'Satur{j}'][i]); q=float(a[f'IntQ{j}'][i])
            carry = j==1 and a['Freq'][i]>58
            if carry: assert file_no==418 and q==a['IntQ1'][129]
            status='EXCLUDE_CARRIED_FORWARD' if carry else ('PRESENT_SCALE_REGIME_UNVERIFIED' if file_no>=501 else 'PRESENT_UNIQUE_FILE')
            row=[float(a['Freq'][i]),j,param,file_no,q,status]
            rows.append(row)
            if carry: masked.append(row)
    headers=['frequency_nominal_MHz','series','Satur_raw_units_unconfirmed','raw_file_id','IntQ_as_supplied','audit_status']
    write_csv(OUT/'july_integrals_audited.csv',headers,rows)
    write_csv(OUT/'july_excluded_carried_forward.csv',headers,masked)
    src=BASE/'email_attachments/SpectrumIntegrals.dat'
    small=np.genfromtxt(src,names=True)
    write_csv(OUT/'july03_cus_integrals.csv',list(small.dtype.names),[list(r) for r in small])
    return a,small

def raw_trace_audit():
    rows=[]
    for p in sorted((BASE/'anilite_raw').glob('*_KN_ADC.dat')):
        a=np.loadtxt(p,skiprows=1)
        assert a.shape==(1024,4),p.name
        assert np.allclose(np.diff(a[:,0]),0.16,rtol=0,atol=1e-10),p.name
        tail=a[a[:,0]>125]
        rows.append([p.name,len(a),float(np.median(np.diff(a[:,0]))),float(np.std(tail[:,1])),
            float(np.std(tail[:,2])),float(np.mean(tail[:,3])),float(np.max(a[:,3])),
            float(np.max(np.abs(np.hypot(a[:,1],a[:,2])-a[:,3]))),hashlib.sha256(p.read_bytes()).hexdigest()])
    assert len(rows)==600
    assert len({row[-1] for row in rows})==600,'Identical raw files detected'
    write_csv(OUT/'july_raw_trace_audit.csv',
        ['filename','sample_count','dt_us','tail_Real_sd','tail_Imag_sd','tail_Module_mean','max_Module',
         'max_Module_consistency_error','sha256'],rows)

def opj_sweep_audit(sheets):
    s=sheets['Spec05122016']; f=s['Freq']
    valid=np.isfinite(f); status=np.full(len(f),'PRESENT_UNIQUE',dtype=object)
    status[~valid]='ORIGIN_MISSING_SENTINEL'
    for frequency in np.unique(f[valid]):
        indexes=np.flatnonzero(f==frequency)
        if len(indexes)>1:valid[indexes]=False;status[indexes]='EXCLUDE_DUPLICATE_FREQUENCY_AND_FILE_CONFLICT'
    rows=[]
    for i in range(len(f)):
        if not np.isfinite(f[i]):continue
        rows.append([i+1,f[i],s['File'][i],s['IntQ'][i],s['File1'][i],s['IntQ1'][i],status[i]])
    write_csv(OUT/'opj_Spec05122016_audit.csv',
              ['origin_row','Freq_MHz','File_20ms','IntQ_20ms','File_500ms','IntQ_500ms','audit_status'],rows)
    return valid

def format_axes(ax):
    ax.spines[['top','right']].set_visible(False)
    ax.grid(axis='y',color='0.9',linewidth=.6)
    ax.tick_params(direction='out')

def savefig(fig,name):
    fig.savefig(OUT/f'{name}.png',dpi=220,bbox_inches='tight')
    fig.savefig(OUT/f'{name}.pdf',bbox_inches='tight')
    plt.close(fig)

def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':11,
                         'axes.titlesize':12,'legend.fontsize':9,'pdf.fonttype':42,'ps.fonttype':42})
    extract_sources()
    raw_trace_audit()
    cols,meta,end=read_opj(INPUT/'Anilite.OPJ')
    sheets=save_opj_tables(cols)
    valid_sweep=opj_sweep_audit(sheets)
    graphmeta=graph_metadata(INPUT/'Anilite.OPJ',list(cols))
    (OUT/'opj_graph_metadata.json').write_text(json.dumps(graphmeta,indent=2,ensure_ascii=False))
    (OUT/'opj_inventory.json').write_text(json.dumps({'source_sha256':hashlib.sha256((INPUT/'Anilite.OPJ').read_bytes()).hexdigest(),
        'dataset_end_offset':end,'parser_reference':'https://github.com/jgonera/openopj/blob/master/docs/opj_format.markdown',
        'columns':meta},indent=2))
    a,small=july_export()
    fig,ax=plt.subplots(figsize=(9.1,4.3))
    for j,c in [(1,'#777777'),(2,'#086f98'),(3,'#bf552c')]:
        valid=np.ones(len(a),dtype=bool) if j!=1 else a['Freq']<=58
        ax.plot(a['Freq'][valid],a[f'IntQ{j}'][valid],'.-',ms=3,lw=.8,color=c,label=f'Satur = {int(a[f"Satur{j}"][0]):,} (stored value)')
    ax.axvline(58,color='0.5',ls=':',lw=.8)
    ax.axvline(62.05,color='0.5',ls='--',lw=.8)
    ax.axvspan(62.05,67.1,color='0.96',zorder=-1)
    ax.set(xlabel='Nominal frequency (MHz)',ylabel='Stored echo integral IntQ (a.u.)',xlim=(45,67.1))
    ax.set_title('July 14 frequency sweep from Anilite14072016.zi',loc='left')
    ax.text(.99,.94,'Series 1: data stop at 58.0 MHz\n90 carried-forward entries excluded',transform=ax.transAxes,ha='right',va='top',fontsize=9,
           bbox={'facecolor':'white','edgecolor':'0.8','boxstyle':'round,pad=0.4'})
    ax.text(.865,.26,'Amplitude regime\nunverified\nabove 62.0 MHz',transform=ax.transAxes,ha='center',va='center',fontsize=9,color='0.35')
    ax.legend(loc='upper left',frameon=False);format_axes(ax);fig.tight_layout()
    savefig(fig,'figure_july_sweep_audited')
    sp=sheets['Spec05122016']
    f=sp['Freq'];n=len(f)
    fig,axs=plt.subplots(2,1,figsize=(9.1,7.0),sharex=True)
    for key,label,c in [('IntQ','20 ms','#bf552c'),('IntQ1','500 ms','#086f98')]:
        y=sp[key][:n].copy();y[~valid_sweep]=np.nan
        axs[0].plot(f,y,'.-',ms=3,lw=.9,color=c,label=label)
    axs[0].set(ylabel='Stored integral IntQ (a.u.)')
    axs[0].set_title('(a) Frequency-stepped integrated echo',loc='left')
    for sheet,label,c in [('spec20ms','20 ms','#bf552c'),('spec500ms','500 ms','#086f98')]:
        s=sheets[sheet]
        axs[1].plot(s['FreqMHz'],s['Real'],color=c,lw=.6,label=label)
    axs[1].set(xlabel='Frequency (MHz)',ylabel='Stored real component (a.u.)',xlim=(71.7,85.1))
    axs[1].set_title('(b) Processed spectra retained in the Origin project',loc='left')
    for ax in axs:format_axes(ax);ax.legend(title='Pulse-sequence repetition',frameon=False)
    fig.suptitle('Anilite.OPJ: spectral data and repetition-time dependence',x=.09,ha='left',fontsize=13)
    fig.tight_layout();savefig(fig,'figure_opj_spectra')
    t=sheets['T2a'];valid=np.isfinite(t['TwoTau'])&np.isfinite(t['IntQ'])
    fig,axs=plt.subplots(1,2,figsize=(9.1,3.8))
    for key,label,c in [('IntRe','Real','#086f98'),('IntIm','Imaginary','#777777'),('IntQ','IntQ','#bf552c')]:
        axs[0].plot(t['TwoTau'][valid],t[key][valid],'o',ms=3,color=c,label=label)
    axs[0].set(xlabel='TwoTau (stored units; conversion unverified)',ylabel='Stored integral (a.u.)')
    axs[0].legend(frameon=False)
    axs[0].set_title('(a) T2a worksheet values',loc='left')
    axs[1].plot(small['Freq'],small['Integral'],'o-',ms=3,lw=.9,color='#086f98')
    axs[1].set(xlabel='Frequency (MHz)',ylabel='Stored integral (a.u.)')
    axs[1].set_title('(b) Comparison spectrum',loc='left',fontsize=11)
    for ax in axs:format_axes(ax)
    fig.tight_layout();savefig(fig,'figure_t2_and_comparison')
    summary={'opj_sheets':{},'july_series':[],'cus03_peak':{},'t2':{}}
    summary['opj_sweep_audit']={'finite_rows':int(np.isfinite(f).sum()),
        'unique_frequencies':len(np.unique(f[np.isfinite(f)])),
        'excluded_origin_rows':[188,189], 'reason':'duplicate frequency 81.4 MHz and file IDs 425/426, conflicting integral values',
        'valid_rows_after_exclusion':int(valid_sweep.sum())}
    summary['july_amplitude_boundary']={'frequency_MHz':62.1,'first_file_new_regime':501,
        'interpretation':'Instrumental or processing origin cannot be resolved from supplied metadata; no renormalization applied.'}
    for name,s in sheets.items():
        summary['opj_sheets'][name]={k:{'finite':int(np.isfinite(v).sum()),'min':float(np.nanmin(v)),'max':float(np.nanmax(v))} for k,v in s.items()}
    for j in [1,2,3]:
        valid=a['Freq']<=58 if j==1 else np.ones(len(a),dtype=bool)
        ai=a[valid];k=int(np.argmax(ai[f'IntQ{j}']))
        summary['july_series'].append({'Satur_raw':int(ai[f'Satur{j}'][0]),'n_unique_traces':len(ai),
            'frequency_range':[float(ai['Freq'][0]),float(ai['Freq'][-1])],
            'max_grid_frequency_MHz':float(ai['Freq'][k]),'max_IntQ':float(ai[f'IntQ{j}'][k]),
            'file_at_max':int(ai[f'File{j}'][k])})
    k=np.argmax(small['Integral']);summary['cus03_peak']={'frequency_MHz':float(small['Freq'][k]),'IntQ':float(small['Integral'][k]),'File':int(small['File'][k])}
    for k in ['TwoTau','IntRe','IntIm','IntQ','File']:summary['t2'][k]=t[k][np.isfinite(t['TwoTau'])].tolist()
    (OUT/'data_summary.json').write_text(json.dumps(summary,indent=2))
    fig,axs=plt.subplots(2,1,figsize=(7.0,5.5),sharex=True)
    for key,label,c in [('IntQ','20 ms','#bf552c'),('IntQ1','500 ms','#086f98')]:
        y=sp[key][:n].copy();y[~valid_sweep]=np.nan
        axs[0].plot(f,y,'.-',ms=2.5,lw=.7,color=c,label=label)
    axs[0].set(ylabel='Echo integral (a.u.)')
    for sheet,label,c in [('spec20ms','20 ms','#bf552c'),('spec500ms','500 ms','#086f98')]:
        s=sheets[sheet]
        axs[1].plot(s['FreqMHz'],s['Real'],color=c,lw=.55,label=label)
    axs[1].set(xlabel='Frequency (MHz)',ylabel='Spectral amplitude (a.u.)',xlim=(71.7,85.1))
    for i,ax in enumerate(axs):
        format_axes(ax)
        ax.text(.02,.96,f'({chr(97+i)})',transform=ax.transAxes,ha='left',va='top',fontsize=11,fontweight='bold')
        ax.legend(title='Repetition time',frameon=False,loc='upper right')
        ax.ticklabel_format(axis='y',style='sci',scilimits=(3,3),useMathText=True)
    axs[1].text(.02,.82,'5 K\nNominal field: 6.720566 T',transform=axs[1].transAxes,fontsize=9,va='top')
    fig.tight_layout();savefig(fig,'figure_opj_manuscript')
    print(json.dumps({'output_directory':str(OUT),'opj_sheets':list(sheets),'columns':len(cols),
        'july_series':summary['july_series']},indent=2))

if __name__=='__main__':main()
