"""Inventory originals, extracted measurements and embedded Origin data.

Read-only with respect to every input. Writes provenance reports and exact
low-temperature NQR T2 table; performs no fit or image digitization.
Run extract_embedded_origin.py before this script.
"""
from pathlib import Path
from collections import Counter
import ast
import csv
import hashlib
import json
import zipfile

BASE=Path(__file__).resolve().parent.parent
OUT=BASE/'provenance'
OUT.mkdir(exist_ok=True)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def record(path,**extra):
    b=path.read_bytes()
    return {'path':str(path.relative_to(BASE)),'filename':path.name,
            'bytes':len(b),'sha256':digest(b),**extra}


ROLES={
 'Anilite.OPJ':'Origin project: finite-field NMR integrals, processed spectra and echo-decay table',
 'Anilite14072016.zi':'ZIP measurement archive: 600 I/Q time traces and two frequency-integral tables',
 'CuS.zi':'ZIP comparison-measurement archive: 51 I/Q time traces, files03--53',
 'Attachments_irek.mukhamedshin@gmail.com_2016-07-04_09-39-32.zip':'Email-attachment ZIP: SpectrumIntegrals.dat and Remarks.pdf',
 'Remarks.pdf':'Two-page experimental notes; duplicate also retained in email-attachment archive',
 'Ilya.pdf':'53-page methodological excerpt on NMR/NQR and spectral reconstruction; not a low-frequency anilite numeric dataset',
 'diploma_magister_geerite_2009_5.doc':'Original thesis, including five embedded Origin numeric projects',
 'diploma_in_brief_Evlampiev_3.rtf':'Short thesis summary; no embedded Origin objects found',
 'geerite to digenite.docx':'Working manuscript with four embedded Origin projects',
 'papers_20.zip':'14-source literature/document archive, including alternate thesis RTF and manuscript DOCX',
}

DERIVED={
 'Anilite.OPJ':['data_analysis/export/opj_inventory.json','data_analysis/export/opj_graph_metadata.json','data_analysis/export/opj_*.csv'],
 'Anilite14072016.zi':['data_analysis/anilite_raw','data_analysis/export/july_integrals_audited.csv','data_analysis/export/july_raw_trace_audit.csv'],
 'CuS.zi':['data_analysis/cus_raw'],
 'Attachments_irek.mukhamedshin@gmail.com_2016-07-04_09-39-32.zip':['data_analysis/email_attachments','data_analysis/export/july03_cus_integrals.csv'],
 'Remarks.pdf':['data_analysis/email_attachments/Remarks.pdf'],
 'Ilya.pdf':[],
 'diploma_magister_geerite_2009_5.doc':['new_evidence/thesis/diploma_magister_geerite_2009_5.pdf','new_evidence/thesis/diploma_magister_geerite_2009_5.txt','data_analysis/thesis_embedded/diploma5_01.opj','data_analysis/thesis_embedded/diploma5_02.opj','data_analysis/thesis_embedded/diploma5_03.opj','data_analysis/thesis_embedded/diploma5_04.opj','data_analysis/thesis_embedded/diploma5_05.opj'],
 'diploma_in_brief_Evlampiev_3.rtf':['new_evidence/thesis/diploma_in_brief_Evlampiev_3.pdf','new_evidence/thesis/diploma_in_brief_Evlampiev_3.txt'],
 'geerite to digenite.docx':['new_evidence/thesis/geerite to digenite.pdf','new_evidence/thesis/geerite_to_digenite_docx.txt','data_analysis/thesis_embedded/geerite_to_digenite_oleObject5_01.opj','data_analysis/thesis_embedded/geerite_to_digenite_oleObject6_01.opj','data_analysis/thesis_embedded/geerite_to_digenite_oleObject7_01.opj','data_analysis/thesis_embedded/geerite_to_digenite_oleObject8_01.opj'],
 'papers_20.zip':['new_evidence/papers/raw','new_evidence/papers/manifest.json','new_evidence/thesis/diploma_magister_geerite_2009_10.pdf','data_analysis/thesis_embedded/embedded_origin_manifest.json'],
}


def trace_shape_audit(folder,pattern):
    counts=Counter();headers=Counter();steps=Counter()
    paths=sorted((BASE/folder).glob(pattern))
    for p in paths:
        lines=p.read_text().splitlines();headers[lines[0]]+=1
        rows=[line.split() for line in lines[1:] if line.strip()]
        counts[(len(rows),tuple(sorted({len(row) for row in rows})))]+=1
        times=[float(row[0]) for row in rows]
        assert all(abs((b-a)-0.16)<1e-8 for a,b in zip(times,times[1:]))
        steps[round(times[1]-times[0],10)]+=1
    assert dict(counts)=={(1024,(4,)):len(paths)}
    return {'files_checked':len(paths),'all_have_1024_rows_and_4_columns':True,
            'headers':dict(headers),'numeric_time_step_counts':dict(steps),
            'time_unit_note':'Anilite headers explicitly say mks; CuS headers say Time without a printed unit.'}


def archive_audit():
    results=[]
    for source,destination in [
       ('Anilite14072016.zi','data_analysis/anilite_raw'),
       ('CuS.zi','data_analysis/cus_raw'),
       ('Attachments_irek.mukhamedshin@gmail.com_2016-07-04_09-39-32.zip','data_analysis/email_attachments'),
       ('papers_20.zip','new_evidence/papers/raw')]:
        p=BASE/'upload'/source;members=[]
        with zipfile.ZipFile(p) as z:
            for info in z.infolist():
                if info.is_dir():continue
                b=z.read(info);target=BASE/destination/info.filename
                identical=target.is_file() and target.read_bytes()==b
                members.append({'member':info.filename,'bytes':len(b),'sha256':digest(b),
                    'extracted_path':str(target.relative_to(BASE)),
                    'extracted_exists':target.is_file(),'extracted_byte_identical':identical})
        results.append({'source':str(p.relative_to(BASE)),'destination':destination,
            'member_count':len(members),'uncompressed_bytes':sum(m['bytes'] for m in members),
            'all_extractions_identical':all(m['extracted_byte_identical'] for m in members),'members':members})
    return results


def table_info(path,delimiter=','):
    with path.open(newline='',encoding='utf-8-sig') as f:
        rows=list(csv.DictReader(f,delimiter=delimiter))
    return rows


def embedded_audit():
    origin=json.loads((BASE/'data_analysis/thesis_embedded/embedded_origin_manifest.json').read_text())
    projects=origin['origin_projects']
    versions=[]
    for i in range(5):
        versions.append({'original':projects[i]['path'],'alternate':projects[i+5]['path'],
                        'project_bytes_identical':projects[i]['sha256']==projects[i+5]['sha256']})
    for i in range(4):
        a=BASE/'data_analysis/thesis_embedded/numeric'/Path(projects[i]['path']).stem/'column_provenance.json'
        b=BASE/'data_analysis/thesis_embedded/numeric'/Path(projects[i+10]['path']).stem/'column_provenance.json'
        ca=json.loads(a.read_text())['columns'];cb=json.loads(b.read_text())['columns']
        same=[(x['name'],x['values']) for x in ca]==[(x['name'],x['values']) for x in cb]
        versions.append({'original':projects[i]['path'],'alternate':projects[i+10]['path'],
            'project_bytes_identical':projects[i]['sha256']==projects[i+10]['sha256'],
            'numeric_columns_identical':same})
    assert all(v.get('numeric_columns_identical',v['project_bytes_identical']) for v in versions)
    p=BASE/'data_analysis/thesis_embedded/numeric/diploma5_02/column_provenance.json'
    columns={c['name']:c['values'] for c in json.loads(p.read_text())['columns']}
    scaling=[]
    for short in ('c','s1','s2','s3','s4'):
        pairs=[(a,b) for a,b in zip(columns['Data1_'+short],columns['Data1_'+short+'i65']) if a is not None and b is not None]
        scaling.append({'column_63':'Data1_'+short,'column_65':'Data1_'+short+'i65',
            'ratio_used_for_comparison':1.081,'paired_rows':len(pairs),
            'max_abs_difference_from_dividing_by_ratio':max(abs(b-a/1.081) for a,b in pairs)})
    return {'project_count':len(projects),'numeric_snapshot_count':5,
        'unique_thesis_numeric_columns':70,'project_equivalence':versions,
        'isotope_temperature_columns_scaling_check':scaling,
        'column_interpretation_note':'Data1_freq contains temperature in K (verified Graph5 x axis); i65 columns exactly match63/1.081 and are not asserted to be independent measurements.',
        'low_frequency_nqr_numeric_spectrum_found':False,
        'scope':'All extracted data-list columns in14 embedded Origin projects; no pixel digitization.'}


def low_temperature_t2():
    project='data_analysis/thesis_embedded/diploma5_05.opj'
    p=BASE/'data_analysis/thesis_embedded/numeric/diploma5_05/column_provenance.json'
    d=json.loads(p.read_text());cols={c['name']:c for c in d['columns']}
    graph=d['graph_first_layers'][0];assert graph['name']=='Graph6'
    assert graph['curves'][10]['x_dataset']=='T2s3_T'
    assert graph['curves'][10]['y_dataset']=='T2s3_B'
    assert graph['curves'][2]['x_dataset']=='T2datd_A'
    assert graph['curves'][2]['y_dataset']=='T2datd_C'
    output=[]
    for phase,group,x,y,err,legend in [
      ('covellite','c','T2datd_A','T2datd_C','T2datd_cerr',3),
      ('historically_geerite_reidentified_component','s3','T2s3_T','T2s3_B','T2s3_C',11)]:
        for row,(t,r,e) in enumerate(zip(cols[x]['values'],cols[y]['values'],cols[err]['values']),1):
            if t is None or r is None or t>10:continue
            output.append({'phase_label':phase,'spectral_group':group,'isotope':'63Cu',
                'temperature_K':t,'T2_inverse_ms_inv':r,'stored_errorbar_ms_inv':e,
                'T2_ms_reciprocal':1/r,'T2_us_reciprocal':1000/r,
                'source_project':project,'source_project_sha256':d['source_sha256'],
                'source_doc_stream':'ObjectPool/_1307219488/Contents','thesis_printed_page':29,'thesis_figure':11,
                'graph_name':'Graph6','legend_curve_index':legend,
                'temperature_column':x,'rate_column':y,'error_column':err,'origin_row':row})
    filename=OUT/'NQR_T2_LOW_T_EXACT.csv'
    with filename.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(output[0]));w.writeheader();w.writerows(output)
    summary={'rows':output,'rate_unit':'ms^-1, verified original Graph6 label and thesisFig11',
        'isotope_source':'Fig11 caption:63Cu','errorbar_convention':'Stored errorbar magnitudes; statistical interpretation not specified.',
        'group_mapping':'Spectroscopic s3, not an atomicCu3 assignment.',
        's1_s2_s4_lowest_saved_plotted_temperature_K':77,
        's3_temperature_interval_K':[6,10],
        's3_T2_ms_interval':[min(r['T2_ms_reciprocal'] for r in output if r['spectral_group']=='s3'),max(r['T2_ms_reciprocal'] for r in output if r['spectral_group']=='s3')]}
    (OUT/'NQR_T2_LOW_T_EXACT.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
    return summary


def script_inventory():
    scripts=[]
    folders=('data_analysis','powder_model','relaxation_analysis','research_literature',
             'new_evidence/structure_model','new_evidence/eta_explanation','draft')
    for folder in folders:
        for p in sorted((BASE/folder).glob('*.py')):
            tree=ast.parse(p.read_text());mods=set()
            for node in ast.walk(tree):
                if isinstance(node,ast.Import):mods.update(a.name.split('.')[0] for a in node.names)
                elif isinstance(node,ast.ImportFrom) and node.module:mods.add(node.module.split('.')[0])
            scripts.append(record(p,import_roots=sorted(mods)))
    return scripts


def main():
    originals=[record(p,role=ROLES[p.name],derived_outputs=DERIVED[p.name]) for p in sorted((BASE/'upload').iterdir()) if p.is_file()]
    assert len(originals)==10
    archives=archive_audit();assert all(a['all_extractions_identical'] for a in archives)
    integrals=table_info(BASE/'data_analysis/email_attachments/SpectrumIntegrals.dat','\t')
    exported=table_info(BASE/'data_analysis/export/july03_cus_integrals.csv')
    assert len(integrals)==len(exported)==51
    for a,b in zip(integrals,exported):
        for k in ('File','Freq','Integral'):assert float(a[k])==float(b[k])
    opj=json.loads((BASE/'data_analysis/export/opj_inventory.json').read_text())
    summary=json.loads((BASE/'data_analysis/export/data_summary.json').read_text())
    inventory={'audit_date':'2026-11-27','originals':originals,'archives':archives,
      'original_count':10,'original_total_bytes':sum(x['bytes'] for x in originals),
      'primary_opj':{'columns':len(opj['columns']),'worksheets':len(summary['opj_sheets']),
        'sheet_summary':summary['opj_sheets'],'source_sha256':opj['source_sha256']},
      'raw_trace_groups':[
        {'path':'data_analysis/anilite_raw','traces':600,'points_per_trace':1024,'time_step_us':0.16,
         'trace_columns':['Time(mks)','Real','Imag','Module'],'summary_tables':2,
         'frequency_associated_trace_ids':[31,600],'frequency_unassigned_trace_ids':[1,30],
         'frequency_range_in_summary_MHz':[45.1,67.0]},
        {'path':'data_analysis/cus_raw','traces':51,'trace_ids':[3,53],
         'trace_columns':['Time','Real','Imag','Module'],'frequency_range_MHz':[57.5,62.5],
         'frequency_mapping':'data_analysis/email_attachments/SpectrumIntegrals.dat'}],
      'raw_trace_shape_checks':[
        trace_shape_audit('data_analysis/anilite_raw','*_KN_ADC.dat'),
        trace_shape_audit('data_analysis/cus_raw','*.dat')],
      'SpectrumIntegrals_dat':{'path':'data_analysis/email_attachments/SpectrumIntegrals.dat',
         'sha256':digest((BASE/'data_analysis/email_attachments/SpectrumIntegrals.dat').read_bytes()),
         'rows':51,'frequency_range_MHz':[57.5,62.5],'frequency_step_MHz':0.1,
         'numeric_export_identical':True,'export_path':'data_analysis/export/july03_cus_integrals.csv',
         'different_from_OPJ_Spec05122016':True,'not_low_frequency_NQR':True,
         'rawSourceSpectrumIntegrals_exact_filename_present':False},
      'embedded_nqr':embedded_audit(),'nqr_low_temperature_T2':low_temperature_t2(),
      'low_frequency_1_5_to_2_2_MHz':{'identified_numeric_measurement_series_found':False,
         'archival_graph':'new_evidence/thesis/diploma5-p23.png','thesis_pages':[22,23,33],
         'not_digitized':True,'scope_note':'No identified1.5--2.2MHz table in the supplied scientific sources and14 embedded Origin projects;30 unassigned2016 traces remain frequency-unassigned.'},
      'known_data_issues':{'opj_sweep_audit':summary['opj_sweep_audit'],
         'july_amplitude_boundary':summary['july_amplitude_boundary'],
         'july_carried_forward_rows':90,
         'T2a_units':'TwoTau interpreted in microseconds by experimental context; legacy graphlabelmS preserved.'},
      'analysis_scripts':script_inventory(),
      'dependency_specification':record(BASE/'requirements.txt'),
      'export_files':[record(p) for p in sorted((BASE/'data_analysis/export').iterdir()) if p.is_file()],
      'duplicate_documents':[
         {'paths':['upload/Remarks.pdf','data_analysis/email_attachments/Remarks.pdf'],
          'byte_identical':(BASE/'upload/Remarks.pdf').read_bytes()==(BASE/'data_analysis/email_attachments/Remarks.pdf').read_bytes()},
         {'paths':['upload/geerite to digenite.docx','new_evidence/papers/raw/geerite to digenite.docx'],
          'byte_identical':(BASE/'upload/geerite to digenite.docx').read_bytes()==(BASE/'new_evidence/papers/raw/geerite to digenite.docx').read_bytes()}]}
    (OUT/'measurement_inventory.json').write_text(json.dumps(inventory,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'originals':10,'archive_members':sum(a['member_count'] for a in archives),
        'all_archive_extractions_byte_identical':True,'embedded_projects':14,
        'numeric_NQR_snapshots':5,'T2_low_temperature_rows':len(inventory['nqr_low_temperature_T2']['rows']),
        'output':'provenance/measurement_inventory.json'},indent=2))


if __name__=='__main__':main()
