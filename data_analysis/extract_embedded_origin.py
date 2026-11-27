"""Read embedded Origin projects as data; never activate OLE applications.

The small CFB reader implements sector chains and directory entries used by the
supplied Word/Origin objects. It is intentionally limited, with explicit bounds
and cycle checks. Original container bytes and provenance are preserved.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import re
import struct
import math
import zipfile

BASE = Path(__file__).resolve().parent.parent
OUT = BASE / 'data_analysis' / 'thesis_embedded'
END, FREE = 0xfffffffe, 0xffffffff


def sha256(b):
    return hashlib.sha256(b).hexdigest()


class CFB:
    def __init__(self, data):
        self.data = data
        assert data[:8] == bytes.fromhex('d0cf11e0a1b11ae1')
        assert data[28:30] == b'\xfe\xff'
        self.major = struct.unpack_from('<H', data, 26)[0]
        self.sector_size = 1 << struct.unpack_from('<H', data, 30)[0]
        self.mini_size = 1 << struct.unpack_from('<H', data, 32)[0]
        assert self.major == 3 and self.sector_size == 512 and self.mini_size == 64
        assert len(data) % self.sector_size == 0
        nfat, first_dir = struct.unpack_from('<II', data, 44)
        self.cutoff, first_minifat, nminifat, first_difat, ndifat = struct.unpack_from('<IIIII', data, 56)
        difat = list(struct.unpack_from('<109I', data, 76))
        sector = first_difat
        visited = set()
        for _ in range(ndifat):
            assert sector not in visited
            visited.add(sector)
            vals = struct.unpack('<128I', self.sector(sector))
            difat.extend(vals[:-1]); sector = vals[-1]
        fat_ids = [i for i in difat if i != FREE]
        assert len(fat_ids) >= nfat
        self.fat = []
        for sector in fat_ids[:nfat]:
            self.fat.extend(struct.unpack('<128I', self.sector(sector)))
        raw_dir = self.chain(first_dir, self.fat, self.sector)
        self.entries = []
        for offset in range(0, len(raw_dir), 128):
            d = raw_dir[offset:offset+128]
            size_name = struct.unpack_from('<H', d, 64)[0]
            kind = d[66]
            name = d[:max(0,size_name-2)].decode('utf-16le') if kind else ''
            self.entries.append(dict(index=offset//128, name=name, kind=kind,
                left=struct.unpack_from('<I',d,68)[0], right=struct.unpack_from('<I',d,72)[0],
                child=struct.unpack_from('<I',d,76)[0], start=struct.unpack_from('<I',d,116)[0],
                size=struct.unpack_from('<Q',d,120)[0]))
        root = self.entries[0]
        assert root['kind'] == 5
        self.mini_stream = self.chain(root['start'],self.fat,self.sector)[:root['size']] if root['size'] else b''
        raw_minifat = self.chain(first_minifat,self.fat,self.sector) if nminifat else b''
        assert len(raw_minifat) == nminifat*self.sector_size
        self.minifat = list(struct.unpack('<'+'I'*(len(raw_minifat)//4),raw_minifat))

    def sector(self, index):
        assert 0 <= index < len(self.data)//self.sector_size-1, index
        start=(index+1)*self.sector_size
        return self.data[start:start+self.sector_size]

    @staticmethod
    def chain(start, table, reader):
        chunks=[]; visited=set(); current=start
        while current != END:
            assert current != FREE and current not in visited and 0 <= current < len(table), current
            visited.add(current); chunks.append(reader(current)); current=table[current]
        return b''.join(chunks)

    def stream(self, entry):
        if not entry['size']: return b''
        if entry['size'] >= self.cutoff:
            b=self.chain(entry['start'],self.fat,self.sector)
        else:
            def mini(i):
                start=i*self.mini_size
                assert start+self.mini_size <= len(self.mini_stream)
                return self.mini_stream[start:start+self.mini_size]
            b=self.chain(entry['start'],self.minifat,mini)
        assert len(b) >= entry['size']
        return b[:entry['size']]

    def streams(self):
        seen=set()
        def walk(index,parent):
            if index == FREE: return
            assert index not in seen and 0 <= index < len(self.entries)
            seen.add(index); entry=self.entries[index]
            yield from walk(entry['left'],parent)
            path=parent+[entry['name']]
            if entry['kind'] == 1:
                yield from walk(entry['child'],path)
            elif entry['kind'] == 2:
                yield '/'.join(path),self.stream(entry)
            yield from walk(entry['right'],parent)
        yield from walk(self.entries[0]['child'],[])


def write_file(path, data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_bytes(data)
    return {'path':str(path.relative_to(BASE)), 'bytes':len(data), 'sha256':sha256(data)}


def inspect_origin_streams(container, label, parent_provenance):
    records=[]
    for stream_path, data in CFB(container).streams():
        if data.startswith(b'CPYA'):
            filename=OUT/f'{label}_{len(records)+1:02}.opj'
            records.append({'container_stream_path':stream_path,
                'signature':data.split(b'\n',1)[0].decode('latin1'),
                **write_file(filename,data), **parent_provenance})
    return records


def extract():
    OUT.mkdir(parents=True,exist_ok=True)
    records=[]; containers=[]
    source=BASE/'upload/diploma_magister_geerite_2009_5.doc'
    doc=source.read_bytes()
    records.extend(inspect_origin_streams(doc,'diploma5',
        {'source':str(source.relative_to(BASE)),'source_sha256':sha256(doc)}))
    source=BASE/'new_evidence/papers/raw/diploma_magister_geerite_2009_10.rtf'
    rtf=source.read_bytes()
    pattern=rb'\\objclass Origin50.Graph\}\{\\\*\\objdata\s+([0-9a-fA-F\s]+)\}'
    for index,match in enumerate(re.finditer(pattern,rtf),1):
        ole1=bytes.fromhex(match.group(1).decode('ascii'))
        assert struct.unpack_from('<II',ole1,0) == (0x501,2)
        p=8; fields=[]
        for _ in range(3):
            n=struct.unpack_from('<I',ole1,p)[0];p+=4
            fields.append(ole1[p:p+n]);p+=n
        assert fields[0].rstrip(b'\0') == b'Origin50.Graph'
        n=struct.unpack_from('<I',ole1,p)[0];p+=4
        native=ole1[p:p+n];assert len(native)==n
        container_meta=write_file(OUT/f'diploma10_object_{index:02}.cfb',native)
        container_meta.update({'source':str(source.relative_to(BASE)),
            'source_sha256':sha256(rtf),'rtf_objclass':'Origin50.Graph',
            'rtf_match_byte_offset':match.start(),'rtf_hex_byte_offset':match.start(1),
            'ole1_native_byte_offset':p,'ole1_native_byte_length':n,
            'ole1_full':write_file(OUT/f'diploma10_object_{index:02}.ole1',ole1)})
        containers.append(container_meta)
        records.extend(inspect_origin_streams(native,f'diploma10_object_{index:02}',
            {'source':str(source.relative_to(BASE)),'source_sha256':sha256(rtf),
             'container_file':container_meta['path'],'container_sha256':container_meta['sha256'],
             'rtf_match_byte_offset':match.start()}))
    docx_containers=[]
    source=BASE/'upload/geerite to digenite.docx'
    source_bytes=source.read_bytes()
    with zipfile.ZipFile(source) as z:
        for name in sorted(z.namelist()):
            if '/embeddings/' not in name:continue
            native=z.read(name)
            # Inspect streams, not the object application's displayed preview.
            streams=list(CFB(native).streams())
            if not any(d.startswith(b'CPYA') for _,d in streams):continue
            label='geerite_to_digenite_'+Path(name).stem
            container_meta=write_file(OUT/f'{label}.cfb',native)
            container_meta.update({'source':str(source.relative_to(BASE)),
                'source_sha256':sha256(source_bytes),'docx_zip_member':name})
            docx_containers.append(container_meta)
            records.extend(inspect_origin_streams(native,label,{
                'source':str(source.relative_to(BASE)),'source_sha256':sha256(source_bytes),
                'container_file':container_meta['path'],'container_sha256':container_meta['sha256'],
                'docx_zip_member':name}))
    result={'mode':'read_only_data_extraction_no_OLE_execution','origin_projects':records,
            'rtf_native_containers':containers,'docx_native_containers':docx_containers}
    (OUT/'embedded_origin_manifest.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    return result


def read_origin_714(path):
    """Parse the 131-byte data-list headers in these Origin 7.5-era objects.

    The common block/data-cell layout follows the public Origin file-format
    description already used by analyze_measurements.py. Unlike that reader,
    this version retains text cells and checks the distinct signature/header
    length actually present in the embedded thesis projects.
    """
    data=path.read_bytes(); signature=data.split(b'\n',1)[0]
    assert signature == b'CPYA 4,2766 714#',signature
    cursor=len(signature)+1
    def block():
        nonlocal cursor
        start=cursor
        assert cursor+5 <= len(data)
        n=struct.unpack_from('<I',data,cursor)[0]
        assert data[cursor+4]==10
        cursor+=5;payload=data[cursor:cursor+n];cursor+=n
        assert len(payload)==n
        if n:
            assert cursor<len(data) and data[cursor]==10
            cursor+=1
        return start,payload
    _,global_header=block();assert len(global_header)==39
    _,empty=block();assert empty==b''
    columns=[]
    while True:
        start,h=block()
        if not h:break
        assert len(h)==131,(start,len(h))
        name=h[0x58:0x71].split(b'\0',1)[0].decode('cp1251')
        value_size=h[0x3d]
        dtype=struct.unpack_from('<H',h,0x16)[0]
        n,first,last=struct.unpack_from('<III',h,0x19)
        data_offset,payload=block();_,mask=block()
        assert value_size==10 and len(payload)==n*10 and not mask,(name,value_size,len(payload),n)
        values=[];text_cells=[];finite=[];missing=[]
        for i in range(n):
            cell=payload[i*10:(i+1)*10]
            flag=struct.unpack_from('<H',cell)[0]
            if flag==0:
                value=struct.unpack_from('<d',cell,2)[0]
                if value == -1.23456789e-300:
                    value=None;missing.append(i+1)
                else:
                    assert math.isfinite(value),(name,i,value)
                    finite.append(value)
            elif flag==1:
                value=cell[2:].split(b'\0',1)[0].decode('cp1251')
                text_cells.append({'origin_row':i+1,'text':value,'raw_cell_hex':cell.hex()})
            else:
                raise ValueError(f'Unsupported cell flag {flag}: {name} row {i+1}')
            values.append(value)
        columns.append({'name':name,'header_offset':start,'content_block_offset':data_offset,
          'numeric_cell_payload_offset':data_offset+5,'dtype_hex':hex(dtype),'value_size':value_size,
          'allocated_rows':n,'header_first_row':first,'header_last_row':last,
          'finite_numeric_rows':len(finite),'missing_rows':len(missing),'text_cells':text_cells,
          'numeric_min':min(finite) if finite else None,'numeric_max':max(finite) if finite else None,
          'values':values})
    return {'source_project':str(path.relative_to(BASE)),'source_sha256':sha256(data),
        'signature':signature.decode('latin1'),'data_list_end_offset':cursor,'columns':columns}


def graph_first_layers(path, names):
    """Recover graph labels and data references, without executing formulas."""
    b=path.read_bytes();candidates=[]
    for match in re.finditer(rb'\x00\x00(Graph[A-Za-z0-9_]{0,18})\x00',b):
        pos=match.start()-5
        if pos<0:continue
        n=struct.unpack_from('<I',b,pos)[0]
        if 195<=n<2048 and pos+5+n<len(b) and b[pos+4]==10 and b[pos+5+n]==10:
            candidates.append((match.group(1).decode('ascii'),pos))
    def block(p):
        n=struct.unpack_from('<I',b,p)[0];assert b[p+4]==10
        d=b[p+5:p+5+n];assert len(d)==n
        if n:assert b[p+5+n]==10
        return p+5+n+(1 if n else 0),d
    result=[]
    for name,offset in candidates:
        p,h=block(offset);p,layer=block(p)
        item={'name':name,'window_header_offset':offset,'parsed_layer':1,
            'x_axis_limits_raw':list(struct.unpack_from('<ddd',layer,0x0f)),
            'annotations':[],'curves':[]}
        while True:
            off=p;p,h=block(p)
            if not h:break
            p,d1=block(p);p,d2=block(p);p,d3=block(p)
            section=h[0x46:0x46+41].split(b'\0',1)[0].decode('cp1251')
            def text(d):
                raw=d.split(b'\0',1)[0]
                return raw.decode('cp1251') if all(c>=32 or c in (9,10,13) for c in raw) else None
            item['annotations'].append({'name':section,'header_offset':off,
                'data1_text':text(d1),'data2_text':text(d2),'data3_text':text(d3)})
        while True:
            off=p;p,c=block(p)
            if not c:break
            p,d=block(p)
            y=struct.unpack_from('<H',c,4)[0];x=struct.unpack_from('<H',c,0x23)[0]
            assert 0<=x<=len(names) and 0<=y<=len(names),(name,x,y)
            item['curves'].append({'curve_index_1_based':len(item['curves'])+1,
                'header_offset':off,'header_bytes':len(c),
                'x_dataset_index_1_based':x,'y_dataset_index_1_based':y,
                'x_dataset':names[x-1] if x else None,'y_dataset':names[y-1] if y else None,
                'formula_or_label_stored_not_executed':d.decode('cp1251')})
        result.append(item)
    return result


def export_numeric_projects(result):
    summaries=[]
    for item in result['origin_projects']:
        path=BASE/item['path'];parsed=read_origin_714(path)
        parsed['graph_first_layers']=graph_first_layers(path,[c['name'] for c in parsed['columns']])
        target=OUT/'numeric'/path.stem;target.mkdir(parents=True,exist_ok=True)
        sheets={}
        for c in parsed['columns']:
            sheet,column=c['name'].rsplit('_',1)
            sheets.setdefault(sheet,{})[column]=c['values']
        csv_files=[]
        for sheet,columns in sheets.items():
            assert re.fullmatch(r'[A-Za-z0-9_-]+',sheet),sheet
            names=list(columns);n=max(len(v) for v in columns.values())
            filename=target/f'{sheet}.csv'
            with filename.open('w',newline='',encoding='utf-8') as f:
                writer=csv.writer(f);writer.writerow(['origin_row']+names)
                for i in range(n):
                    row=[columns[name][i] if i<len(columns[name]) else None for name in names]
                    if any(v is not None for v in row):
                        writer.writerow([i+1]+['' if v is None else format(v,'.17g') if isinstance(v,float) else v for v in row])
            csv_files.append({'sheet':sheet,'path':str(filename.relative_to(BASE)),
                'bytes':filename.stat().st_size,'sha256':sha256(filename.read_bytes())})
        parsed['csv_files']=csv_files
        (target/'column_provenance.json').write_text(json.dumps(parsed,ensure_ascii=False,indent=2)+'\n')
        summaries.append({'source':item['path'],'source_sha256':item['sha256'],
            'columns':[{k:v for k,v in c.items() if k!='values'} for c in parsed['columns']],
            'csv_files':csv_files})
    (OUT/'numeric_inventory.json').write_text(json.dumps(summaries,ensure_ascii=False,indent=2)+'\n')
    return summaries


if __name__ == '__main__':
    result=extract()
    summaries=export_numeric_projects(result)
    print(json.dumps([{'path':r['source'],'columns':len(r['columns']),
        'sheets':[f['sheet'] for f in r['csv_files']]} for r in summaries],indent=2))
