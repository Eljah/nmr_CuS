from pathlib import Path
import re, subprocess, copy, json, csv, math
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_CELL_VERTICAL_ALIGNMENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

BASE=Path(__file__).resolve().parent.parent
SOURCE=BASE/'draft/manuscript_v2.md'
OUT=BASE/'output'
OUT.mkdir(exist_ok=True)
text=SOURCE.read_text(encoding='utf-8')
text=text.replace('[[FIT_RESULTS_SECTION]]',(BASE/'draft/powder_results_v2.md').read_text(encoding='utf-8'))
(BASE/'draft/manuscript_v2_complete.md').write_text(text,encoding='utf-8')
equations=re.findall(r'^\$\$(.*?)\$\$$',text,re.M)
equation_numbers=[re.search(r'\\tag\{(\d+)\}',e).group(1) for e in equations]
equations=[re.sub(r'\\tag\{\d+\}','',e) for e in equations]
math_source=BASE/'tmp/docx/equations.md'
math_source.parent.mkdir(parents=True,exist_ok=True)
math_source.write_text('\n\n'.join('$$'+e+'$$' for e in equations),encoding='utf-8')
math_doc=math_source.with_suffix('.docx')
subprocess.run(['pandoc',str(math_source),'-f','markdown+tex_math_dollars','-t','docx','-o',str(math_doc)],check=True)
eqdoc=Document(math_doc)
eq_paras=[p for p in eqdoc.element.body if p.xpath('.//m:oMath')]
assert len(eq_paras)==len(equations),(len(eq_paras),len(equations))

doc=Document()
section=doc.sections[0]
section.page_width=Inches(8.5); section.page_height=Inches(11)
section.top_margin=Inches(.78);section.bottom_margin=Inches(.72)
section.left_margin=Inches(.95);section.right_margin=Inches(.95)
section.header_distance=Inches(.28);section.footer_distance=Inches(.3)
styles=doc.styles
for name in ['Normal','Body Text','Title','Heading 1','Heading 2','Heading 3','Caption']:
    s=styles[name];s.font.name='Times New Roman';s.font.color.rgb=RGBColor(0,0,0)
    s.element.get_or_add_rPr().rFonts.set(qn('w:eastAsia'),'Times New Roman')
normal=styles['Normal']
normal.font.size=Pt(11)
normal.paragraph_format.line_spacing=1.12
normal.paragraph_format.space_after=Pt(5)
normal.paragraph_format.widow_control=True
styles['Title'].font.size=Pt(17)
styles['Title'].font.bold=True
styles['Title'].paragraph_format.space_after=Pt(12)
styles['Title'].paragraph_format.line_spacing=1.05
for name,size in [('Heading 1',13),('Heading 2',11.5),('Heading 3',11)]:
    st=styles[name];st.font.size=Pt(size);st.font.bold=True
    st.paragraph_format.space_before=Pt(13);st.paragraph_format.space_after=Pt(6)
    st.paragraph_format.keep_with_next=True
styles['Caption'].font.size=Pt(10)
styles['Caption'].font.italic=False
styles['Caption'].paragraph_format.space_before=Pt(4)
styles['Caption'].paragraph_format.space_after=Pt(8)
styles['Caption'].paragraph_format.line_spacing=1.03
styles['Caption'].paragraph_format.keep_with_next=False

def add_inline(p,s):
    for token in re.split(r'(ν_Q|ν_L|R_W)',s):
        if token in ('ν_Q','ν_L','R_W'):
            stem,sub=token.split('_')
            p.add_run(stem)
            run=p.add_run(sub);run.font.subscript=True
        else: p.add_run(token)

def set_cell_shade(cell,fill):
    pr=cell._tc.get_or_add_tcPr();node=OxmlElement('w:shd');node.set(qn('w:fill'),fill);pr.append(node)

def add_table(rows,index):
    widths=([1.55,4.05,1.0] if len(rows[0])==3 else [.8,1.6,1.6,2.6])
    table=doc.add_table(rows=1,cols=len(rows[0]));table.alignment=WD_TABLE_ALIGNMENT.CENTER
    table.autofit=False
    for col,w in zip(table.columns,widths):col.width=Inches(w)
    for i,row in enumerate(rows):
        cells=table.rows[0].cells if i==0 else table.add_row().cells
        for j,s in enumerate(row):
            c=cells[j];c.width=Inches(widths[j]);c.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER
            pr=c._tc.get_or_add_tcPr()
            mar=OxmlElement('w:tcMar')
            for side,value in [('top','90'),('bottom','90'),('left','100'),('right','100')]:
                z=OxmlElement('w:'+side);z.set(qn('w:w'),value);z.set(qn('w:type'),'dxa');mar.append(z)
            pr.append(mar)
            borders=OxmlElement('w:tcBorders')
            for side in ['top','left','bottom','right','insideH','insideV']:
                z=OxmlElement('w:'+side);z.set(qn('w:val'),'single');z.set(qn('w:sz'),'4');z.set(qn('w:color'),'D9D9D9');borders.append(z)
            pr.append(borders)
            p=c.paragraphs[0];p.paragraph_format.space_after=Pt(0);p.paragraph_format.line_spacing=1.02
            p.alignment=WD_ALIGN_PARAGRAPH.CENTER if (len(row)==4 or j!=1) else WD_ALIGN_PARAGRAPH.LEFT
            add_inline(p,s)
            for run in p.runs:run.font.size=Pt(10)
            if i==0:
                set_cell_shade(c,'3A3A3A')
                for run in p.runs:run.bold=True;run.font.color.rgb=RGBColor(255,255,255)
            else:set_cell_shade(c,'F3F3F3' if i%2==0 else 'FFFFFF')
        trpr=table.rows[i]._tr.get_or_add_trPr();cant=OxmlElement('w:cantSplit');trpr.append(cant)
        if i==0:
            z=OxmlElement('w:tblHeader');trpr.append(z)
        if i<len(rows)-1:
            for c in cells:
                for p in c.paragraphs:p.paragraph_format.keep_with_next=True
    p=doc.add_paragraph();p.paragraph_format.space_after=Pt(0);p.paragraph_format.space_before=Pt(0)
    p.paragraph_format.line_spacing=1;p.add_run().font.size=Pt(3)
    return table

def add_numbered_equation(el,number):
    if doc.paragraphs:
        doc.paragraphs[-1].paragraph_format.keep_with_next=True
    for math_para in list(el.xpath('.//m:oMathPara')):
        parent=math_para.getparent();position=parent.index(math_para)
        for math in list(math_para.findall(qn('m:oMath'))):
            math_para.remove(math);parent.insert(position,math);position+=1
        parent.remove(math_para)
    table=doc.add_table(rows=1,cols=3)
    table.alignment=WD_TABLE_ALIGNMENT.CENTER
    table.autofit=False
    for col,width in zip(table.columns,[.4,5.8,.4]):
        col.width=Inches(width)
    for cell,width in zip(table.rows[0].cells,[.4,5.8,.4]):
        cell.width=Inches(width)
        cell.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER
        pr=cell._tc.get_or_add_tcPr()
        borders=OxmlElement('w:tcBorders')
        for side in ['top','left','bottom','right','insideH','insideV']:
            border=OxmlElement('w:'+side);border.set(qn('w:val'),'nil');borders.append(border)
        pr.append(borders)
        mar=OxmlElement('w:tcMar')
        for side in ['left','right']:
            margin=OxmlElement('w:'+side);margin.set(qn('w:w'),'0');margin.set(qn('w:type'),'dxa');mar.append(margin)
        pr.append(mar)
        cell.paragraphs[0].paragraph_format.space_before=Pt(0)
        cell.paragraphs[0].paragraph_format.space_after=Pt(0)
    middle=table.cell(0,1)
    middle._tc.remove(middle.paragraphs[0]._p)
    middle._tc.append(el)
    p=middle.paragraphs[0]
    p.style='Normal'
    p.alignment=WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before=Pt(4)
    p.paragraph_format.space_after=Pt(4)
    p.paragraph_format.keep_together=True
    for math_para in el.xpath('.//m:oMathPara'):
        props=math_para.find(qn('m:oMathParaPr'))
        if props is None:
            props=OxmlElement('m:oMathParaPr');math_para.insert(0,props)
        jc=props.find(qn('m:jc'))
        if jc is None:jc=OxmlElement('m:jc');props.append(jc)
        jc.set(qn('m:val'),'center')
    label=table.cell(0,2).paragraphs[0]
    label.alignment=WD_ALIGN_PARAGRAPH.RIGHT
    label.add_run('('+number+')')
    trpr=table.rows[0]._tr.get_or_add_trPr();trpr.append(OxmlElement('w:cantSplit'))

lines=text.splitlines();i=0;eq_i=0;table_i=0;refs=False
while i<len(lines):
    line=lines[i].strip()
    if not line:i+=1;continue
    if line.startswith('# '):
        p=doc.add_paragraph(line[2:],style='Title')
        status=doc.add_paragraph('Редакция 2 для обсуждения с соавторами · 26 ноября 2026')
        status.paragraph_format.space_after=Pt(10)
        for run in status.runs:run.font.size=Pt(9)
    elif line.startswith('## '):
        heading=line[3:];p=doc.add_paragraph(heading,style='Heading 1')
        if heading=='Литература':refs=True
        elif heading=='English title and abstract':refs=False
    elif line.startswith('### '):
        doc.add_paragraph(line[4:],style='Heading 2')
    elif line in ('[[FIGURE_MAIN]]','[[FIGURE_POWDER]]','[[FIGURE_T2]]'):
        p=doc.add_paragraph();p.paragraph_format.keep_with_next=True
        p.paragraph_format.space_before=Pt(4);p.paragraph_format.space_after=Pt(1)
        figures={
            '[[FIGURE_MAIN]]':('data_analysis/export/figure_opj_manuscript.png','Исходные спектры ЯМР меди при периодах повторения 20 и 500 миллисекунд.'),
            '[[FIGURE_POWDER]]':('powder_model/powder_model_comparison.png','Совместная подгонка широкого порошкового отклика при 20 и 500 миллисекундах. Области металлической меди датчика маскированы. Показаны два сопоставимых условных решения.'),
            '[[FIGURE_T2]]':('relaxation_analysis/export/figure_t2_reanalysis.png','Затухание модуля комплексного интеграла эха и изменение фазы при 5 K, 51 МГц и 4,57 Тл. Эффективное время около 0,95 миллисекунды.')}
        filename,description=figures[line]
        p.add_run().add_picture(str(BASE/filename),width=Inches(6.3))
        p.alignment=WD_ALIGN_PARAGRAPH.CENTER
        for el in p._p.xpath('.//wp:docPr'):
            el.set('descr',description)
    elif line.startswith('$$'):
        el=copy.deepcopy(eq_paras[eq_i])
        add_numbered_equation(el,equation_numbers[eq_i]);eq_i+=1
    elif line.startswith('|'):
        rows=[]
        while i<len(lines) and lines[i].strip().startswith('|'):
            l=lines[i].strip()
            if not re.fullmatch(r'[|:\s-]+',l):rows.append([c.strip() for c in l.strip('|').split('|')])
            i+=1
        assert all(len(r)==len(rows[0]) for r in rows),rows
        add_table(rows,table_i);table_i+=1;continue
    else:
        if line.startswith('Рисунок '):style='Caption'
        else:style='Normal'
        p=doc.add_paragraph(style=style);add_inline(p,line)
        if line.startswith('Таблица '):
            p.paragraph_format.keep_with_next=True
            p.paragraph_format.space_before=Pt(7)
            p.paragraph_format.space_after=Pt(6)
            for run in p.runs:run.font.size=Pt(10)
        elif style=='Normal' and not refs:
            p.alignment=WD_ALIGN_PARAGRAPH.JUSTIFY
        if refs:
            p.paragraph_format.line_spacing=1.04;p.paragraph_format.space_after=Pt(5)
            p.paragraph_format.keep_together=True
            for run in p.runs:run.font.size=Pt(9.5)
    i+=1
footer=section.footer.paragraphs[0];footer.alignment=WD_ALIGN_PARAGRAPH.CENTER
run=footer.add_run();fld=OxmlElement('w:fldSimple');fld.set(qn('w:instr'),'PAGE');run._r.addnext(fld)
doc.core_properties.title='Порошковые спектры ЯМР меди в анилите Cu7S4: квадрупольные параметры, затухание эха и связь с прежними данными ЯКР'
doc.core_properties.subject='Рукопись для научного обсуждения и проверки экспериментального протокола'
doc.core_properties.author=''
doc.core_properties.keywords='anilite, Cu7S4, NMR, NQR'

comments=[('Порошковые спектры ЯМР меди', 'Редакция 2 включает новые авторские материалы, совместную подгонку порошкового спектра и обработку T2a. Состав авторов, аффилиации, журнал и благодарности пока не назначены.'), ('Перед измерениями ЯМР выполнена повторная рентгенография', 'Повторная идентификация как анилита и общность образца ЯКР/ЯМР сообщены автором исследования. Общность образца больше не является открытым вопросом. Для окончательной статьи добавить сам повторный XRD, дату и результат фазового анализа.'), ('Основной набор включает спектры', 'Добавить аппаратные параметры основной серии: спектрометр, импульсы, tau, накопления, усиление, калибровку поля и метод интегрирования. Параметры ранней отдельной серии из Remarks.pdf не переносились автоматически.'), ('Для численного сравнения используется таблица', 'Подписи 5 K и 67.20566 кЭ связаны непосредственно с графиком обработанных спектров. Для табличных интегралов Spec05122016 отдельный журнал не найден. Подтвердить соответствие таблицы и графика при подготовке окончательного экспериментального раздела; исходные файлы 425/426 нужны только для разрешения конфликта 81.40 МГц.'), ('Основной набор ограничений взят', 'Частоты диплома использованы как основной собственный набор. Статья 2013 содержит иной набор, а ранняя английская рукопись относит его во введении к Anashkin 1994. Проверка чувствительности выполнена; при публикации следует согласовать окончательные частоты с исходными ЯКР спектрами и учесть разницу 4.2/5 K.'), ('Расчёт различает модель всех переходов', 'Модели, критерии и численные результаты сохранены в powder_model. Фит выполнен, но его параметры остаются условными: анизотропия магнитного сдвига и полный импульсный отклик пока не идентифицированы. Сводная таблица не является присвоением η позициям Cu1–Cu5.'), ('В настоящем анализе TwoTau интерпретируется', 'Микросекундный масштаб восстановлен по диплому и пояснению автора. В формуле используется полное время эха t=2tau. Перед окончательной подачей желательно сверить технический масштаб экспорта; прежняя подпись mS сохранена в метаданных исходного проекта.')]

for needle,comment in comments:
    found=[p for p in doc.paragraphs if needle in p.text]
    assert found,needle
    doc.add_comment(found[0].runs,text=comment,author='Редакционные вопросы',initials='РВ')
path=OUT/'Anilite_NMR_draft_ru.docx'
doc.save(path)
print(json.dumps({'docx':str(path),'equations':eq_i,'tables':table_i,'comments':len(comments),'paragraphs':len(doc.paragraphs)},ensure_ascii=False))
