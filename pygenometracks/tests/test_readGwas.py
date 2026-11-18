import os
from tempfile import NamedTemporaryFile

from pygenometracks.readGwas import ReadGwas

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")


def test_read_gwas_4col():
    gwas = ReadGwas(os.path.join(ROOT, 'gwas_1.gwas'))
    assert not gwas.has_header
    assert gwas.used_fields == ['chromosome', 'position', 'variant_id', 'pvalue']
    assert gwas.req_f_pos == [0, 1, 3]
    records = [rec for rec in gwas]
    assert records[0].position == 3002145
    assert records[1].variant_id == 'rs6610293'
    assert records[2].pvalue == 9.2e-5
    assert records[-1].variant_id == 'rs8840123'
    assert set([r.chromosome for r in records]) == {'X'}


def test_read_gwas_header():
    gwas = ReadGwas(os.path.join(ROOT, 'gwas_2.gwas'), has_header=True)
    assert gwas.has_header
    assert gwas.used_fields == ['chromosome', 'position', 'variant_id', 'pvalue', 'BETA', 'SE', 'MAF']
    assert gwas.req_f_pos == [0, 1, 3]
    records = [rec for rec in gwas]
    assert records[0].position == 3001200
    assert records[1].variant_id == 'rs100002'
    assert records[2].pvalue == 1.1e-6
    assert records[3].BETA == '0.33'
    assert records[4].SE == '0.11'
    assert records[-1].MAF == '0.21'
    assert set([r.chromosome for r in records]) == {'X'}


def test_read_gwas_glm_linear():
    gwas = ReadGwas(os.path.join(ROOT, 'head_all_hg38_qcd_LE1_simgwas_quant1a.simgwas_quant1.glm.linear'), has_header=True)
    assert gwas.has_header
    assert gwas.used_fields == ['chromosome', 'position', 'variant_id', 'REF', 'ALT', 'PROVISIONAL_REF_', 'A1', 'OMITTED', 'A1_FREQ', 'TEST', 'OBS_CT', 'BETA', 'SE', 'T_STAT', 'pvalue', 'ERRCODE']
    assert gwas.req_f_pos == [0, 1, 14]
    records = [rec for rec in gwas]
    assert records[0].position == 54490
    assert records[1].variant_id == '1:58176G,A'
    assert records[2].pvalue == 1.70492e-14
    assert records[3].BETA == '-0.285574'
    assert records[4].SE == '0.0433716'
    assert records[-1].REF == 'G'
    assert set([r.chromosome for r in records]) == {'1'}


def test_read_gwas_OR_colname():
    infile = NamedTemporaryFile(suffix='.gwas', prefix='gwas_test_',
                                delete=False)
    with open(infile.name, 'w') as fo:
        fo.write('chrom\tbp\tp\tOR\tor\nchr1\t12\t0.001\twhatever\tsomething\n')
    gwas = ReadGwas(infile.name, has_header=True)
    assert "OR" in gwas.used_fields
    assert "or" not in gwas.used_fields
    records = [rec for rec in gwas]
    assert records[0].OR == "whatever"
    # The 'or' colnames is changed to _4:
    assert records[0]._4 == "something"
