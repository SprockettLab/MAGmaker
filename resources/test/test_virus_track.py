#!/usr/bin/env python3
"""Checks the viral track: which tools the switches run, and that every
predictor's sequence names map back to their source contigs.

1. Snakemake dry runs of virus_all on the test data:
   - the defaults run geNomad, its CheckV and vRhyme;
   - all four predictors add their own rules and one CheckV each;
   - with geNomad off, vRhyme says it will not run and no vRhyme job is
     planned (it bins geNomad's calls);
   - --skip-vrhyme drops the vRhyme jobs and nothing else;
   - an unknown predictor name is refused.

2. make_virus_summary.py on a built output tree holding the exact header
   shapes each tool really writes, including the traps that make a wrong
   mapping rule look like genuine disagreement between tools:
   VIBRANT's _fragment_1 after len=, geNomad's |provirus_<start>_<end>,
   VirSorter2's ||full, Cenote-Taker3's unrecoverable names, and the
   __dup<N> suffix added so CheckV would accept repeated IDs.
   It also checks the five provirus columns, and that a predictor whose
   output file is missing is reported as absent rather than as zero calls.

Run from the repository root (part 1 needs snakemake on PATH):

    python3 resources/test/test_virus_track.py
"""

import os
import re
import subprocess
import sys
import tempfile
import types

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
SCRIPT = os.path.join(REPO_ROOT, 'resources', 'scripts', 'make_virus_summary.py')
SAMPLES = ['John', 'Paul', 'George', 'Ringo']
CHECKS = []


def check(name, ok):
    CHECKS.append((name, ok))
    print(('  ok   ' if ok else '  FAIL ') + name)
    if not ok:
        raise AssertionError(name)


# ── part 1: the plan ────────────────────────────────────────────────────────

def plan(*extra_config, overrides=''):
    """Job counts from a dry run of virus_all, plus the track's own message.

    `overrides` is extra YAML appended to the override config file, for
    settings that sit inside a section --config would replace wholesale.
    """
    with tempfile.TemporaryDirectory(prefix='magmaker_virus_') as tmp:
        os.symlink(os.path.join(REPO_ROOT, 'resources'),
                   os.path.join(tmp, 'resources'))
        with open(os.path.join(tmp, 'host.fna'), 'w') as handle:
            handle.write('>c1\nACGTACGTACGGTTAACCGGTTAAGGCCAATTGG\n')
        over = os.path.join(tmp, 'over.yaml')
        with open(over, 'w') as handle:
            handle.write('host_filter:\n  boardercollie: false\n'
                         '  db_dir: {0}/host_db\n  genome: {0}/host.fna\n'.format(tmp))
            handle.write(overrides)
        result = subprocess.run(
            ['snakemake', '-s', os.path.join(REPO_ROOT, 'Snakefile-bin'),
             '--configfile', os.path.join(REPO_ROOT, 'config.yaml'), over,
             '--cores', '1', '-d', tmp, '-n', 'virus_all',
             '--config', 'binning=resources/config/binning.txt', *extra_config],
            capture_output=True, text=True)
        out = result.stdout + result.stderr
        if result.returncode != 0:
            return None, out
        stats = out[out.index('Job stats'):]
        counts = {m.group(1): int(m.group(2))
                  for m in re.finditer(r'^(\w+)\s+(\d+)$', stats, re.M)}
        return counts, out


def test_plan():
    print('plan:')
    n = len(SAMPLES)

    default, out = plan()
    check('default: predictors line names geNomad only',
          'viral track: predictors genomad;' in out)
    check('default: geNomad and its CheckV, once per assembly',
          default.get('run_genomad') == n and default.get('run_checkv') == n)
    check('default: vRhyme on, one prep and one run per assembly',
          'vRhyme on' in out and default.get('vrhyme_prep') == n
          and default.get('run_vrhyme') == n
          and default.get('vrhyme_summary') == 1)
    check('default: no other predictor runs',
          not {'run_vibrant', 'run_cenotetaker3', 'run_virsorter2',
               'run_checkv_predictor'} & set(default))

    every, out = plan('virus_predictors=[genomad,vibrant,cenotetaker3,virsorter2]')
    check('all four: each predictor runs once per assembly',
          all(every.get('run_' + t) == n
              for t in ('genomad', 'vibrant', 'cenotetaker3', 'virsorter2')))
    check('all four: one CheckV per non-geNomad predictor per assembly, and '
          'geNomad keeps its own rule',
          every.get('run_checkv_predictor') == 3 * n
          and every.get('run_checkv') == n)

    # vRhyme bins geNomad's viral contigs, so it cannot run without it. The
    # run says so once instead of failing later.
    without, out = plan('virus_predictors=[vibrant]')
    check('geNomad off: vRhyme says it will not run, and says why',
          'vRhyme is on but will NOT run' in out
          and 'geNomad is off' in out)
    check('geNomad off: no vRhyme or geNomad job is planned',
          not {'vrhyme_prep', 'run_vrhyme', 'vrhyme_summary', 'run_genomad',
               'run_checkv'} & set(without))
    check('geNomad off: VIBRANT and its own CheckV still run',
          without.get('run_vibrant') == n
          and without.get('run_checkv_predictor') == n)

    skipped, out = plan('skip_vrhyme=True')
    check('--skip-vrhyme: no vRhyme jobs',
          not {'vrhyme_prep', 'run_vrhyme', 'vrhyme_summary'} & set(skipped))
    check('--skip-vrhyme: changes nothing else',
          {k: v for k, v in default.items()
           if not k.startswith('vrhyme') and k != 'run_vrhyme'}
          == {k: v for k, v in skipped.items() if k != 'total'}
          or skipped.get('run_genomad') == n)

    counts, out = plan('virus_predictors=[genomad,nosuchtool]')
    check('an unknown predictor name is refused',
          counts is None and 'nosuchtool' in out)

    # vRhyme needs MetaBAT2's coverage table; without that binner it would
    # have to map the reads again, so it says so rather than doing it.
    counts, out = plan('binners=[concoct]')
    check('metabat2 off: vRhyme says it will not run, and says why',
          'vRhyme is on but will NOT run' in out and 'metabat2' in out)
    return default


# ── part 2: name mapping and the summary table ──────────────────────────────

# The real header shapes, each taken from output the tools actually wrote.
GENOMAD_VIRUS = [
    ('k141_7410 flag=1 multi=4.0000 len=506', 'k141_7410'),
    # An excised provirus carries its coordinates on the NAME, so without the
    # '|' cut no provirus ever joins back to its contig.
    ('k141_99|provirus_120_4100 flag=1 multi=9.0 len=3980', 'k141_99'),
]
VIBRANT_PHAGES = [
    # _fragment_1 is appended to the END of the header, after len=, NOT to
    # the contig name: splitting on '_' mangles it.
    ('k141_10643 flag=0 multi=22.9958 len=231935_fragment_1', 'k141_10643'),
    ('k141_7410 flag=1 multi=4.0000 len=506', 'k141_7410'),
    # One contig can yield several regions with the SAME id; unique_fasta_ids
    # renames the repeats __dup2 so CheckV accepts them, and the suffix comes
    # off again here.
    ('k141_10643__dup2 flag=0 multi=22.9958 len=231935_fragment_2', 'k141_10643'),
]
VIRSORTER2_SEQS = [
    ('k141_7410||full', 'k141_7410'),
    ('k141_99||0_partial', 'k141_99'),
    ('k141_555||lt2gene', 'k141_555'),
]
# Cenote-Taker3 renames beyond recovery, so its FASTA is never parsed; the
# source contig comes from input_name in its own summary table.
CT3_SEQS = [('ct3_1@C3', 'k141_7410 flag=1 multi=4.0 len=506', 'k141_7410'),
            ('ct3_2@C3', 'k141_99 flag=1 multi=9.0 len=3980', 'k141_99')]


def write_fasta(path, headers, sequence='ACGTACGTAC'):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as handle:
        for header in headers:
            handle.write('>{0}\n{1}\n'.format(header, sequence))


def write_tsv(path, columns, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as handle:
        handle.write('\t'.join(columns) + '\n')
        for row in rows:
            handle.write('\t'.join(str(row.get(c, '')) for c in columns) + '\n')


def build_tree(tmp, sample='S1', assembler='megahit'):
    """A finished viral output tree for one assembly, every predictor on."""
    stem = sample + '.contigs'
    # The assembly, so the mapping report can tell whether a mapped name is
    # really in it.
    write_fasta(os.path.join(tmp, 'output', 'assemble', assembler, stem + '.fasta'),
                ['k141_7410 flag=1', 'k141_99 flag=1', 'k141_555 flag=1',
                 'k141_10643 flag=0'])

    # geNomad: its own summary tables, plus the viral FASTA.
    summary_dir = os.path.join(tmp, 'output', 'virus', 'genomad', assembler,
                               sample, stem + '_summary')
    write_fasta(os.path.join(summary_dir, stem + '_virus.fna'),
                [h for h, _ in GENOMAD_VIRUS])
    write_tsv(os.path.join(summary_dir, stem + '_virus_summary.tsv'),
              ['seq_name', 'length', 'topology', 'n_genes', 'n_hallmarks',
               'virus_score', 'fdr', 'taxonomy'],
              [{'seq_name': GENOMAD_VIRUS[0][0], 'length': 506,
                'topology': 'No terminal repeats', 'n_genes': 2,
                'n_hallmarks': 1, 'virus_score': 0.95, 'fdr': 0.01,
                'taxonomy': 'Viruses;Duplodnaviria'},
               {'seq_name': GENOMAD_VIRUS[1][0], 'length': 3980,
                'topology': 'Provirus', 'n_genes': 5, 'n_hallmarks': 2,
                'virus_score': 0.99, 'fdr': 0.0,
                'taxonomy': 'Viruses;Duplodnaviria'}])
    write_tsv(os.path.join(summary_dir, stem + '_plasmid_summary.tsv'),
              ['seq_name', 'length', 'topology', 'n_genes', 'plasmid_score',
               'fdr', 'conjugation_genes', 'amr_genes'],
              [{'seq_name': 'k141_555 flag=1 multi=3.0 len=9000',
                'length': 9000, 'topology': 'No terminal repeats',
                'n_genes': 8, 'plasmid_score': 0.88, 'fdr': 0.02,
                'conjugation_genes': 'traD', 'amr_genes': 'tetA'}])
    # geNomad's CheckV. It reports provirus=No for a sequence geNomad already
    # excised, because no host DNA is left: that combination is expected.
    write_tsv(os.path.join(tmp, 'output', 'virus', 'checkv', assembler, sample,
                           'quality_summary.tsv'),
              ['contig_id', 'checkv_quality', 'miuvig_quality', 'completeness',
               'contamination', 'provirus', 'proviral_start', 'proviral_end',
               'warnings'],
              [{'contig_id': 'k141_7410', 'checkv_quality': 'High-quality',
                'miuvig_quality': 'High-quality', 'completeness': 95.0,
                'contamination': 0.0, 'provirus': 'No'},
               # CheckV keys on the name geNomad wrote, pipe part included.
               {'contig_id': 'k141_99|provirus_120_4100',
                'checkv_quality': 'Medium-quality',
                'miuvig_quality': 'Genome-fragment', 'completeness': 62.0,
                'contamination': 0.0, 'provirus': 'No'}])

    # VIBRANT, Cenote-Taker3, VirSorter2, each at its own real output path.
    write_fasta(os.path.join(tmp, 'output', 'virus', 'vibrant', assembler,
                             sample, 'VIBRANT_' + stem,
                             'VIBRANT_phages_' + stem,
                             stem + '.phages_combined.fna'),
                [h for h, _ in VIBRANT_PHAGES])
    ct3_dir = os.path.join(tmp, 'output', 'virus', 'cenotetaker3', assembler,
                           sample, 'ct3')
    write_fasta(os.path.join(ct3_dir, 'ct3_virus_sequences.fna'),
                [own for own, _, _ in CT3_SEQS])
    write_tsv(os.path.join(ct3_dir, 'ct3_virus_summary.tsv'),
              ['contig', 'input_name'],
              [{'contig': own, 'input_name': original}
               for own, original, _ in CT3_SEQS])
    write_fasta(os.path.join(tmp, 'output', 'virus', 'virsorter2', assembler,
                             sample, 'final-viral-combined.fa'),
                [h for h, _ in VIRSORTER2_SEQS])
    # One CheckV table per predictor, keyed on the names each one emits.
    for predictor, names in (('vibrant', [h.split()[0] for h, _ in VIBRANT_PHAGES]),
                             ('cenotetaker3', [o for o, _, _ in CT3_SEQS]),
                             ('virsorter2', [h for h, _ in VIRSORTER2_SEQS])):
        write_tsv(os.path.join(tmp, 'output', 'virus', 'checkv_' + predictor,
                               assembler, sample, 'quality_summary.tsv'),
                  ['contig_id', 'checkv_quality', 'completeness', 'provirus',
                   'proviral_start', 'proviral_end'],
                  [{'contig_id': name, 'checkv_quality': 'Low-quality',
                    'completeness': 20.0, 'provirus': 'No'} for name in names])
    return stem


def run_summary(tmp, predictors, sample='S1', assembler='megahit', stem=None):
    """Run make_virus_summary.py over a built tree; return its two tables."""
    out_dir = os.path.join(tmp, 'output', 'virus')
    os.makedirs(os.path.join(tmp, 'output', 'logs', 'virus'), exist_ok=True)
    stub = types.SimpleNamespace(
        params=types.SimpleNamespace(
            contig_samples=[sample], assemblers=[assembler],
            pairs=[(assembler, sample)], predictors=predictors,
            row_sample={sample: sample},
            row_stem={sample: stem or (sample + '.contigs')},
            long_read_rows=[],
            genomad_base=os.path.join(out_dir, 'genomad'),
            checkv_base=os.path.join(out_dir, 'checkv'),
            virus_base=out_dir),
        output=types.SimpleNamespace(
            summary=os.path.join(out_dir, 'virus_summary.tsv'),
            mapping=os.path.join(out_dir, 'virus_name_mapping.tsv')),
        log=[os.path.join(tmp, 'output', 'logs', 'virus', 'summary.log')])
    saved_out, saved_err, saved_cwd = sys.stdout, sys.stderr, os.getcwd()
    os.chdir(tmp)
    try:
        with open(SCRIPT) as handle:
            code = compile(handle.read(), SCRIPT, 'exec')
        exec(code, {'snakemake': stub, '__name__': '__main__'})
    finally:
        sys.stdout, sys.stderr = saved_out, saved_err
        os.chdir(saved_cwd)
    return (read_tsv(stub.output.summary), read_tsv(stub.output.mapping),
            open(stub.log[0]).read())


def read_tsv(path):
    with open(path) as handle:
        lines = handle.read().splitlines()
    columns = lines[0].split('\t')
    return [dict(zip(columns, line.split('\t'))) for line in lines[1:]]


def test_score_calibration():
    """geNomad computes FDR only with calibrated scores, and calibration also
    turns on its --max-fdr filter, so the flag is a switch and the run records
    which way it went."""
    print('geNomad score calibration:')
    _, out = plan('-p')
    check('on by default, the flag reaches the geNomad command',
          '--enable-score-calibration' in out)
    _, out = plan('-p', overrides='params:\n  genomad:\n    score_calibration: false\n')
    check('switched off, the flag is absent',
          'genomad end-to-end' in out
          and '--enable-score-calibration' not in out)


def test_summary():
    print('summary and name mapping:')
    with tempfile.TemporaryDirectory(prefix='magmaker_vsum_') as tmp:
        stem = build_tree(tmp)
        every = ['genomad', 'vibrant', 'cenotetaker3', 'virsorter2']
        rows, mapping, log = run_summary(tmp, every, stem=stem)

        by_predictor = {}
        for row in rows:
            by_predictor.setdefault(row['Predictor'], []).append(row)
        check('every predictor has rows, labelled with its own name',
              sorted(by_predictor) == sorted(every))

        # Each tool's mapping rule, against the header shapes it really writes.
        expected = {
            'genomad': {source for _, source in GENOMAD_VIRUS} | {'k141_555'},
            'vibrant': {source for _, source in VIBRANT_PHAGES},
            'virsorter2': {source for _, source in VIRSORTER2_SEQS},
            'cenotetaker3': {source for _, _, source in CT3_SEQS},
        }
        for predictor in every:
            got = {row['Source_Contig'] for row in by_predictor[predictor]}
            check('{0}: names map to their source contigs'.format(predictor),
                  got == expected[predictor])

        check('VIBRANT: _fragment_N after len= does not reach the contig name',
              all('fragment' not in row['Source_Contig']
                  for row in by_predictor['vibrant']))
        check('a __dup<N> suffix is stripped',
              all('__dup' not in row['Source_Contig'] for row in rows))
        check("geNomad: a provirus' coordinates are cut off the name",
              'k141_99' in {row['Source_Contig']
                            for row in by_predictor['genomad']})

        # The mapping report is what makes a wrong rule visible.
        check('every predictor maps every sequence, with nothing unmapped',
              all(row['Unmapped'] == '0' for row in mapping))
        check('no mapped name is missing from the assembly',
              all(row['Not_In_Assembly'] == '0' for row in mapping))
        check('the report counts the sequences each predictor wrote',
              {row['Predictor']: row['Sequences'] for row in mapping}
              == {'genomad': '3', 'vibrant': '3', 'cenotetaker3': '2',
                  'virsorter2': '3'})

        # The five provirus columns.
        provirus = [row for row in by_predictor['genomad']
                    if row['Source_Contig'] == 'k141_99'][0]
        check("geNomad's excised provirus: predictor Yes with its coordinates "
              "in the source contig, CheckV No, Any Yes",
              provirus['Provirus_Predictor'] == 'Yes'
              and provirus['Provirus_Predictor_Coords'] == '120-4100'
              and provirus['provirus'] == 'No'
              and provirus['Provirus_Any'] == 'Yes')
        plain = [row for row in by_predictor['genomad']
                 if row['Source_Contig'] == 'k141_7410'][0]
        check('a virus that is not a provirus: all three No',
              plain['Provirus_Predictor'] == 'No'
              and plain['provirus'] == 'No'
              and plain['Provirus_Any'] == 'No')

        # Plasmids: geNomad only, from the same run.
        plasmids = [row for row in rows if row['Element_Type'] == 'plasmid']
        check('plasmids come from geNomad alone, with its gene columns',
              len(plasmids) == 1 and plasmids[0]['Predictor'] == 'genomad'
              and plasmids[0]['Conjugation_Genes'] == 'traD'
              and plasmids[0]['AMR_Genes'] == 'tetA')

        # geNomad's own columns survive, and CheckV is joined on.
        check("geNomad's score, topology and taxonomy are carried",
              plain['Score'] == '0.95' and plain['Topology'] == 'No terminal repeats'
              and plain['Taxonomy'].startswith('Viruses')
              and plain['checkv_quality'] == 'High-quality')

        # With only geNomad on, the table is the old one plus the new columns.
        genomad_only, mapping_only, _ = run_summary(tmp, ['genomad'], stem=stem)
        check('geNomad alone gives exactly its own rows',
              len(genomad_only) == len(by_predictor['genomad'])
              and {row['Predictor'] for row in genomad_only} == {'genomad'})
        check('geNomad alone reports only itself in the mapping table',
              {row['Predictor'] for row in mapping_only} == {'genomad'})


def test_missing_output_is_not_zero():
    """A predictor that wrote no file is reported as absent, not as 0 calls."""
    print('a missing predictor output:')
    with tempfile.TemporaryDirectory(prefix='magmaker_vmiss_') as tmp:
        stem = build_tree(tmp)
        fasta = os.path.join(tmp, 'output', 'virus', 'vibrant', 'megahit', 'S1',
                             'VIBRANT_' + stem, 'VIBRANT_phages_' + stem,
                             stem + '.phages_combined.fna')
        os.remove(fasta)
        rows, mapping, log = run_summary(tmp, ['genomad', 'vibrant'], stem=stem)
        vibrant = [row for row in mapping if row['Predictor'] == 'vibrant'][0]
        check('the mapping table says the output was not present',
              vibrant['Output_Present'] == 'no' and vibrant['Sequences'] == '0')
        check('the log says so too, naming the file',
              'has no output file' in log)
        check('geNomad is unaffected',
              any(row['Predictor'] == 'genomad' for row in rows))


def test_wrong_rule_would_be_caught():
    """If a mapping rule were wrong, the report must say so rather than let a
    small overlap look like genuine disagreement. Modelled by a name that is
    not in the assembly."""
    print('a name that is not in the assembly:')
    with tempfile.TemporaryDirectory(prefix='magmaker_vbad_') as tmp:
        stem = build_tree(tmp)
        write_fasta(os.path.join(tmp, 'output', 'virus', 'virsorter2', 'megahit',
                                 'S1', 'final-viral-combined.fa'),
                    ['k141_NOTREAL||full'])
        rows, mapping, log = run_summary(tmp, ['virsorter2'], stem=stem)
        entry = mapping[0]
        check('the count of names absent from the assembly is reported',
              entry['Not_In_Assembly'] == '1'
              and entry['Example_Not_In_Assembly'] == 'k141_NOTREAL')
        check('the log warns that the rule is probably wrong',
              'probably wrong' in log)


if __name__ == '__main__':
    default = test_plan()
    test_score_calibration()
    test_summary()
    test_missing_output_is_not_zero()
    test_wrong_rule_would_be_caught()
    print('test_virus_track: {0} checks passed'.format(len(CHECKS)))
