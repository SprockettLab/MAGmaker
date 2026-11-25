#!/usr/bin/env python3
"""Every rule shell that runs `python` puts its own conda env first on PATH.

On demon another environment's bin can reach a job ahead of the rule's own
env, so a bare `python` in a rule shell runs the wrong interpreter. It fails
only when the script imports a package that other interpreter lacks: on
WF18_FullTest_2026 (2026-11-24) all 14 polish_long_read_assembly jobs ran
Polypolish and pypolca to completion, then died in polish_changes.py on
`import edlib`, which the rule's polish.yaml env has (python-edlib 1.3.9 on
python 3.13.18). Same failure as run_isescan (PR #53).

For every shell block in resources/snakefiles/*.smk that runs python (a line
starting with `python`, or `| python`), checks that
`export PATH="$CONDA_PREFIX/bin:$PATH"` comes before the first such line,
and that the rule has a conda: env for $CONDA_PREFIX to point at.

    python3 resources/test/test_rule_python_env.py
"""

import glob
import os
import re
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
EXPORT = 'export PATH="$CONDA_PREFIX/bin:$PATH"'
PY_CALL = re.compile(r'^\s*(?:\|\s*)?python3?\s')
FAILURES = []


def rules(text):
    """(rule name, rule text) for every rule in a snakefile."""
    starts = [(m.start(), m.group(1)) for m in
              re.finditer(r'^(?:local)?rule (\w+):', text, re.M)]
    for i, (pos, name) in enumerate(starts):
        end = starts[i + 1][0] if i + 1 < len(starts) else len(text)
        yield name, text[pos:end]


def main():
    checked = 0
    for path in sorted(glob.glob(os.path.join(REPO_ROOT, 'resources/snakefiles/*.smk'))):
        for name, body in rules(open(path).read()):
            m = re.search(r'^\s+shell:\s*\n(.*?)"""(.*?)"""', body, re.M | re.S)
            if not m:
                continue
            lines = m.group(2).split('\n')
            calls = [i for i, l in enumerate(lines) if PY_CALL.match(l)]
            if not calls:
                continue
            checked += 1
            exports = [i for i, l in enumerate(lines) if l.strip() == EXPORT]
            where = '%s: rule %s' % (os.path.basename(path), name)
            ok = bool(exports) and exports[0] < calls[0]
            print(("PASS  " if ok else "FAIL  ") + where + ': env bin first before python')
            if not ok:
                FAILURES.append(where)
            has_conda = re.search(r'^\s+conda:\s*$', body, re.M) is not None
            print(("PASS  " if has_conda else "FAIL  ") + where + ': has a conda env')
            if not has_conda:
                FAILURES.append(where + ' (no conda env)')
    print("\n%d rule(s) run python in their shell" % checked)
    if checked == 0:
        FAILURES.append('no rule found that runs python; the pattern is broken')


if __name__ == '__main__':
    main()
    print("%d failure(s)" % len(FAILURES))
    sys.exit(1 if FAILURES else 0)
