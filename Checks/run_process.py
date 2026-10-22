#!/usr/bin/env python3
"""Run every step of one process headlessly, each step in a fresh kernel.

    python Checks/run_process.py Processes/LeptonQuark --route FORM --set spin=Polarized
    python Checks/run_process.py Processes/LeptonQuark --steps 4,5 --plan

By default the steps run in place and write into <process>/Results/.  With
--copy DIR a source-only copy of the repository (no Results/, no Paths.wl
unless --paths is given) is made first and the run happens there, which is the
cleanest test that the shipped sources are complete.  With --plan the steps
that would run are listed (with the receipt each must write) and nothing runs.

The process folder is given relative to the repository root (an absolute path
is accepted only if it lies inside this repository; with --copy the same
folder of the copy is run, never the original).  The run stops with a non-zero exit code when the
folder does not exist, when no step file matches the selection, or when a step
does not end with a completed receipt: run_step.wls writes
<process>/Results/<step>[-<override values>].receipt.json after every cell and
sets "completed": true after the last one; the receipt must have been written
by this run (not left over from an earlier one) and name this step file.

Prerequisites of other processes are run first when their results are
missing: Processes/Hadronic reads the Born table of Processes/PartonicBorn.
The FORM and Tracer routes of step 2 compare their squares with the FeynCalc
ones at their end, so run the FeynCalc route of a class first (without it the
comparison is skipped and reported as such).
"""
import argparse, json, os, shutil, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
# process -> [(prerequisite process, file it must have written, its steps to run)]
PREREQUISITES = {'Hadronic': [('PartonicBorn', 'Results/Born.wl', ['1'])]}

def fail(message):
    print(f'run_process: {message}', file=sys.stderr)
    return 2

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('process', help='process folder relative to the repository root, e.g. Processes/LeptonQuark')
    ap.add_argument('--route', choices=['FeynCalc', 'FORM', 'Tracer'], default='FeynCalc',
                    help='algebra route for the 2-Algebra step and the route variable of later steps')
    ap.add_argument('--set', action='append', default=[], metavar='NAME=VALUE',
                    help='override a top-level string choice in the notebooks, e.g. spin=Polarized')
    ap.add_argument('--steps', help='comma-separated step numbers to run, e.g. 1,2 (default: all)')
    ap.add_argument('--copy', type=Path, help='run in a fresh source-only copy at DIR, outside the repository')
    ap.add_argument('--paths', type=Path, help='Paths.wl to use in the copy')
    ap.add_argument('--plan', action='store_true', help='list the steps that would run, and run nothing')
    ap.add_argument('--wolfram', default=shutil.which('wolframscript'))
    a = ap.parse_args()
    for item in a.set:
        if '=' not in item or not item.split('=', 1)[0]:
            ap.error(f'--set {item!r}: expected NAME=VALUE')

    # Copying into the source tree would recursively copy the growing copy itself.
    # resolve() also handles a destination whose existing parent is a symlink.
    if a.copy:
        a.copy = a.copy.resolve()
        if a.copy == ROOT or ROOT in a.copy.parents:
            ap.error('--copy destination must be outside the source repository')

    # the process folder must exist inside the repository (before any copy is made)
    source = (ROOT / a.process).resolve()
    try:
        relative = source.relative_to(ROOT)
    except ValueError:
        return fail(f'{a.process} is not inside the repository {ROOT}; give the process folder relative to it')
    if not source.is_dir() or relative.parts[:1] != ('Processes',) or len(relative.parts) != 2:
        return fail(f'{a.process} is not a process folder (expected Processes/<name>, an existing folder under {ROOT})')

    # the steps: numbered .wl files, the 2-Algebra step only for the chosen route
    wanted = None if not a.steps else [n.strip() for n in a.steps.split(',') if n.strip()]
    if wanted is not None and not wanted:
        return fail('--steps selects no step')
    candidates = sorted(source.glob('[0-9]*.wl'))
    chosen = [p for p in candidates
              if (wanted is None or p.name.split('-')[0] in wanted)
              and ('-Algebra-' not in p.name or p.name.endswith(f'-Algebra-{a.route}.wl'))]
    if wanted is not None:
        missing = [n for n in wanted if not any(p.name.split('-')[0] == n for p in chosen)]
        if missing:
            return fail(f'no step file numbered {", ".join(missing)} in {relative} for route {a.route}')
    if not chosen:
        return fail(f'no step files selected in {relative}')

    overrides = [f'route={a.route}'] + a.set
    suffix = '-' + '-'.join(o.split('=', 1)[1] for o in overrides)
    if a.plan:
        print(f'{relative}  FEYNFACET_SET={";".join(overrides)}' + (f'  (in a copy at {a.copy})' if a.copy else ''))
        for p in chosen:
            print(f'  {p.name:32s} receipt Results/{p.stem}{suffix}.receipt.json')
        return 0

    root = ROOT
    if a.copy:
        a.copy = a.copy.resolve()
        if a.copy.exists():
            ap.error(f'{a.copy} already exists; choose a new directory')
        shutil.copytree(ROOT, a.copy, ignore=shutil.ignore_patterns('.git', 'Results', '__pycache__', 'Paths.wl', '*.log'))
        if a.paths:
            shutil.copyfile(a.paths, a.copy / 'Paths.wl')
        root = a.copy
    if not a.wolfram:
        return fail('wolframscript not found; give it with --wolfram')
    for name, needed, steps in PREREQUISITES.get(relative.name, []):
        if not (root / 'Processes' / name / needed).is_file():
            print(f'Processes/{name}/{needed} is missing: running Processes/{name} (steps {",".join(steps)}) first', flush=True)
            pre = root / 'Processes' / name
            pre_steps = [q for q in sorted(pre.glob('[0-9]*.wl')) if q.name.split('-')[0] in steps]
            status = run_steps(root, pre, [q.name for q in pre_steps], ['route=FeynCalc'], a.wolfram)
            if status:
                return status
    if a.route != 'FeynCalc' and any('-Algebra-' in p.name for p in chosen):
        print(f'note: step 2 of the {a.route} route compares with the FeynCalc squares of the same class when they exist', flush=True)
    return run_steps(root, root / relative, [p.name for p in chosen], overrides, a.wolfram)

def run_steps(root, proc, names, overrides, wolfram):
    """Run the step files `names` of the process folder `proc` in order, each in a fresh kernel;
    0 when every step ends with a completed receipt written by this run, 1 otherwise."""
    suffix = '-' + '-'.join(o.split('=', 1)[1] for o in overrides)
    env = dict(os.environ)
    env['FEYNFACET_SET'] = ';'.join(overrides)
    for p in (proc / n for n in names):
        log = proc / 'Results' / (p.stem + suffix + '.log')
        receipt = proc / 'Results' / f'{p.stem}{suffix}.receipt.json'
        log.parent.mkdir(exist_ok=True)
        t0 = time.time()
        with log.open('w') as stream:
            r = subprocess.run([wolfram, '-file', str(root / 'Checks' / 'run_step.wls'), str(p)],
                               cwd=proc, env=env, stdout=stream, stderr=subprocess.STDOUT)
        print(f'{p.name:32s} exit {r.returncode}  {time.time() - t0:8.1f} s', flush=True)
        problem = None
        if r.returncode != 0:
            problem = f'exit code {r.returncode}'
        elif not receipt.is_file():
            problem = f'no receipt {receipt}'
        elif receipt.stat().st_mtime < t0 - 1:
            problem = f'the receipt {receipt} was not written by this run'
        else:
            try:
                data = json.loads(receipt.read_text())
            except ValueError:
                data = {}
            if data.get('completed') is not True:
                problem = f'the receipt {receipt} does not say "completed": true'
            elif Path(str(data.get('file', ''))).resolve() != p.resolve():
                problem = f'the receipt {receipt} belongs to {data.get("file")}'
        if problem:
            print(f'stopped: {p.name}: {problem}; see {log}', file=sys.stderr)
            print(log.read_text()[-3000:], file=sys.stderr)
            return 1
    print(f'all {len(names)} steps of {proc.name} completed, each with a completed receipt')
    return 0

if __name__ == '__main__':
    sys.exit(main())
