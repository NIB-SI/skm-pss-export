#!/usr/bin/env python3
'''
Make the exports of the latest PSS kept in this repository (pss/): every format of the registry
(pss_export.formats) except those in EXCLUDE, public data, Arabidopsis thaliana; the models also with model fixes.

The files have fixed names and are overwritten at each update; older versions are in the git history.

    python scripts/export_pss.py 2026-11-09              # the date of the PSS dump in the database
    python scripts/export_pss.py 2026-11-09 --out pss

The database connection comes from the environment or a .env file (MY_NEO4J_URI, ...), as for PSSAdapter.
'''

import argparse
import datetime
import importlib.metadata
import shutil
import subprocess
from pathlib import Path

import pss_export
from pss_export import PSSAdapter
from pss_export.formats import FORMATS

ACCESS = 'public'
SPECIES = 'ath'
EXCLUDE = {'faidare'}       # every species, large and not of interest here
MODEL_FIXES_SUFFIX = '-with-model-fixes'


def pss_export_version():
    ''' The version of the pss_export package used: from git when it is a checkout (e.g. v0.2.1, or
    v0.2.1-3-gabc1234 after the tag; the metadata of an editable install can be out of date), else the installed
    version '''
    try:
        return subprocess.run(['git', 'describe', '--tags', '--always'], cwd=Path(pss_export.__file__).parent,
                              capture_output=True, text=True, check=True).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return importlib.metadata.version('pss-export')


def export_all(dump_date, out):
    formats = [fmt for key, fmt in FORMATS.items() if key not in EXCLUDE]

    # remove the format folders (only those), so files of an excluded format or a dropped file do not stay behind
    for key in FORMATS:
        if (out / key).is_dir():
            shutil.rmtree(out / key)
    out.mkdir(parents=True, exist_ok=True)

    adapter = PSSAdapter(model_id='pss',
                         model_name='Plant Stress Signalling (PSS)',
                         model_version=dump_date)
    adapter.collect_reactions(access=ACCESS, species=SPECIES)

    written = []
    for fmt in formats:
        folder = out / fmt.key
        folder.mkdir()
        variants = [(False, '')] + ([(True, MODEL_FIXES_SUFFIX)] if fmt.model_fixes else [])
        for model_fixes, suffix in variants:
            files = {f.argument: folder / f'{f.name}{suffix}.{f.extension}' for f in fmt.files}
            adapter.export(fmt.key, model_fixes=model_fixes, **{a: str(p) for a, p in files.items()})
            written.extend(files.values())

    (out / 'version.txt').write_text(
        f'PSS dump\t{dump_date}\n'
        f'pss-export\t{pss_export_version()}\n'
        f'exported\t{datetime.date.today().isoformat()}\n'
        f'access\t{ACCESS}\n'
        f'species\t{SPECIES}\n')
    return written


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('dump_date', type=datetime.date.fromisoformat,
                        help='date of the PSS dump in the database (YYYY-MM-DD)')
    parser.add_argument('--out', type=Path, default=Path('pss'), help='output folder (default: pss)')
    args = parser.parse_args()

    print(f'pss_export from {pss_export.__file__}')
    for path in export_all(args.dump_date.isoformat(), args.out):
        print(f'{path.stat().st_size:>12,}  {path}')


if __name__ == '__main__':
    main()
