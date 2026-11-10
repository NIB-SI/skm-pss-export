# PSS exports (reference)

The exports of the latest PSS made by pss-export, kept here to check how the exports run: a change to pss-export or
to PSS shows as a diff of these files. To use PSS, download the exports from https://skm.nib.si (every species,
public and restricted data).

- Public data, *Arabidopsis thaliana* (`ath`); every format of `pss_export.formats` except FAIDARE.
- `<format>/<file>.<ext>`; the models (SBML, TabularQual, BoolNet, SBGN) also as `<file>-with-model-fixes.<ext>`.
- `version.txt`: the date of the PSS dump, the pss-export version that made them, and the export date.
- The files keep their names and are overwritten at each update; older versions are in the git history.

Remake them (the database connection from the environment or `.env`, as for `PSSAdapter`):

    python scripts/export_pss.py <dump date, YYYY-MM-DD>

License: CC BY-NC 4.0 (see [LICENSE](LICENSE)); commercial use by individual agreement. The code is MIT.
