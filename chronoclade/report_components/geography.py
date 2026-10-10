"""Country counts shared by dense and large profile reports."""
from collections import Counter
from pathlib import Path
import re


def write_country_figure(output, cohort_id, records):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    if not re.fullmatch(r'[A-Za-z0-9_-]+', cohort_id):
        raise ValueError('Country figure requires a safe cohort identifier')
    counts = Counter(str(row.get('country') or 'Unknown') for row in records)
    fig, ax = plt.subplots(figsize=(7, max(3, len(counts) * .3)))
    keys = sorted(counts, key=lambda key: (-counts[key], key))
    ax.barh(keys, [counts[key] for key in keys])
    ax.invert_yaxis()
    ax.set(xlabel='Genomes in analysed cohort', title='Available country metadata')
    fig.tight_layout()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    path = output / f'{cohort_id}_countries.svg'
    fig.savefig(path)
    plt.close(fig)
    return str(path)
