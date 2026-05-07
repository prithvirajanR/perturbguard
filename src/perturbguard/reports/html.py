from __future__ import annotations

from pathlib import Path

import pandas as pd
from jinja2 import Template


TEMPLATE = """
<!doctype html>
<html lang="en">
<head><meta charset="utf-8"><title>PerturbGuard Report</title>
<style>
body{font-family:Arial,sans-serif;margin:0;line-height:1.45;color:#202124;background:#f8fafc}
header{background:#ffffff;border-bottom:1px solid #d9e2ec;padding:1.25rem 2rem;position:sticky;top:0;z-index:2}
main{padding:1.5rem 2rem}.toolbar{display:flex;gap:.75rem;flex-wrap:wrap;align-items:center}
input,select{border:1px solid #b8c2cc;border-radius:6px;padding:.5rem .65rem;background:white}
.summary{display:flex;gap:.75rem;flex-wrap:wrap;margin:1rem 0}.summary-card{background:#fff;border:1px solid #d9e2ec;border-radius:8px;padding:.75rem 1rem;min-width:7rem}
section{background:#fff;border:1px solid #d9e2ec;border-radius:8px;margin:1rem 0;padding:1rem;overflow:auto}
table{border-collapse:collapse;margin:1rem 0;width:100%;font-size:.92rem}td,th{border:1px solid #e3e8ef;padding:.45rem;text-align:left;vertical-align:top}th{background:#f1f5f9;position:sticky;top:4.8rem}
.FAIL{color:#b00020;font-weight:700}.WARNING{color:#8a5a00;font-weight:700}.PASS,.SUPPORTED{color:#0a6b2b;font-weight:700}.UNSUPPORTED{color:#8b1a1a;font-weight:700}.INSUFFICIENT_METADATA{color:#5f6368;font-weight:700}
a{color:#0b57d0}
</style></head>
<body>
<header>
<h1>PerturbGuard Report</h1>
<div class="toolbar">
  <label>Search tables <input data-table-search type="search" placeholder="status, column, message"></label>
  <label>Status <select data-status-filter><option value="">All</option>{% for status in statuses %}<option value="{{ status }}">{{ status }}</option>{% endfor %}</select></label>
</div>
<div class="summary">
{% for status, count in status_counts %}
  <div class="summary-card"><strong class="{{ status }}">{{ status }}</strong><br>{{ count }}</div>
{% endfor %}
</div>
</header>
<main>
{% if plots %}
<h2>Plots</h2>
<ul>
{% for name, path in plots %}
  <li><a href="{{ path }}">{{ name }}</a></li>
{% endfor %}
</ul>
{% endif %}
{% for title, table in sections %}
<section>
<h2>{{ title }}</h2>
{{ table }}
</section>
{% endfor %}
<script>
const search = document.querySelector('[data-table-search]');
const statusFilter = document.querySelector('[data-status-filter]');
function applyFilters(){
  const q = (search.value || '').toLowerCase();
  const status = statusFilter.value;
  document.querySelectorAll('tbody tr').forEach(row => {
    const text = row.innerText.toLowerCase();
    const hasStatus = !status || row.innerText.includes(status);
    row.style.display = text.includes(q) && hasStatus ? '' : 'none';
  });
}
search.addEventListener('input', applyFilters);
statusFilter.addEventListener('change', applyFilters);
</script>
</main></body></html>
"""


def _title(key: str) -> str:
    return key.replace("_", " ").title()


def write_html_report(
    path: str | Path,
    sections: dict[str, pd.DataFrame],
    plot_paths: dict[str, Path] | None = None,
) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rendered_sections = []
    status_counts: dict[str, int] = {}
    for key, df in sections.items():
        if "status" in df:
            for status, count in df["status"].astype(str).value_counts().items():
                status_counts[status] = status_counts.get(status, 0) + int(count)
        rendered_sections.append(
            (_title(key), df.to_html(index=False, escape=True, classes="table", na_rep=""))
        )
    plots = []
    for name, plot_path in (plot_paths or {}).items():
        try:
            rel = plot_path.relative_to(path.parent).as_posix()
        except ValueError:
            rel = plot_path.as_posix()
        plots.append((_title(name), rel))
    ordered_counts = sorted(status_counts.items())
    path.write_text(
        Template(TEMPLATE).render(
            sections=rendered_sections,
            plots=plots,
            statuses=[status for status, _ in ordered_counts],
            status_counts=ordered_counts,
        ),
        encoding="utf-8",
    )
    return path
