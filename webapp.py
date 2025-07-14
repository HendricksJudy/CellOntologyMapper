"""Minimal Flask app to browse example notebooks."""

from __future__ import annotations

from pathlib import Path

from flask import Flask, render_template_string, abort, send_file
import nbformat
from nbconvert import HTMLExporter

app = Flask(__name__)

# Directory containing the Jupyter notebooks bundled with this repository
TUTORIAL_DIR = Path(__file__).parent / "tutorials"

INDEX_HTML = """
<h2>Available Notebooks</h2>
<ul>
{% for nb in notebooks %}
  <li><a href="{{ url_for('view_notebook', name=nb.name) }}">{{ nb.name }}</a> (<a href="{{ url_for('download_notebook', name=nb.name) }}">download</a>)</li>
{% endfor %}
</ul>
"""

@app.route('/')
def index():
    """List notebooks available in the tutorials directory."""
    notebooks = sorted(TUTORIAL_DIR.glob('*.ipynb'))
    return render_template_string(INDEX_HTML, notebooks=notebooks)

@app.route('/notebook/<path:name>')
def view_notebook(name: str):
    """Render the selected notebook as HTML."""
    path = TUTORIAL_DIR / name
    if not path.exists() or path.suffix != '.ipynb':
        abort(404)
    nb = nbformat.read(path, as_version=4)
    html_exporter = HTMLExporter()
    body, _ = html_exporter.from_notebook_node(nb)
    return body

@app.route('/download/<path:name>')
def download_notebook(name: str):
    """Allow users to download the raw notebook."""
    path = TUTORIAL_DIR / name
    if not path.exists() or path.suffix != '.ipynb':
        abort(404)
    return send_file(path, as_attachment=True)

if __name__ == '__main__':
    app.run(debug=True)
