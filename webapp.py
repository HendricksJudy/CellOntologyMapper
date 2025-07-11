"""Simple Flask web app for running the trophoblast pipeline"""

from __future__ import annotations

import tempfile
from pathlib import Path

from flask import (
    Flask,
    flash,
    redirect,
    render_template_string,
    request,
    session,
    url_for,
)
from werkzeug.utils import secure_filename

import colormaps as cmaps

from interactive_trophoblast_pipeline import (
    load_adata,
    clean_celltype_predictions,
    palette_from_categories,
    run_mapper,
)


app = Flask(__name__)
app.secret_key = "change-me"
UPLOAD_DIR = Path(tempfile.gettempdir()) / "com_webapp"
UPLOAD_DIR.mkdir(parents=True, exist_ok=True)


INDEX_HTML = """
<h2>CellOntologyMapper Web Interface</h2>
<p><a href='{{ url_for("upload_h5ad") }}'>Start mapping pipeline</a></p>
"""

UPLOAD_H5AD_HTML = """
<h2>Step 1: Upload h5ad file</h2>
<form method=post enctype=multipart/form-data>
  <input type=file name=h5ad required>
  <input type=submit value='Upload'>
</form>
"""

UPLOAD_RESOURCES_HTML = """
<h2>Step 2: Upload ontology resources</h2>
<form method=post enctype=multipart/form-data>
  <p>Cell Ontology JSON: <input type=file name=cl_json required></p>
  <p>Cell Taxonomy TXT: <input type=file name=taxonomy required></p>
  <p>Model cache directory: <input type=text name=model_dir value='model_cache'></p>
  <input type=submit value='Run pipeline'>
</form>
"""

RESULT_HTML = """
<h2>Pipeline completed</h2>
<p>Annotated metadata saved to: {{ csv_path }}</p>
"""


def save_upload(upload_file):
    filename = secure_filename(upload_file.filename)
    path = UPLOAD_DIR / filename
    upload_file.save(path)
    return path


@app.route("/")
def index():
    return render_template_string(INDEX_HTML)


@app.route("/upload_h5ad", methods=["GET", "POST"])
def upload_h5ad():
    if request.method == "POST":
        file = request.files.get("h5ad")
        if file:
            path = save_upload(file)
            session["h5ad"] = str(path)
            return redirect(url_for("upload_resources"))
        flash("File required")
    return render_template_string(UPLOAD_H5AD_HTML)


@app.route("/upload_resources", methods=["GET", "POST"])
def upload_resources():
    if "h5ad" not in session:
        return redirect(url_for("upload_h5ad"))
    if request.method == "POST":
        cl_json = request.files.get("cl_json")
        taxonomy = request.files.get("taxonomy")
        model_dir = request.form.get("model_dir", "model_cache")
        if cl_json and taxonomy:
            cl_json_path = save_upload(cl_json)
            taxonomy_path = save_upload(taxonomy)
            session["cl_json"] = str(cl_json_path)
            session["taxonomy"] = str(taxonomy_path)
            session["model_dir"] = model_dir
            return redirect(url_for("run_pipeline_route"))
        flash("All files required")
    return render_template_string(UPLOAD_RESOURCES_HTML)


@app.route("/run_pipeline")
def run_pipeline_route():
    if not {"h5ad", "cl_json", "taxonomy"} <= session.keys():
        return redirect(url_for("upload_h5ad"))
    h5ad = Path(session["h5ad"])
    cl_json = Path(session["cl_json"])
    taxonomy = Path(session["taxonomy"])
    model_dir = Path(session.get("model_dir", "model_cache"))

    adata = load_adata(h5ad)
    clean_celltype_predictions(adata)
    cats = list(adata.obs["celltype_predictions"].cat.categories)
    adata.uns["celltype_predictions_colors"] = palette_from_categories(
        cats, cmaps.cet_g_bw_minc_minl
    )
    adata = run_mapper(adata, cl_json, taxonomy, model_dir)

    out_csv = h5ad.with_name(h5ad.stem + "_annotated.csv")
    adata.obs.to_csv(out_csv)
    csv_path = str(out_csv.resolve())
    return render_template_string(RESULT_HTML, csv_path=csv_path)


if __name__ == "__main__":
    app.run(debug=True)
