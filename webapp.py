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
    Response,
    send_from_directory,
    stream_with_context,
)
from threading import Thread
from queue import Queue
import contextlib
import io
import logging
import json
from werkzeug.utils import secure_filename

import colormaps as cmaps

from interactive_trophoblast_pipeline import (
    load_adata,
    clean_celltype_predictions,
    palette_from_categories,
    run_mapper,
    plot_umap,
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

PIPELINE_HTML = """
<h2>Running pipeline...</h2>
<pre id='log' style='height:300px; overflow:auto; background:#f0f0f0; padding:5px;'></pre>
<div id='results'></div>
<script>
const logElem = document.getElementById('log');
const results = document.getElementById('results');
const es = new EventSource('{{ url_for("events") }}');
es.onmessage = (e) => {
  const data = JSON.parse(e.data);
  if (data.type === 'log') {
    logElem.textContent += data.message + "\n";
    logElem.scrollTop = logElem.scrollHeight;
  } else if (data.type === 'csv') {
    const a = document.createElement('a');
    a.href = '/files/' + data.filename;
    a.textContent = 'Download annotated CSV';
    results.appendChild(a);
  } else if (data.type === 'figure') {
    const img = document.createElement('img');
    img.src = '/files/' + data.filename;
    img.style.width = '300px';
    results.appendChild(img);
  } else if (data.type === 'done') {
    es.close();
  }
};
</script>
"""

event_queue: Queue[str] = Queue()
pipeline_started = False

class QueueWriter(io.TextIOBase):
    def __init__(self, queue: Queue[str]):
        self.queue = queue
        self.buf = ""

    def write(self, msg: str) -> int:
        self.buf += msg
        while "\n" in self.buf:
            line, self.buf = self.buf.split("\n", 1)
            line = line.strip()
            if line:
                self.queue.put_nowait(json.dumps({"type": "log", "message": line}))
        return len(msg)

    def flush(self) -> None:
        if self.buf.strip():
            self.queue.put_nowait(json.dumps({"type": "log", "message": self.buf.strip()}))
        self.buf = ""

def pipeline_worker(h5ad: Path, cl_json: Path, taxonomy: Path, model_dir: Path) -> None:
    figdir = UPLOAD_DIR
    writer = QueueWriter(event_queue)
    handler = logging.StreamHandler(writer)
    root_logger = logging.getLogger()
    root_logger.addHandler(handler)
    figs: list[Path] = []
    out_csv = h5ad.with_name(h5ad.stem + "_annotated.csv")
    try:
        with contextlib.redirect_stdout(writer):
            adata = load_adata(h5ad)
            clean_celltype_predictions(adata)
            cats = list(adata.obs["celltype_predictions"].cat.categories)
            adata.uns["celltype_predictions_colors"] = palette_from_categories(cats, cmaps.cet_g_bw_minc_minl)
            plot_umap(adata, "celltype_predictions", "Trophoblast (author labels)", adata.uns["celltype_predictions_colors"], figdir / "umap_trophoblast_author.png", size=6)
            adata = run_mapper(adata, cl_json, taxonomy, model_dir)
            map_targets = {
                "enhanced_cell_ontology": "Cell Ontology term",
                "enhanced_cell_ontology_taxonomy_match": "Cell Taxonomy term",
                "enhanced_cell_ontology_ct_id": "Cell Taxonomy ID",
                "enhanced_cell_ontology_cl_id": "Cell Ontology ID",
            }
            figs = [figdir / "umap_trophoblast_author.png"]
            for key, label in map_targets.items():
                mapping = (
                    adata.obs[["celltype_predictions", key]]
                    .drop_duplicates()
                    .set_index("celltype_predictions")[key]
                    .to_dict()
                )
                palette = {
                    mapping[k]: adata.uns["celltype_predictions_colors"][i]
                    for i, k in enumerate(mapping.keys())
                    if k in mapping
                }
                figpath = figdir / f"umap_trophoblast_{key}.png"
                plot_umap(adata, key, f"Trophoblast – {label}", palette, figpath, size=6)
                figs.append(figpath)
            adata.obs.to_csv(out_csv)
    finally:
        writer.flush()
        root_logger.removeHandler(handler)
    event_queue.put(json.dumps({"type": "csv", "filename": out_csv.name}))
    for f in figs:
        event_queue.put(json.dumps({"type": "figure", "filename": f.name}))
    event_queue.put(json.dumps({"type": "done"}))



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
    return render_template_string(PIPELINE_HTML)


@app.route("/events")
def events():
    global pipeline_started
    if not pipeline_started:
        pipeline_started = True
        h5ad = Path(session["h5ad"])
        cl_json = Path(session["cl_json"])
        taxonomy = Path(session["taxonomy"])
        model_dir = Path(session.get("model_dir", "model_cache"))
        Thread(target=pipeline_worker, args=(h5ad, cl_json, taxonomy, model_dir)).start()

    def stream():
        global pipeline_started
        while True:
            msg = event_queue.get()
            yield f"data: {msg}\n\n"
            if msg == json.dumps({"type": "done"}):
                break
        pipeline_started = False

    return Response(stream_with_context(stream()), mimetype="text/event-stream")


@app.route("/files/<path:filename>")
def files(filename: str):
    return send_from_directory(UPLOAD_DIR, filename)


if __name__ == "__main__":
    app.run(debug=True)
