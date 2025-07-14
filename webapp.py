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
import json
import subprocess
import sys
from werkzeug.utils import secure_filename


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
  } else if (data.type === 'error') {
    const div = document.createElement('div');
    div.style.color = 'red';
    div.textContent = data.message;
    results.appendChild(div);
  } else if (data.type === 'done') {
    es.close();
  }
};
</script>
"""

event_queue: Queue[str] = Queue()
pipeline_started = False


def pipeline_worker(h5ad: Path, cl_json: Path, taxonomy: Path, model_dir: Path) -> None:
    cmd = [
        sys.executable,
        "monitor_pipeline.py",
        str(h5ad),
        str(cl_json),
        str(taxonomy),
        str(model_dir),
        "--figdir",
        str(UPLOAD_DIR),
    ]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
    for line in process.stdout:
        event_queue.put(json.dumps({"type": "log", "message": line.rstrip()}))
    returncode = process.wait()
    if returncode != 0:
        event_queue.put(json.dumps({"type": "error", "message": f'Pipeline exited with code {returncode}'}))
    out_csv = h5ad.with_name(h5ad.stem + "_annotated.csv")
    event_queue.put(json.dumps({"type": "csv", "filename": out_csv.name}))
    figs = ["umap_trophoblast_author.png"]
    for key in [
        "enhanced_cell_ontology",
        "enhanced_cell_ontology_taxonomy_match",
        "enhanced_cell_ontology_ct_id",
        "enhanced_cell_ontology_cl_id",
    ]:
        figs.append(f"umap_trophoblast_{key}.png")
    for fig in figs:
        event_queue.put(json.dumps({"type": "figure", "filename": fig}))
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
