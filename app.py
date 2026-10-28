#!/usr/bin/env python3
"""Streamlit app for msap_sim.py: sequence-only EcoRI-HpaII/MspI MSAP primer screening.

Keep this file next to msap_sim.py, then:
    pip install -r requirements.txt
    streamlit run app.py

.streamlit/config.toml lifts the upload limit to 1 GB. On your own machine the
"Local path" source reads a whole-genome FASTA in place, without uploading it.
Runs on old Streamlit installs (>= 1.11) and on current releases.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import inspect
import io
import itertools
import json
import os
import random
import re
import zipfile
import zlib
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from msap_sim import SITES, candidates

APP_VERSION = "1.0"
HERE = Path(__file__).resolve().parent
DEMO_FILE = HERE / "synthetic_demo.fasta"
MAX_N = 3          # selective-extension length kept per end, as in msap_sim.analyze
SIZE_CAP = 1500    # fragments are generated once up to this length, then filtered
FRAGMENT_FIELDS = ["contig", "start_1based", "end_1based", "bp",
                   "E_extension", "HM_extension", "orientation"]
RANK_FIELDS = ["E_extension", "HM_extension", "fragments", "occupied_size_bins",
               "co_migrating_excess", "unique_contigs", "min_bp", "max_bp"]
RANK_LABELS = {"rank": "#", "E_extension": "E +", "HM_extension": "HM +", "score": "Score",
               "fragments": "Fragments", "occupied_size_bins": "Bins",
               "co_migrating_excess": "Co-migrating", "unique_contigs": "Contigs",
               "min_bp": "Min bp", "max_bp": "Max bp"}
METRICS = {"score": "Rank score", "fragments": "Fragments",
           "occupied_size_bins": "Occupied size bins",
           "co_migrating_excess": "Co-migrating excess"}
SOURCES = ["Demo (synthetic)", "Upload FASTA", "Local path"]
UPLOAD_LIMIT_MB = 1024   # per-widget on new Streamlit; .streamlit/config.toml on old
ORGANELLE_IDS = {"mt", "pt", "chrm", "chrc", "chrmt", "chrpt"}
ORGANELLE_DESC = re.compile(r"mitochondri|chloroplast|plastid", re.IGNORECASE)
FASTA_MEMBER = re.compile(r"\.(fa|fasta|fna|fas)(\.gz)?$", re.IGNORECASE)
GEL_LEVELS = {1: (0.45, 2.0), 2: (0.70, 3.0), 3: (0.88, 4.0), 4: (1.0, 5.0)}  # alpha, px
BAND_RGB = "170,255,205"
LADDER_COLOUR = "rgba(215,215,215,0.9)"
GEL_BACKGROUND = "#101418"
TRACE_COLOUR = "#1f9e89"
PEAK_COLOUR = "#e4572e"


# ------------------------------------------------------------- compatibility ----

def supported(fn, **kwargs):
    """Keep only the keyword arguments this Streamlit version accepts."""
    params = inspect.signature(fn).parameters
    kept = {}
    for name, value in kwargs.items():
        if name in params:
            kept[name] = value
    return kept


def stretch(fn):
    """Full-width keyword: width='stretch' on new Streamlit, use_container_width on old."""
    params = inspect.signature(fn).parameters
    if "width" in params and isinstance(params["width"].default, str):
        return {"width": "stretch"}
    if "use_container_width" in params:
        return {"use_container_width": True}
    return {}


def show_chart(fig, key):
    kwargs = stretch(st.plotly_chart)
    kwargs.update(supported(st.plotly_chart, key=key, config={"displaylogo": False}))
    st.plotly_chart(fig, **kwargs)


def show_table(df, height=None, hide_index=True):
    kwargs = stretch(st.dataframe)
    kwargs.update(supported(st.dataframe, hide_index=hide_index))
    if height is not None:
        kwargs["height"] = int(height)
    st.dataframe(df, **kwargs)


def download(target, label, data, file_name, mime):
    target.download_button(label, data, file_name=file_name, mime=mime,
                           **stretch(st.download_button))


# --------------------------------------------------------------------- input ----

def demo_fasta():
    """synthetic_demo.fasta if present, else regenerated exactly as make_demo.py does."""
    if DEMO_FILE.exists():
        return DEMO_FILE.read_bytes()
    rng = random.Random(20261016)
    parts = []
    for i in range(120):
        parts.append("GAATTC" if i % 2 == 0 else "CCGG")
        parts.append("".join(rng.choices("ACGT", k=rng.randint(130, 310))))
    seq = "".join(parts)
    lines = [">synthetic_demo_not_rice"]
    for i in range(0, len(seq), 80):
        lines.append(seq[i:i + 80])
    return ("\n".join(lines) + "\n").encode()


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        block = handle.read(1 << 20)
        while block:
            digest.update(block)
            block = handle.read(1 << 20)
    return digest.hexdigest()


@contextmanager
def fasta_lines(raw):
    """Text lines from plain, gzip/bgzip or zip input (largest FASTA member of a zip)."""
    head = raw.read(4)
    raw.seek(0)
    opened = []
    stream = raw
    if head == b"PK\x03\x04":
        archive = zipfile.ZipFile(raw)
        opened.append(archive)
        member = None
        for info in archive.infolist():
            if info.is_dir() or not FASTA_MEMBER.search(info.filename):
                continue
            if member is None or info.file_size > member.file_size:
                member = info
        if member is None:
            raise ValueError("no .fa/.fasta/.fna file inside the zip")
        stream = archive.open(member)
        opened.append(stream)
        head = stream.peek(2)[:2]
    if head[:2] == b"\x1f\x8b":
        stream = gzip.GzipFile(fileobj=stream)
        opened.append(stream)
    text = io.TextIOWrapper(stream, encoding="latin-1")
    try:
        yield text
    finally:
        text.detach()
        for item in reversed(opened):
            item.close()


def iter_fasta(lines):
    """Yield (id, description, SEQUENCE); ids parsed as in msap_sim.fasta_records."""
    name = None
    desc = ""
    chunks = []
    unnamed = 0
    for line in lines:
        if line.startswith(">"):
            if name is not None:
                yield name, desc, "".join(chunks).upper()
            head = line[1:].split(None, 1)
            if head:
                name = head[0]
                desc = head[1].strip() if len(head) > 1 else ""
            else:
                unnamed += 1
                name = f"unnamed_{unnamed}"
                desc = ""
            chunks = []
        else:
            chunks.append(line.strip())
    if name is not None:
        yield name, desc, "".join(chunks).upper()


def empty_fragments():
    columns = {}
    for name in FRAGMENT_FIELDS:
        dtype = "int64" if name in ("start_1based", "end_1based", "bp") else "object"
        columns[name] = pd.Series(dtype=dtype)
    return pd.DataFrame(columns)


def scan_genome(raw, size, tick=None):
    """One pass: per-contig stats and every EcoRI-CCGG fragment of 1..SIZE_CAP bp.

    Fragments come from msap_sim.candidates with max_n=3, as in msap_sim.analyze,
    so any size window inside 1..SIZE_CAP is a plain filter of this table.
    """
    rows = []
    frames = []
    with fasta_lines(raw) as lines:
        for name, desc, seq in iter_fasta(lines):
            counts = {}
            for base in "ACGTN":
                counts[base] = seq.count(base)
            rows.append({
                "contig": name,
                "description": desc,
                "length_bp": len(seq),
                "acgt": counts["A"] + counts["C"] + counts["G"] + counts["T"],
                "gc": counts["G"] + counts["C"],
                "n": counts["N"],
                "EcoRI_sites": seq.count(SITES["E"]),
                "CCGG_sites": seq.count(SITES["H"]),
                "organelle": (name.lower() in ORGANELLE_IDS
                              or ORGANELLE_DESC.search(desc) is not None),
            })
            found = list(candidates(name, seq, 1, SIZE_CAP, MAX_N))
            if found:
                frames.append(pd.DataFrame.from_records(found, columns=FRAGMENT_FIELDS))
            if tick is not None:
                tick(min(raw.tell() / max(size, 1), 1.0), name)
    if not rows:
        raise ValueError("no FASTA records found")
    contigs = pd.DataFrame(rows)
    frags = pd.concat(frames, ignore_index=True) if frames else empty_fragments()
    frags["contig"] = pd.Categorical(frags["contig"], categories=contigs["contig"].unique())
    return contigs, frags


# ------------------------------------------------------------------- ranking ----

def kmers(k):
    """Selective-base strings of length k, in msap_sim's itertools.product order."""
    out = []
    for letters in itertools.product("ACGT", repeat=k):
        out.append("".join(letters))
    return out


def pair_options(e_len, h_len):
    options = []
    for e in kmers(e_len):
        for h in kmers(h_len):
            options.append((e, h))
    return options


def pair_label(pair):
    return f"E-{pair[0]}/HM-{pair[1]}"


def top_pairs(ranking, n):
    lanes = []
    for e, h in zip(ranking["E_extension"].iloc[:n], ranking["HM_extension"].iloc[:n]):
        lanes.append((e, h))
    return lanes


def rank_pairs(frags, min_bp, e_len, h_len, bin_bp):
    """Same counts, tie order and sort as msap_sim.analyze, plus rank and score columns."""
    work = pd.DataFrame({
        "E_extension": frags["E_extension"].str.slice(0, e_len),
        "HM_extension": frags["HM_extension"].str.slice(0, h_len),
        "bin": (frags["bp"] - min_bp) // bin_bp,
        "contig": frags["contig"],
        "bp": frags["bp"],
    })
    stats = work.groupby(["E_extension", "HM_extension"], sort=False).agg(
        fragments=("bp", "size"),
        occupied_size_bins=("bin", "nunique"),
        unique_contigs=("contig", "nunique"),
        min_bp=("bp", "min"),
        max_bp=("bp", "max"),
    )
    grid = pd.MultiIndex.from_product([kmers(e_len), kmers(h_len)],
                                      names=["E_extension", "HM_extension"])
    table = stats.reindex(grid).reset_index()
    for col in ("fragments", "occupied_size_bins", "unique_contigs"):
        table[col] = table[col].fillna(0).astype("int64")
    table["co_migrating_excess"] = table["fragments"] - table["occupied_size_bins"]
    for col in ("min_bp", "max_bp"):
        table[col] = table[col].astype("Int64")
    table["score"] = table["occupied_size_bins"] - 2 * table["co_migrating_excess"]
    # stable descending sort on (score, bins) == msap_sim's sorted(..., reverse=True)
    order = np.lexsort((-table["occupied_size_bins"].to_numpy(), -table["score"].to_numpy()))
    table = table.iloc[order].reset_index(drop=True)
    table.insert(0, "rank", np.arange(1, len(table) + 1))
    return table[["rank"] + RANK_FIELDS + ["score"]]


# ------------------------------------------------------------------- exports ----

def csv_bytes(records, fieldnames):
    """CSV exactly as msap_sim.py writes it (csv.DictWriter: comma, minimal quoting, CRLF)."""
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(records)
    return buffer.getvalue().encode()


def fragments_csv(frags):
    return csv_bytes(frags[FRAGMENT_FIELDS].to_dict("records"), FRAGMENT_FIELDS)


def ranking_csv(ranking):
    records = []
    for row in ranking[RANK_FIELDS].itertuples(index=False):
        record = {}
        for name, value in zip(RANK_FIELDS, row):
            record[name] = "" if pd.isna(value) else value
        records.append(record)
    return csv_bytes(records, RANK_FIELDS)


def results_zip(run_info, frags, ranking):
    """CLI-format CSVs + run_info.json, rebuilt only when the run changes."""
    memo = st.session_state.get("_zip")
    if memo is not None and memo[0] == run_info:
        return memo[1]
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("predicted_fragments.csv", fragments_csv(frags))
        archive.writestr("primer_ranking.csv", ranking_csv(ranking))
        archive.writestr("run_info.json", run_info)
    data = buffer.getvalue()
    st.session_state["_zip"] = (run_info, data)
    return data


def make_run_info(src, genome, included, settings, frags, ranking):
    top = ranking.iloc[0]
    info = {
        "app": f"MSAP Streamlit app {APP_VERSION}",
        "model": "sequence-only, unmethylated, complete EcoRI + HpaII/MspI digest "
                 "(msap_sim.candidates)",
        "input": {"file": src["label"], "sha256": src["sha256"], "bytes": src["size"],
                  "loaded_utc": genome["loaded_utc"]},
        "contigs": {"included": included, "excluded": list(settings["excluded"])},
        "parameters": {"e_bases": settings["e_len"], "hm_bases": settings["h_len"],
                       "min_bp": settings["min_bp"], "max_bp": settings["max_bp"],
                       "bin_bp": settings["bin_bp"], "max_n": MAX_N},
        "results": {"fragments": int(len(frags)), "primer_pairs": int(len(ranking)),
                    "top_pair": (pair_label((top["E_extension"], top["HM_extension"]))
                                 if top["fragments"] > 0 else None)},
        "cli_equivalent": (f"python msap_sim.py --fasta {src['label']} "
                           f"--e-bases {settings['e_len']} --hm-bases {settings['h_len']} "
                           f"--min-bp {settings['min_bp']} --max-bp {settings['max_bp']} "
                           f"--bin-bp {settings['bin_bp']} --out results"),
    }
    if settings["excluded"]:
        info["cli_note"] = ("The CLI keeps every contig; drop the excluded ones from the "
                            "FASTA to reproduce this run.")
    return json.dumps(info, indent=2)


# ------------------------------------------------------------------- figures ----

def finish(fig, height, **layout):
    """Tight margins that still grow to fit tick labels and axis titles."""
    fig.update_layout(height=height, showlegend=False, margin=dict(l=10, r=10, t=10, b=10),
                      **layout)
    fig.update_xaxes(automargin=True)
    fig.update_yaxes(automargin=True)
    return fig


def heatmap_figure(ranking, e_len, h_len, metric):
    rows = kmers(e_len)
    cols = kmers(h_len)
    grid = pd.MultiIndex.from_product([rows, cols])
    values = ranking.set_index(["E_extension", "HM_extension"])[metric].reindex(grid)
    z = values.to_numpy(dtype=float).reshape(len(rows), len(cols))
    dtick = max(1, int(np.ceil((np.nanmax(z) - np.nanmin(z)) / 6)))
    top = ranking.iloc[0]
    fig = go.Figure()
    fig.add_trace(go.Heatmap(
        z=z, x=cols, y=rows, colorscale="Viridis",
        colorbar=dict(title=dict(text=METRICS[metric]), thickness=12, dtick=dtick),
        hovertemplate="E-%{y} / HM-%{x}<br>" + METRICS[metric] + ": %{z}<extra></extra>",
    ))
    fig.add_trace(go.Scatter(
        x=[top["HM_extension"]], y=[top["E_extension"]], mode="markers",
        marker=dict(symbol="star", size=14, color="white", line=dict(color="black", width=1)),
        hovertemplate="Top-ranked pair<extra></extra>", showlegend=False,
    ))
    font = 9 if len(rows) > 16 or len(cols) > 16 else 12
    finish(fig, int(np.clip(15 * len(rows) + 150, 340, 900)))
    fig.update_xaxes(title_text="HpaII/MspI-end selective bases", type="category",
                     side="top", tickfont=dict(size=font))
    fig.update_yaxes(title_text="EcoRI-end selective bases", type="category",
                     autorange="reversed", tickfont=dict(size=font))
    return fig


def ladder_sizes(lo, hi):
    step = 50 if hi - lo <= 600 else 100
    sizes = []
    size = -(-lo // step) * step
    while size <= hi:
        sizes.append(size)
        size += step
    return sizes


def add_bands(fig, x, sizes, colour, width, half=0.36):
    if len(sizes) == 0:
        return
    xs = []
    ys = []
    for size in sizes:
        xs.extend([x - half, x + half, None])
        ys.extend([size, size, None])
    fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", line=dict(color=colour, width=width),
                             hoverinfo="skip", showlegend=False))


def fragment_hover(sel, offset, crowd=None, bin_bp=None):
    texts = []
    for i, (contig, start, end, bp, orient) in enumerate(zip(
            sel["contig"], sel["start_1based"], sel["end_1based"], sel["bp"], sel["orientation"])):
        text = f"{bp} bp"
        if offset:
            text += f" → {bp + offset} bp amplicon"
        text += f"<br>{contig}:{start:,}–{end:,} ({orient})"
        if crowd is not None:
            text += f"<br>{crowd[i]} fragment(s) in this {bin_bp}-bp bin"
        texts.append(text)
    return texts


def gel_figure(view, lanes, min_bp, max_bp, bin_bp, offset):
    """Virtual gel: log-size axis, ladders at both edges, brighter bands where fragments co-migrate."""
    lo = min_bp + offset
    hi = max_bp + offset
    ladder = ladder_sizes(lo, hi)
    columns = [None] + list(lanes) + [None]
    crowded = len(lanes) > 6
    fig = go.Figure()
    labels = []
    for x, lane in enumerate(columns):
        if lane is None:
            add_bands(fig, x, ladder, LADDER_COLOUR, 2.5)
            labels.append("Ladder")
            continue
        sel = view[(view["E_sel"] == lane[0]) & (view["HM_sel"] == lane[1])]
        crowd = sel["bin"].map(sel["bin"].value_counts())
        level = crowd.clip(upper=4)
        for key, (alpha, width) in GEL_LEVELS.items():
            sizes = (sel.loc[level == key, "bp"] + offset).tolist()
            add_bands(fig, x, sizes, f"rgba({BAND_RGB},{alpha})", width)
        fig.add_trace(go.Scatter(
            x=[x] * len(sel), y=(sel["bp"] + offset).tolist(), mode="markers",
            marker=dict(size=12, color="rgba(0,0,0,0)"),
            text=fragment_hover(sel, offset, crowd.tolist(), bin_bp),
            hoverinfo="text", showlegend=False))
        labels.append(pair_label(lane) if crowded else f"E-{lane[0]}<br>HM-{lane[1]}")
    ticks = []
    for size in ladder:
        ticks.append(str(size))
    finish(fig, 600 if crowded else 560, plot_bgcolor=GEL_BACKGROUND, hovermode="closest")
    fig.update_xaxes(tickvals=list(range(len(columns))), ticktext=labels,
                     tickangle=-90 if crowded else 0, range=[-0.6, len(columns) - 0.4],
                     showgrid=False, zeroline=False, fixedrange=True)
    fig.update_yaxes(type="log", range=[np.log10(lo * 0.94), np.log10(hi * 1.05)],
                     tickvals=ladder, ticktext=ticks, showgrid=False, zeroline=False,
                     title_text="Amplicon length (bp)" if offset else "Fragment length (bp)")
    return fig


def profile_figure(sel, min_bp, max_bp, bin_bp, offset):
    """Electropherogram-style trace: one Gaussian per fragment, FWHM = bin width."""
    sigma = bin_bp / (2 * np.sqrt(2 * np.log(2)))
    step = 0.25 if sigma >= 1 else 0.1
    start = min_bp + offset - 3 * bin_bp
    stop = max_bp + offset + 3 * bin_bp
    centres = np.arange(start, stop + step / 2, step)
    edges = np.append(centres - step / 2, centres[-1] + step / 2)
    sizes = sel["bp"].to_numpy() + offset
    counts = np.histogram(sizes, bins=edges)[0].astype(float)
    half = int(np.ceil(4 * sigma / step))
    kernel = np.exp(-0.5 * (np.arange(-half, half + 1) * step / sigma) ** 2)
    signal = np.convolve(counts, kernel, mode="same")
    index = np.clip(np.round((sizes - start) / step).astype(int), 0, len(signal) - 1)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=centres, y=signal, mode="lines", fill="tozeroy",
                             line=dict(color=TRACE_COLOUR, width=1.6),
                             fillcolor="rgba(31,158,137,0.18)", hoverinfo="skip"))
    fig.add_trace(go.Scatter(x=sizes, y=signal[index], mode="markers",
                             marker=dict(size=7, color=PEAK_COLOUR),
                             text=fragment_hover(sel, offset), hoverinfo="text"))
    finish(fig, 340)
    fig.update_xaxes(title_text="Amplicon length (bp)" if offset else "Fragment length (bp)",
                     range=[start, stop])
    fig.update_yaxes(title_text="Relative signal", rangemode="tozero")
    return fig


def size_histogram(frags, min_bp, max_bp, bin_bp):
    width = max(bin_bp, (max_bp - min_bp) / 80)
    fig = go.Figure(go.Histogram(x=frags["bp"], marker_color=TRACE_COLOUR,
                                 xbins=dict(start=min_bp, end=max_bp + 1, size=width)))
    finish(fig, 280, bargap=0.05)
    fig.update_xaxes(title_text="Fragment length (bp)", range=[min_bp, max_bp + 1])
    fig.update_yaxes(title_text="Fragments")
    return fig


# ------------------------------------------------------------------------ UI ----

def local_path_allowed():
    """Server-side paths only make sense on your own machine.

    Hidden when Streamlit runs headless (Community Cloud, Linux servers);
    the environment variable MSAP_LOCAL_PATH=1 or =0 overrides.
    """
    flag = os.environ.get("MSAP_LOCAL_PATH")
    if flag in ("0", "1"):
        return flag == "1"
    return not st.get_option("server.headless")


def choose_source():
    """Sidebar source picker -> dict(label, sha256, size, open) or None."""
    digests = st.session_state.get("_digests")
    if digests is None:
        digests = {}
        st.session_state["_digests"] = digests
    sources = SOURCES if local_path_allowed() else SOURCES[:2]
    mode = st.radio("Source", sources, key="source")
    if mode == SOURCES[0]:
        data = demo_fasta()
        return {"label": "synthetic_demo.fasta", "sha256": hashlib.sha256(data).hexdigest(),
                "size": len(data), "open": lambda: io.BytesIO(data)}
    if mode == SOURCES[1]:
        upload = st.file_uploader(
            "Genome FASTA", type=["fa", "fasta", "fna", "fas", "gz", "zip"],
            help="Plain or gzip FASTA, or an NCBI Datasets .zip. Use an unmasked nuclear assembly.",
            **supported(st.file_uploader, max_upload_size=UPLOAD_LIMIT_MB))
        if upload is None:
            return None
        data = upload.getvalue()
        tag = (upload.name, len(data), getattr(upload, "file_id", getattr(upload, "id", "")))
        if tag not in digests:
            digests[tag] = hashlib.sha256(data).hexdigest()
        return {"label": upload.name, "sha256": digests[tag], "size": len(data),
                "open": lambda: io.BytesIO(data)}
    typed = st.text_input("FASTA path on this machine", key="path",
                          help="Plain, .gz or NCBI .zip; read in place, nothing is uploaded.")
    if not typed.strip():
        return None
    path = Path(typed.strip()).expanduser()
    if not path.is_file():
        st.error("File not found.")
        return None
    stat = path.stat()
    tag = (str(path.resolve()), stat.st_size, stat.st_mtime_ns)
    if tag not in digests:
        with st.spinner("Computing SHA-256…"):
            digests[tag] = sha256_file(path)
    return {"label": path.name, "sha256": digests[tag], "size": stat.st_size,
            "open": lambda: open(path, "rb")}


def load_genome(src):
    """Scan the FASTA once per file (keyed by SHA-256) and keep the result for this session."""
    memo = st.session_state.get("_genome")
    if memo is not None and memo["sha256"] == src["sha256"]:
        return memo
    message = f"Reading {src['label']}…"
    bar = st.progress(0.0, **supported(st.progress, text=message))

    def tick(frac, name):
        bar.progress(frac, **supported(st.progress, text=f"{message} {name}"))

    try:
        with src["open"]() as raw:
            contigs, frags = scan_genome(raw, src["size"], tick)
    except (OSError, EOFError, ValueError, zipfile.BadZipFile, zlib.error) as exc:
        bar.empty()
        st.error(f"Could not read {src['label']}: {exc}")
        st.stop()
    bar.empty()
    memo = {"sha256": src["sha256"], "contigs": contigs, "fragments": frags,
            "loaded_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    st.session_state["_genome"] = memo
    return memo


def sidebar_settings(src, contigs):
    names = contigs["contig"].unique().tolist()
    default = contigs.loc[contigs["organelle"], "contig"].unique().tolist()
    horizontal = supported(st.radio, horizontal=True)
    with st.sidebar:
        excluded = st.multiselect(
            "Exclude contigs", names, default=default, key="exclude_" + src["sha256"][:12],
            help="Organelle-like sequences (Mt/Pt, mitochondrion, chloroplast/plastid) "
                 "start excluded for a nuclear-only analysis.")
        st.header("Selective bases")
        e_len = st.radio("EcoRI end", [1, 2, 3], index=1, key="e_len", **horizontal)
        h_len = st.radio("HpaII/MspI end", [1, 2, 3], index=1, key="h_len", **horizontal)
        st.header("Fragment window")
        min_bp = st.number_input("Min bp", min_value=20, max_value=SIZE_CAP, value=100,
                                 step=10, key="min_bp")
        max_bp = st.number_input("Max bp", min_value=20, max_value=SIZE_CAP, value=500,
                                 step=10, key="max_bp")
        bin_bp = st.number_input(
            "Bin width (bp)", min_value=1, max_value=50, value=5, step=1, key="bin_bp",
            help="Approximate gel/capillary resolution; fragments in one bin co-migrate.")
        st.header("Display")
        offset = st.number_input(
            "Amplicon offset (bp)", min_value=0, max_value=200, value=0, step=1, key="offset",
            help="Added to sizes in the gel and profile views only. For primer cores ending "
                 "…AATTC and …CGG: len(E core) + len(HM core) − 6.")
    return {"excluded": excluded, "e_len": int(e_len), "h_len": int(h_len),
            "min_bp": int(min_bp), "max_bp": int(max_bp), "bin_bp": int(bin_bp),
            "offset": int(offset)}


def genome_warnings(contigs):
    total = int(contigs["length_bp"].sum())
    if total == 0:
        st.error("The included sequences are empty.")
        return
    n_pct = 100 * contigs["n"].sum() / total
    other_pct = 100 * (total - contigs["acgt"].sum() - contigs["n"].sum()) / total
    if other_pct > 10:
        st.error(f"{other_pct:.0f}% of characters are not A/C/G/T/N. Is this a protein FASTA?")
    if n_pct > 5:
        st.warning(f"{n_pct:.1f}% N: looks hard-masked or gappy. Use an unmasked nuclear assembly.")
    dupes = contigs.loc[contigs["contig"].duplicated(), "contig"].unique().tolist()
    if dupes:
        st.warning("Duplicate sequence IDs are filtered together: " + ", ".join(dupes[:5]))


def summary_metrics(frags, ranking):
    top = ranking.iloc[0]
    cols = st.columns(4)
    cols[0].metric("Fragments in window", f"{len(frags):,}")
    cols[1].metric("Primer pairs", f"{len(ranking):,}")
    cols[2].metric("Median fragments / pair", f"{ranking['fragments'].median():g}")
    label = pair_label((top["E_extension"], top["HM_extension"])) if top["fragments"] > 0 else "—"
    note = (f"Score {top['score']} = {top['occupied_size_bins']} bins − 2 × "
            f"{top['co_migrating_excess']} co-migrating")
    cols[3].metric("Top pair", label, **supported(st.metric, help=note))


def render_ranking(ranking, frags, run_info):
    shown = ranking
    most = int(ranking["fragments"].max())
    if most > 0:
        low, high = st.slider("Fragments per pair", 0, most, (0, most), key=f"frange_{most}")
        shown = ranking[ranking["fragments"].between(low, high)]
    order = ["rank", "E_extension", "HM_extension", "score", "fragments", "occupied_size_bins",
             "co_migrating_excess", "unique_contigs", "min_bp", "max_bp"]
    show_table(shown[order].rename(columns=RANK_LABELS).set_index("#"), height=420,
               hide_index=False)
    st.caption(f"{len(shown):,} of {len(ranking):,} pairs. Score = bins − 2 × co-migrating "
               "fragments, ties broken by bins. Click a column header to re-sort.")
    left, right = st.columns(2)
    download(left, "⬇ primer_ranking.csv", ranking_csv(ranking), "primer_ranking.csv",
             "text/csv")
    download(right, "⬇ All results (.zip)", results_zip(run_info, frags, ranking),
             "msap_results.zip", "application/zip")


def render_heatmap(ranking, e_len, h_len):
    metric = st.selectbox("Colour by", list(METRICS), format_func=METRICS.get,
                          key="heat_metric")
    show_chart(heatmap_figure(ranking, e_len, h_len, metric), "heatmap")
    st.caption("★ top-ranked pair. Rows: EcoRI-end selective bases; "
               "columns: HpaII/MspI-end selective bases.")


def render_gel(view, ranking, settings):
    e_len = settings["e_len"]
    h_len = settings["h_len"]
    mode = st.radio("Lanes", ["Top ranked", "Choose pairs"], key="lane_mode",
                    **supported(st.radio, horizontal=True))
    if mode == "Top ranked":
        most = min(16, len(ranking))
        count = st.slider("Number of lanes", 1, most, min(8, most), key="lane_n")
        lanes = top_pairs(ranking, count)
    else:
        lanes = st.multiselect("Primer pairs (up to 16)", pair_options(e_len, h_len),
                               default=top_pairs(ranking, 6), format_func=pair_label,
                               key=f"lanes_{e_len}{h_len}",
                               **supported(st.multiselect, max_selections=16))[:16]
    show_chart(gel_figure(view, lanes, settings["min_bp"], settings["max_bp"],
                          settings["bin_bp"], settings["offset"]), "gel")
    st.caption(f"Ladder lanes at both edges. Brighter, thicker bands = more fragments sharing a "
               f"{settings['bin_bp']}-bp bin. Under this unmethylated model the HpaII and MspI "
               "lanes are identical.")


def render_pair(view, ranking, settings):
    e_len = settings["e_len"]
    h_len = settings["h_len"]
    offset = settings["offset"]
    options = pair_options(e_len, h_len)
    rows = {}
    for record in ranking.to_dict("records"):
        rows[(record["E_extension"], record["HM_extension"])] = record

    def describe(pair):
        return f"{pair_label(pair)} · #{rows[pair]['rank']} · {rows[pair]['fragments']} fragments"

    top = (ranking["E_extension"].iloc[0], ranking["HM_extension"].iloc[0])
    pair = st.selectbox("Primer pair", options, index=options.index(top), format_func=describe,
                        key=f"pair_{e_len}{h_len}")
    row = rows[pair]
    st.caption(f"Rank {row['rank']} · {row['fragments']} fragments in "
               f"{row['occupied_size_bins']} bins · co-migrating {row['co_migrating_excess']} · "
               f"score {row['score']}")
    sel = view[(view["E_sel"] == pair[0]) & (view["HM_sel"] == pair[1])]
    if len(sel) == 0:
        st.info("No fragments for this pair in the current window.")
        return
    show_chart(profile_figure(sel, settings["min_bp"], settings["max_bp"], settings["bin_bp"],
                              offset), "profile")
    table = sel[FRAGMENT_FIELDS].copy()
    if offset:
        table["amplicon_bp"] = table["bp"] + offset
    show_table(table, height=min(38 + 35 * len(table), 400))
    download(st, "⬇ Fragments for this pair (.csv)", fragments_csv(sel),
             f"fragments_E-{pair[0]}_HM-{pair[1]}.csv", "text/csv")


def render_genome(src, contigs, frags_all, frags, settings, run_info):
    excluded = settings["excluded"]
    keep = ~contigs["contig"].isin(excluded)
    inc = contigs[keep]
    total = int(inc["length_bp"].sum())
    acgt = int(inc["acgt"].sum())
    cols = st.columns(4)
    cols[0].metric("Included sequence",
                   f"{total / 1e6:,.2f} Mb" if total >= 100_000 else f"{total:,} bp")
    cols[1].metric("GC", f"{100 * inc['gc'].sum() / acgt:.1f} %" if acgt else "—")
    cols[2].metric("EcoRI sites", f"{int(inc['EcoRI_sites'].sum()):,}")
    cols[3].metric("CCGG sites", f"{int(inc['CCGG_sites'].sum()):,}")
    window = frags_all[frags_all["bp"].between(settings["min_bp"], settings["max_bp"])]
    per_contig = window.groupby("contig", observed=False).size()
    table = pd.DataFrame({
        "Contig": contigs["contig"],
        "Included": keep,
        "Description": contigs["description"],
        "Length (bp)": contigs["length_bp"],
        "GC (%)": (100 * contigs["gc"] / contigs["acgt"].where(contigs["acgt"] > 0)).round(1),
        "N (%)": (100 * contigs["n"] / contigs["length_bp"].where(contigs["length_bp"] > 0)).round(2),
        "EcoRI sites": contigs["EcoRI_sites"],
        "CCGG sites": contigs["CCGG_sites"],
        "Fragments in window": contigs["contig"].map(per_contig).fillna(0).astype("int64"),
    })
    show_table(table.set_index("Contig"), height=min(38 + 35 * len(table), 420),
               hide_index=False)
    st.caption("Fragments in window counts every contig; excluded ones are left out of the ranking.")
    if len(frags):
        show_chart(size_histogram(frags, settings["min_bp"], settings["max_bp"],
                                  settings["bin_bp"]), "sizes")
    st.caption(f"{src['label']} · {src['size']:,} bytes · SHA-256 {src['sha256']}")
    with st.expander("Run record (run_info.json)"):
        st.code(run_info, language="json")


def model_notes():
    with st.expander("Model, assumptions and scoring"):
        st.markdown(
            "- Complete digestion of **unmethylated** DNA by EcoRI (G^AATTC) and HpaII/MspI "
            "(C^CGG), so HpaII and MspI lanes are identical; there is no methylation model.\n"
            "- A fragment is a pair of **adjacent** cut sites of different enzymes; any cut in "
            "between, by either enzyme, breaks it.\n"
            "- Selective bases are read inward from each site in the primer's 5′→3′ direction "
            "(reverse-complemented at right-hand sites).\n"
            "- Lengths run cut to cut, without adapters; the amplicon offset only shifts the gel "
            "and profile views.\n"
            "- No partial digestion, polymorphism, ligation/PCR bias or tissue effects. "
            "A high score is a hypothesis to test in the lab, not a validated assay.")
        st.markdown("For one primer pair with fragment lengths $L_i$, window start "
                    "$L_{\\min}$ and bin width $w$:")
        st.latex(r"b_i=\left\lfloor \frac{L_i-L_{\min}}{w}\right\rfloor,\qquad "
                 r"B=\left|\{b_i\}\right|,\qquad X=F-B,\qquad S=B-2X")
        st.markdown("$F$ = fragments, $B$ = occupied bins, $X$ = co-migrating excess, "
                    "$S$ = rank score (ties broken by the larger $B$). Heuristic, not calibrated.")
        st.markdown("The profile view sums one Gaussian per fragment, so each peak's FWHM "
                    "equals $w$; $\\delta$ is the amplicon offset:")
        st.latex(r"y(x)=\sum_i \exp\!\left[-\frac{(x-L_i-\delta)^2}{2\sigma^2}\right],\qquad "
                 r"\sigma=\frac{w}{2\sqrt{2\ln 2}}")
        st.markdown("For primer cores of length $\\ell_E$ (ending …AATTC) and $\\ell_{HM}$ "
                    "(ending …CGG), selective bases excluded:")
        st.latex(r"\delta=\ell_E+\ell_{HM}-6")


def main():
    st.set_page_config(page_title="MSAP primer screen", page_icon="🧬", layout="wide")
    st.title("In silico MSAP primer screen")
    st.caption("Potential EcoRI–HpaII/MspI fragments from sequence alone · complete-digest, "
               "unmethylated model · heuristic ranking, not calibrated")
    with st.sidebar:
        st.header("Genome")
        src = choose_source()
    if src is None:
        where = ", or a FASTA path on this machine" if local_path_allowed() else ""
        st.info("Choose a genome in the sidebar: the synthetic demo, an uploaded FASTA "
                f"(plain, .gz or NCBI Datasets .zip){where}.")
        model_notes()
        return
    genome = load_genome(src)
    contigs = genome["contigs"]
    frags_all = genome["fragments"]
    settings = sidebar_settings(src, contigs)
    if settings["min_bp"] > settings["max_bp"]:
        st.error("Min bp must not exceed max bp.")
        st.stop()

    included = contigs.loc[~contigs["contig"].isin(settings["excluded"]), "contig"]
    included = included.unique().tolist()
    keep = (frags_all["bp"].between(settings["min_bp"], settings["max_bp"])
            & frags_all["contig"].isin(included))
    frags = frags_all.loc[keep].reset_index(drop=True)
    ranking = rank_pairs(frags, settings["min_bp"], settings["e_len"], settings["h_len"],
                         settings["bin_bp"])
    view = frags.assign(
        E_sel=frags["E_extension"].str.slice(0, settings["e_len"]),
        HM_sel=frags["HM_extension"].str.slice(0, settings["h_len"]),
        bin=(frags["bp"] - settings["min_bp"]) // settings["bin_bp"],
    )
    run_info = make_run_info(src, genome, included, settings, frags, ranking)

    summary_metrics(frags, ranking)
    genome_warnings(contigs[contigs["contig"].isin(included)])
    tabs = st.tabs(["Ranking", "Heatmap", "Virtual gel", "Pair detail", "Genome"])
    with tabs[0]:
        render_ranking(ranking, frags, run_info)
    with tabs[1]:
        render_heatmap(ranking, settings["e_len"], settings["h_len"])
    with tabs[2]:
        render_gel(view, ranking, settings)
    with tabs[3]:
        render_pair(view, ranking, settings)
    with tabs[4]:
        render_genome(src, contigs, frags_all, frags, settings, run_info)
    model_notes()


if __name__ == "__main__":
    main()
