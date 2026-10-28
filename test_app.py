"""App <-> CLI parity and app smoke tests:  python test_app.py"""
import contextlib
import gzip
import io
import os
import random
import tempfile
import zipfile
from pathlib import Path

import app
from msap_sim import analyze

HERE = Path(__file__).resolve().parent
DEMO = HERE / "synthetic_demo.fasta"
CASES = [(100, 500, 2, 2, 5), (100, 500, 3, 3, 5), (50, 800, 1, 3, 1),
         (200, 300, 3, 1, 10), (20, 1500, 2, 3, 7)]


def write_random_genome(path, seed=7, contigs=4, length=60000):
    """Contigs with extra sites, Ns, a soft-masked stretch, descriptions and a 'Mt' organelle."""
    rng = random.Random(seed)
    lines = []
    for k in range(contigs + 1):
        seq = rng.choices("ACGT", weights=[28, 22, 22, 28], k=length)
        for _ in range(length // 270):
            motif = rng.choice(["GAATTC", "CCGG", "CCGG"])
            pos = rng.randrange(length - len(motif))
            seq[pos:pos + len(motif)] = list(motif)
        for _ in range(40):
            seq[rng.randrange(length)] = "N"
        text = "".join(seq)
        text = text[:5000] + text[5000:9000].lower() + text[9000:]
        if k == contigs:
            lines.append(">Mt Oryza-like mitochondrion, complete genome")
        else:
            lines.append(f">chr{k + 1} chromosome {k + 1}, test assembly")
        for i in range(0, len(text), 70):
            lines.append(text[i:i + 70])
    path.write_text("\n".join(lines) + "\n")


def app_outputs(path, min_bp, max_bp, e_len, h_len, bin_bp):
    with open(path, "rb") as raw:
        contigs, frags = app.scan_genome(raw, path.stat().st_size)
    window = frags[frags["bp"].between(min_bp, max_bp)].reset_index(drop=True)
    ranking = app.rank_pairs(window, min_bp, e_len, h_len, bin_bp)
    return app.fragments_csv(window), app.ranking_csv(ranking)


def cli_outputs(path, min_bp, max_bp, e_len, h_len, bin_bp):
    with tempfile.TemporaryDirectory() as out:
        with contextlib.redirect_stdout(io.StringIO()):
            analyze(path, min_bp, max_bp, e_len, h_len, bin_bp, out)
        folder = Path(out)
        return ((folder / "predicted_fragments.csv").read_bytes(),
                (folder / "primer_ranking.csv").read_bytes())


def test_demo_regeneration():
    original = app.DEMO_FILE
    app.DEMO_FILE = HERE / "does_not_exist.fasta"
    try:
        assert app.demo_fasta() == DEMO.read_bytes()
    finally:
        app.DEMO_FILE = original


def test_cli_parity(tmp):
    genome = tmp / "random.fa"
    write_random_genome(genome)
    zipped = tmp / "random.fa.gz"
    zipped.write_bytes(gzip.compress(genome.read_bytes()))
    empty = tmp / "no_sites.fa"
    empty.write_text(">blank\n" + "ACGT" * 500 + "\n>also_blank\nTTTTAAAA\n")
    checked = 0
    for fasta in (DEMO, genome, zipped, empty):
        for case in CASES:
            assert app_outputs(fasta, *case) == cli_outputs(fasta, *case), (fasta.name, case)
            checked += 1
    return checked


def test_archives(tmp):
    genome = tmp / "random.fa"
    reference = app_outputs(genome, *CASES[0])
    bundle = tmp / "ncbi_dataset.zip"
    with zipfile.ZipFile(bundle, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("README.md", "not a genome")
        archive.writestr("ncbi_dataset/data/GCF_x/cds_from_genomic.fna", ">cds\nGAATTCAAAA\n")
        archive.write(genome, "ncbi_dataset/data/GCF_x/GCF_x_genomic.fna")
    assert app_outputs(bundle, *CASES[0]) == reference
    nested = tmp / "nested.zip"
    with zipfile.ZipFile(nested, "w") as archive:
        archive.writestr("genome.fna.gz", gzip.compress(genome.read_bytes()))
    assert app_outputs(nested, *CASES[0]) == reference
    with open(genome, "rb") as raw:
        contigs, _ = app.scan_genome(raw, genome.stat().st_size)
    assert contigs.loc[contigs["organelle"], "contig"].tolist() == ["Mt"]


def test_app_runs(tmp):
    try:
        from streamlit.testing.v1 import AppTest
    except ImportError:
        print("AppTest needs Streamlit >= 1.28; skipped")
        return
    os.environ["MSAP_LOCAL_PATH"] = "0"        # as on Community Cloud: no server paths
    at = AppTest.from_file(str(HERE / "app.py"), default_timeout=120)
    at.run()
    assert "Local path" not in at.radio(key="source").options
    if hasattr(AppTest, "file_uploader"):      # end-to-end upload of a gzipped genome
        at.radio(key="source").set_value("Upload FASTA").run()
        gz = (tmp / "random.fa.gz").read_bytes()
        at.file_uploader[0].upload("random.fa.gz", gz, "application/gzip").run()
        assert not at.exception, at.exception
        assert at.multiselect[0].value == ["Mt"], at.multiselect[0].value
        assert at.metric[0].value != "80", "upload should replace the demo genome"
    os.environ["MSAP_LOCAL_PATH"] = "1"        # own machine
    at = AppTest.from_file(str(HERE / "app.py"), default_timeout=120)
    at.run()
    assert not at.exception, at.exception
    assert "Local path" in at.radio(key="source").options
    assert at.metric[0].value == "80", at.metric[0].value
    at.radio(key="e_len").set_value(3).run()
    at.radio(key="h_len").set_value(3).run()
    at.number_input(key="bin_bp").set_value(1).run()
    at.number_input(key="offset").set_value(26).run()
    at.radio(key="lane_mode").set_value("Choose pairs").run()
    assert not at.exception, at.exception
    at.selectbox(key="heat_metric").set_value("fragments").run()
    at.selectbox(key="pair_33").set_value(("TGG", "CAG")).run()
    assert not at.exception, at.exception
    at.number_input(key="min_bp").set_value(600).run()
    at.number_input(key="max_bp").set_value(700).run()
    assert not at.exception, at.exception
    assert at.metric[0].value == "0", at.metric[0].value
    at.number_input(key="max_bp").set_value(300).run()
    assert at.error, "min > max should show an error"
    at.number_input(key="max_bp").set_value(1500).run()
    at.number_input(key="min_bp").set_value(20).run()
    gz = tmp / "random.fa.gz"
    at.radio(key="source").set_value("Local path").run()
    at.text_input(key="path").input(str(gz)).run()
    assert not at.exception, at.exception
    excluded = at.multiselect[0].value
    assert excluded == ["Mt"], excluded
    at.multiselect[0].set_value([]).run()
    at.radio(key="e_len").set_value(1).run()
    at.radio(key="h_len").set_value(1).run()
    assert not at.exception, at.exception
    at.text_input(key="path").input(str(tmp / "missing.fa")).run()
    assert at.error, "missing path should show an error"


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as folder:
        workdir = Path(folder)
        test_demo_regeneration()
        n = test_cli_parity(workdir)
        test_archives(workdir)
        test_app_runs(workdir)
    print(f"app tests passed ({n} byte-identical CLI comparisons)")
