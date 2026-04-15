import json
from pathlib import Path


def _all_code_cells_compile(nb_path: str):
    nb = json.loads(Path(nb_path).read_text())
    for i, cell in enumerate(nb.get("cells", [])):
        if cell.get("cell_type") == "code":
            compile("".join(cell.get("source", [])), f"{nb_path}:cell{i}", "exec")


def test_ch2_notebooks_parse_and_compile():
    _all_code_cells_compile("notebooks/chapter2_bounds_demo.ipynb")
    _all_code_cells_compile("notebooks/chapter2_design_demo.ipynb")


def test_report_exists():
    p = Path("reports/chapter2_demo_report.md")
    assert p.exists()
    assert "Chapter-2" in p.read_text()
