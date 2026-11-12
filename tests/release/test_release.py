"""Validate runtime resources, graph search, scoring, and the configurable entry point."""

import importlib
import importlib.util
from importlib import resources
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.sparse import lil_matrix


ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize("name", ["basic.constant", "basic.vdmer", "database.core", "database.database_cluster",
                                  "search.search", "search.search_selfcenter", "combs.search_lig_indep_2ndshell"])
def test_public_modules_import(name):
    assert importlib.import_module("metalprot." + name)


def test_runtime_resources_present():
    root = resources.files("metalprot") / "constants"
    for name in ["APBLE.txt", "abple_dict.pkl", "tetrahydral_geo.pdb", "Zn_Tet_O.pdb", "ideal_ala.pdb"]:
        assert (root / name).is_file(), name


def test_graph_requires_complete_triplet():
    from metalprot.search.find_path_by_matrix import calc_adj_matrix_paths
    graph = lil_matrix((4, 4), dtype=bool)
    for i, j in [(0, 1), (0, 2), (1, 2), (2, 3)]:
        graph[i, j] = True
    assert calc_adj_matrix_paths(graph, 3) == [[0, 1, 2]]


def test_native_contact_atoms():
    import prody as pr
    from metalprot.basic.vdmer import get_metal_contact_atoms
    target = pr.parsePDB(str(ROOT / "tests/test_data/5od1_zn.pdb"))
    contacts = get_metal_contact_atoms(target)
    assert [a.getName() for a in contacts] == ["ZN", "ND1", "NE2", "NE2"]
    assert [a.getResnum() for a in contacts[1:]] == [35, 61, 65]


def test_combination_score_reference():
    from metalprot.search.search import Search_vdM
    search = Search_vdM.__new__(Search_vdM)
    search.density_radius = 1.0
    search.aa_vdm_info_dict = {"H": [60, 1, 60, 60]}
    search.vdms = [SimpleNamespace(score=0)]
    search.cluster_centroid_dict = {("H", 0): 0}
    key = ((0, 1, 2), (("H", 0), ("H", 0), ("H", 0)))
    info = SimpleNamespace(query_dict={0: [0]*2, 1: [0]*4, 2: [0]*8},
                           overlap_query_id_dict=None,
                           centroid_dict={i: SimpleNamespace(clu_num=60) for i in range(3)})
    search.neighbor_calc_comb_score({key: info})
    assert info.cluScore == pytest.approx(0)
    assert info.overlapScore == pytest.approx(3 / math.pi)


def test_separate_open_site_geometry(tmp_path, monkeypatch):
    import prody as pr
    from metalprot.basic import constant
    from metalprot.basic.filter import Search_filter
    from metalprot.search.search import Search_vdM
    geometry = constant.tetrahydra_geo_o.copy()
    geometry.setCoords(geometry.getCoords() + np.array([10., 0., 0.]))
    path = tmp_path / "open_site.pdb"
    pr.writePDB(str(path), geometry)
    monkeypatch.setattr(Search_vdM, "setup", lambda self: None)
    search = Search_vdM(str(ROOT / "tests/test_data/5od1_zn.pdb"), str(tmp_path / "run"), [], {}, None,
                        search_filter=Search_filter(), geometry_o_path=str(path))
    np.testing.assert_allclose(search.geo_o_struct.getCoords(), geometry.getCoords(), atol=0.001)


def load_runner():
    spec = importlib.util.spec_from_file_location("release_runner", ROOT / "scripts/run_search.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("radius", [0, -0.1, float("nan")])
def test_bad_density_radius_rejected(tmp_path, radius):
    config = json.loads((ROOT / "examples/basic_search/config.json").read_text())
    config["density_radius"] = radius
    path = tmp_path / "invalid.json"
    path.write_text(json.dumps(config))
    with pytest.raises(ValueError, match="density_radius"):
        load_runner().load_config(path)


def test_config_paths_do_not_depend_on_cwd(monkeypatch, tmp_path):
    # Check path resolution and file existence without unpickling a database.
    runner = load_runner()
    config_dir = tmp_path / "config"
    database = tmp_path / "database"
    config_dir.mkdir()
    database.mkdir()
    for name in runner.DATABASE_FILES:
        (database / name).write_bytes(b"path-resolution-fixture")
    settings = json.loads((ROOT / "examples/basic_search/config.json").read_text())
    settings["target"] = str(ROOT / "examples/basic_search/input.pdb")
    settings["database"] = "../database"
    path = config_dir / "config.json"
    path.write_text(json.dumps(settings))
    monkeypatch.chdir(ROOT)
    config, inputs = runner.load_config(path)
    assert Path(config["target"]) == ROOT / "examples/basic_search/input.pdb"
    assert Path(config["database"]) == database
    assert all(p.is_file() for p in inputs)


def test_missing_external_database_has_actionable_error(tmp_path):
    settings = json.loads((ROOT / "examples/basic_search/config.json").read_text())
    settings["target"] = str(ROOT / "examples/basic_search/input.pdb")
    settings["database"] = "missing_database"
    path = tmp_path / "config.json"
    path.write_text(json.dumps(settings))
    with pytest.raises(FileNotFoundError, match="Database files are not distributed with the code"):
        load_runner().load_config(path)
