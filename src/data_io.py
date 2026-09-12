"""Configuration, source integrity and generated-artifact helpers."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_config():
    return json.loads((ROOT / "config/model_config.json").read_text(encoding="utf-8"))


def verify_sources():
    manifest = json.loads((ROOT / "reports/input_manifest.json").read_text(encoding="utf-8"))
    for item in manifest["files"]:
        if digest(item["path"]) != item["sha256"]:
            raise ValueError("Input differs from stage0 manifest: " + item["path"])
    return manifest


def save_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def source_fingerprint():
    paths = sorted((ROOT / "src").glob("*.py"))
    paths += [ROOT / "config/model_config.json"]
    return {str(p.relative_to(ROOT)): digest(p) for p in paths}
