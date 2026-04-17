#!/usr/bin/env python3
"""Run MRILAB3 tests through matlab.engine.

Usage examples:
  python test_pipeline_with_matlab_engine.py --mode smoke
  python test_pipeline_with_matlab_engine.py --mode full
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _matlab_str(path: Path) -> str:
    return str(path).replace("\\", "/").replace("'", "''")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Test MRILAB3 with matlab.engine")
    parser.add_argument(
        "--mode",
        choices=("smoke", "full"),
        default="smoke",
        help="smoke: launch MATLAB + load config; full: run complete pipeline",
    )
    parser.add_argument(
        "--repo-root",
        default=str(Path(__file__).resolve().parent),
        help="Absolute path to MRILAB3 repository root",
    )
    parser.add_argument(
        "--print-config-json",
        action="store_true",
        help="Print selected cfg.templates fields as JSON in smoke mode",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    if not repo_root.is_dir():
        print(f"[error] repo root not found: {repo_root}", file=sys.stderr)
        return 1

    try:
        import matlab.engine  # type: ignore
    except Exception as exc:  # pragma: no cover
        print(
            "[error] matlab.engine is not available. "
            "Install from your MATLAB: <matlabroot>/extern/engines/python",
            file=sys.stderr,
        )
        print(f"[detail] {exc}", file=sys.stderr)
        return 2

    eng = None
    try:
        print("[info] starting MATLAB engine ...")
        eng = matlab.engine.start_matlab()
        repo_matlab = _matlab_str(repo_root)
        eng.eval(f"cd('{repo_matlab}');", nargout=0)
        eng.eval(f"addpath(genpath('{repo_matlab}'));", nargout=0)
        print("[info] MATLAB started, repository path added.")

        if args.mode == "smoke":
            eng.eval("cfg = task_fmri_pipeline_config();", nargout=0)
            eng.eval("disp(cfg.templates);", nargout=0)
            if args.print_config_json:
                payload = {
                    "installRoots": eng.eval("cfg.templates.installRoots", nargout=1),
                    "installRootEnvVars": eng.eval("cfg.templates.installRootEnvVars", nargout=1),
                    "searchDirs": eng.eval("cfg.templates.searchDirs", nargout=1),
                }
                print(json.dumps(payload, ensure_ascii=False, default=str, indent=2))
            print("[ok] smoke test passed.")
            return 0

        print("[info] running full pipeline (run_task_fmri_pipeline) ...")
        eng.run_task_fmri_pipeline(nargout=0)
        print("[ok] full pipeline run finished.")
        return 0
    except Exception as exc:  # pragma: no cover
        print(f"[error] matlab.engine run failed: {exc}", file=sys.stderr)
        return 3
    finally:
        if eng is not None:
            try:
                eng.quit()
            except Exception:
                pass


if __name__ == "__main__":
    raise SystemExit(main())
