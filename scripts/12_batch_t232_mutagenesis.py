"""
Batch FoldX mutagenesis for one residue position.

Default task:
    5JZT T232 -> all standard amino acids except T

The script is intentionally staged:
1. Always run WT site checks.
2. Prepare FoldX BuildModel workspaces and individual_list.txt by default.
3. Only run FoldX when --run-foldx is passed.
4. Optionally run the structural feature pipeline and WT-mutant comparison.

Examples
--------
Prepare all T232 mutants without running FoldX:
    python scripts/12_batch_t232_mutagenesis.py

Run FoldX for all T232 mutants:
    python scripts/12_batch_t232_mutagenesis.py --run-foldx

Run only a subset:
    python scripts/12_batch_t232_mutagenesis.py --targets K R D E --run-foldx

Run FoldX, then extract features and compare against WT:
    python scripts/12_batch_t232_mutagenesis.py --run-foldx --run-pipeline --compare-wt
"""

import argparse
import copy
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from npstructfeat.foldx import (
    collect_foldx_mutant_model,
    prepare_foldx_buildmodel_workspace,
    run_foldx_buildmodel,
    validate_mutation_set_in_pdb,
)
from npstructfeat.io import load_config, resolve_path
from npstructfeat.mutation import (
    AA1_TO_AA3,
    check_mutation_set_sites,
    normalize_mutation_label,
    save_mutation_set_site_check,
    summarize_mutation_set_site_check,
)
from npstructfeat.utils import require_config_keys


DEFAULT_CHAINS = ["A", "B", "C", "D", "E", "F", "G"]
DEFAULT_AA_ORDER = list("ARNDCQEGHILKMFPSTWYV")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Batch prepare/run FoldX mutants at one residue position."
    )

    parser.add_argument(
        "--config",
        type=str,
        default="config/default.yaml",
        help="WT YAML config. Default: config/default.yaml",
    )
    parser.add_argument(
        "--original-aa",
        type=str,
        default="T",
        help="Original one-letter amino acid. Default: T",
    )
    parser.add_argument(
        "--position",
        type=int,
        default=232,
        help="Residue number to mutate. Default: 232",
    )
    parser.add_argument(
        "--targets",
        type=str,
        nargs="+",
        default=None,
        help=(
            "Target one-letter amino acids. Default: all 20 standard amino acids "
            "except the original amino acid."
        ),
    )
    parser.add_argument(
        "--include-original",
        action="store_true",
        help="Include the original amino acid as a no-op target.",
    )
    parser.add_argument(
        "--chains",
        type=str,
        nargs="+",
        default=DEFAULT_CHAINS,
        help="Chains to mutate. Default: A B C D E F G",
    )
    parser.add_argument(
        "--skip-workspace",
        action="store_true",
        help="Only run site checks and write configs; do not prepare FoldX workspaces.",
    )
    parser.add_argument(
        "--overwrite-individual-list",
        action="store_true",
        help="Overwrite existing individual_list.txt files.",
    )
    parser.add_argument(
        "--run-foldx",
        action="store_true",
        help="Actually run FoldX BuildModel for each mutant.",
    )
    parser.add_argument(
        "--run-pipeline",
        action="store_true",
        help="Run scripts/run_pipeline.py for successfully collected mutant models.",
    )
    parser.add_argument(
        "--compare-wt",
        action="store_true",
        help="Run scripts/10_compare_wt_mutant_features.py after --run-pipeline.",
    )
    parser.add_argument(
        "--wt-feature",
        type=str,
        default="data/processed/nanopore_features/5JZT_nanopore_structure_features.csv",
        help="WT nanopore feature CSV used by --compare-wt.",
    )
    parser.add_argument(
        "--config-output-dir",
        type=str,
        default="config",
        help="Directory for generated mutant YAML configs. Default: config",
    )
    parser.add_argument(
        "--summary-output",
        type=str,
        default="data/processed/mutation_sites/5JZT_T232_batch_mutagenesis_summary.csv",
        help="Batch summary CSV path.",
    )
    parser.add_argument(
        "--no-write-configs",
        action="store_true",
        help="Do not write mutant YAML config files.",
    )
    parser.add_argument(
        "--stop-on-error",
        action="store_true",
        help="Stop at the first failed mutant instead of continuing.",
    )

    return parser.parse_args()


def normalize_targets(original_aa: str, targets: List[str] | None, include_original: bool) -> List[str]:
    original_aa = original_aa.strip().upper()
    if original_aa not in AA1_TO_AA3:
        raise ValueError(f"Unknown original amino acid: {original_aa}")

    if targets is None:
        target_list = DEFAULT_AA_ORDER.copy()
    else:
        target_list = [target.strip().upper() for target in targets]

    invalid = [target for target in target_list if target not in AA1_TO_AA3]
    if invalid:
        raise ValueError(f"Unknown target amino acid(s): {invalid}")

    deduped = []
    for target in target_list:
        if target == original_aa and not include_original:
            continue
        if target not in deduped:
            deduped.append(target)

    return deduped


def build_mutation_ids(original_aa: str, position: int, targets: List[str]) -> List[str]:
    original_aa = original_aa.strip().upper()
    return [f"{original_aa}{position}{target}" for target in targets]


def get_base_pdb_id(config: Dict[str, Any]) -> str:
    template_pdb_id = config["input"].get("template_pdb_id")
    if template_pdb_id:
        return str(template_pdb_id).upper()
    return str(config["input"]["pdb_id"]).upper()


def build_mutant_config(config: Dict[str, Any], mutation_id: str) -> Dict[str, Any]:
    mutant_config = copy.deepcopy(config)
    base_pdb_id = get_base_pdb_id(config)
    mutation_label = normalize_mutation_label(mutation_id)
    mutant_pdb_id = f"{base_pdb_id}_{mutation_label}"

    modeled_root = config.get("foldx", {}).get("modeled_dir", "data/modeled")
    model_file = Path(modeled_root) / mutant_pdb_id / f"{mutant_pdb_id}_model.pdb"

    mutant_config["input"]["pdb_id"] = mutant_pdb_id
    mutant_config["input"]["nanopore_id"] = f"aerolysin_{mutation_label}"
    mutant_config["input"]["structure_file"] = model_file.as_posix()
    mutant_config["input"]["file_format"] = "pdb"
    mutant_config["input"]["assembly_type"] = "modeled_biological_assembly"
    mutant_config["input"]["structure_source_type"] = "template_based_mutant_model"
    mutant_config["input"]["template_pdb_id"] = base_pdb_id
    mutant_config["input"]["modeling_method"] = "FoldX"
    mutant_config["input"]["mutation_list"] = [mutation_id]

    return mutant_config


def write_mutant_config(
    config: Dict[str, Any],
    mutation_id: str,
    output_dir: str | Path,
) -> Path:
    output_dir = resolve_path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    base_pdb_id = get_base_pdb_id(config)
    mutation_label = normalize_mutation_label(mutation_id)
    output_file = output_dir / f"{base_pdb_id}_{mutation_label}_foldx.yaml"

    mutant_config = build_mutant_config(config=config, mutation_id=mutation_id)
    with output_file.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(mutant_config, handle, sort_keys=False, allow_unicode=False)

    return output_file


def run_command(cmd: List[str]) -> Dict[str, Any]:
    completed = subprocess.run(cmd, capture_output=True, text=True)
    return {
        "cmd": cmd,
        "returncode": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
    }


def run_pipeline_for_config(config_file: Path) -> Dict[str, Any]:
    return run_command(
        [sys.executable, "scripts/run_pipeline.py", "--config", str(config_file)]
    )


def compare_wt_mutant(
    wt_feature: str,
    mutation_id: str,
    mutant_config: Dict[str, Any],
) -> Dict[str, Any]:
    mutant_pdb_id = mutant_config["input"]["pdb_id"]
    mutant_id = mutant_config["input"]["nanopore_id"]
    template_pdb_id = mutant_config["input"]["template_pdb_id"]
    mutant_feature = (
        f"data/processed/nanopore_features/"
        f"{mutant_pdb_id}_nanopore_structure_features.csv"
    )

    return run_command(
        [
            sys.executable,
            "scripts/10_compare_wt_mutant_features.py",
            "--wt-feature",
            wt_feature,
            "--mutant-feature",
            mutant_feature,
            "--wt-id",
            "aerolysin_WT",
            "--mutant-id",
            mutant_id,
            "--template-pdb-id",
            template_pdb_id,
            "--mutation",
            mutation_id,
            "--modeling-method",
            "FoldX",
        ]
    )


def append_status(row: Dict[str, Any], status: str, message: str = "") -> Dict[str, Any]:
    row["status"] = status
    row["message"] = message
    return row


def main() -> None:
    args = parse_args()

    config = load_config(args.config)
    require_config_keys(
        config,
        required_keys=["project", "input", "structure", "pore", "output", "save"],
    )

    original_aa = args.original_aa.strip().upper()
    targets = normalize_targets(
        original_aa=original_aa,
        targets=args.targets,
        include_original=args.include_original,
    )
    mutation_ids = build_mutation_ids(
        original_aa=original_aa,
        position=args.position,
        targets=targets,
    )

    print("Batch mutagenesis settings")
    print(f"Config: {args.config}")
    print(f"Position: {original_aa}{args.position}")
    print(f"Targets: {' '.join(targets)}")
    print(f"Mutants: {len(mutation_ids)}")
    print(f"Chains: {' '.join(args.chains)}")
    print(f"Run FoldX: {args.run_foldx}")
    print(f"Run pipeline: {args.run_pipeline}")
    print(f"Compare WT: {args.compare_wt}")

    rows = []

    for index, mutation_id in enumerate(mutation_ids, start=1):
        mutation_label = normalize_mutation_label(mutation_id)
        row: Dict[str, Any] = {
            "mutation_id": mutation_id,
            "mutation_label": mutation_label,
            "target_aa": mutation_id[-1],
            "site_check_file": "",
            "site_check_passed": False,
            "foldx_workspace": "",
            "individual_list_file": "",
            "config_file": "",
            "model_file": "",
            "foldx_returncode": "",
            "validation_file": "",
            "validation_passed": False,
            "pipeline_returncode": "",
            "compare_returncode": "",
        }

        print(f"\n[{index}/{len(mutation_ids)}] {mutation_id}")

        try:
            site_df = check_mutation_set_sites(
                config=config,
                mutation_id=mutation_id,
                chains=args.chains,
            )
            summary = summarize_mutation_set_site_check(site_df)
            site_check_file = save_mutation_set_site_check(
                config=config,
                mutation_id=mutation_id,
                chains=args.chains,
            )
            row["site_check_file"] = str(site_check_file)
            row["site_check_passed"] = bool(
                summary["all_sites_found"] and summary["all_expected_match"]
            )
            print(f"Site check: {row['site_check_passed']} -> {site_check_file}")

            if not args.no_write_configs:
                config_file = write_mutant_config(
                    config=config,
                    mutation_id=mutation_id,
                    output_dir=args.config_output_dir,
                )
                row["config_file"] = str(config_file)
                print(f"Config: {config_file}")
            else:
                config_file = None

            if not args.skip_workspace:
                workspace = prepare_foldx_buildmodel_workspace(
                    config=config,
                    mutation_id=mutation_id,
                    chains=args.chains,
                    overwrite_individual_list=args.overwrite_individual_list,
                )
                row["foldx_workspace"] = str(workspace["mutant_dir"])
                row["individual_list_file"] = str(workspace["individual_list_file"])
                print(f"Workspace: {workspace['mutant_dir']}")

            if args.run_foldx:
                foldx_result = run_foldx_buildmodel(
                    config=config,
                    mutation_id=mutation_id,
                    chains=args.chains,
                    overwrite_individual_list=args.overwrite_individual_list,
                )
                row["foldx_returncode"] = foldx_result["returncode"]
                print(f"FoldX return code: {foldx_result['returncode']}")

                if foldx_result["returncode"] != 0:
                    raise RuntimeError(
                        "FoldX BuildModel failed. "
                        f"stderr: {foldx_result['stderr'][:500]}"
                    )

                model_file = collect_foldx_mutant_model(
                    config=config,
                    mutation_id=mutation_id,
                )
                row["model_file"] = str(model_file)
                print(f"Model: {model_file}")

                validation_df = validate_mutation_set_in_pdb(
                    pdb_file=model_file,
                    mutation_id=mutation_id,
                    chains=args.chains,
                )
                validation_passed = bool(validation_df["is_target_match"].all())
                row["validation_passed"] = validation_passed

                validation_file = (
                    Path(row["foldx_workspace"])
                    / f"{mutation_label}_mutant_validation.csv"
                )
                validation_df.to_csv(validation_file, index=False, encoding="utf-8-sig")
                row["validation_file"] = str(validation_file)
                print(f"Validation: {validation_passed} -> {validation_file}")

                if not validation_passed:
                    raise RuntimeError("Mutant validation failed.")

            if args.run_pipeline:
                if config_file is None:
                    raise RuntimeError("--run-pipeline requires generated config files.")

                pipeline_result = run_pipeline_for_config(config_file)
                row["pipeline_returncode"] = pipeline_result["returncode"]
                print(f"Pipeline return code: {pipeline_result['returncode']}")

                if pipeline_result["returncode"] != 0:
                    raise RuntimeError(
                        "Feature pipeline failed. "
                        f"stderr: {pipeline_result['stderr'][:500]}"
                    )

            if args.compare_wt:
                if not args.run_pipeline:
                    raise RuntimeError("--compare-wt requires --run-pipeline.")
                if config_file is None:
                    raise RuntimeError("--compare-wt requires generated config files.")

                mutant_config = build_mutant_config(config=config, mutation_id=mutation_id)
                compare_result = compare_wt_mutant(
                    wt_feature=args.wt_feature,
                    mutation_id=mutation_id,
                    mutant_config=mutant_config,
                )
                row["compare_returncode"] = compare_result["returncode"]
                print(f"Compare return code: {compare_result['returncode']}")

                if compare_result["returncode"] != 0:
                    raise RuntimeError(
                        "WT-mutant comparison failed. "
                        f"stderr: {compare_result['stderr'][:500]}"
                    )

            append_status(row, "ok")

        except Exception as exc:
            append_status(row, "failed", str(exc))
            print(f"FAILED: {exc}")
            rows.append(row)
            if args.stop_on_error:
                break
            continue

        rows.append(row)

    summary_output = resolve_path(args.summary_output)
    summary_output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(summary_output, index=False, encoding="utf-8-sig")

    ok_count = sum(1 for row in rows if row.get("status") == "ok")
    failed_count = sum(1 for row in rows if row.get("status") == "failed")

    print("\nBatch summary")
    print(f"OK: {ok_count}")
    print(f"Failed: {failed_count}")
    print(f"Summary CSV: {summary_output}")


if __name__ == "__main__":
    main()
