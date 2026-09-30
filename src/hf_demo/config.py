"""Configuration loading and validation for the public demonstration."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml


REQUIRED_KEYS = {
    "domain": ("length_x", "length_y", "cells_x", "cells_y"),
    "material": (
        "young_modulus",
        "poisson_ratio",
        "fracture_toughness",
        "regularization_length",
        "residual_stiffness",
        "biot_coefficient",
    ),
    "flow": (
        "storage",
        "matrix_permeability",
        "fracture_permeability",
        "viscosity",
        "source_amplitude",
        "source_width",
    ),
    "time": (
        "step_size",
        "number_of_steps",
        "staggered_tolerance",
        "staggered_max_iterations",
        "output_interval",
    ),
    "precrack": ("half_length", "center_x", "center_y"),
    "output": ("directory",),
}


def load_config(path: str | Path) -> dict[str, Any]:
    """Load a YAML configuration and return a validated independent mapping."""
    config_path = Path(path)
    if not config_path.is_file():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    with config_path.open("r", encoding="utf-8") as stream:
        raw = yaml.safe_load(stream)

    if not isinstance(raw, dict):
        raise ValueError("The configuration root must be a YAML mapping.")

    config = deepcopy(raw)
    validate_config(config)
    return config


def validate_config(config: dict[str, Any]) -> None:
    """Raise ValueError when a required or physical constraint is invalid."""
    if not str(config.get("case_name", "")).strip():
        raise ValueError("case_name must be a non-empty string.")

    for section, keys in REQUIRED_KEYS.items():
        if section not in config or not isinstance(config[section], dict):
            raise ValueError(f"Missing configuration section: {section}")
        missing = [key for key in keys if key not in config[section]]
        if missing:
            raise ValueError(f"Missing keys in {section}: {', '.join(missing)}")

    domain = config["domain"]
    material = config["material"]
    flow = config["flow"]
    time = config["time"]
    crack = config["precrack"]

    positive = {
        "domain.length_x": domain["length_x"],
        "domain.length_y": domain["length_y"],
        "domain.cells_x": domain["cells_x"],
        "domain.cells_y": domain["cells_y"],
        "material.young_modulus": material["young_modulus"],
        "material.fracture_toughness": material["fracture_toughness"],
        "material.regularization_length": material["regularization_length"],
        "material.residual_stiffness": material["residual_stiffness"],
        "flow.storage": flow["storage"],
        "flow.matrix_permeability": flow["matrix_permeability"],
        "flow.fracture_permeability": flow["fracture_permeability"],
        "flow.viscosity": flow["viscosity"],
        "flow.source_width": flow["source_width"],
        "time.step_size": time["step_size"],
        "time.number_of_steps": time["number_of_steps"],
        "time.staggered_tolerance": time["staggered_tolerance"],
        "time.staggered_max_iterations": time["staggered_max_iterations"],
        "time.output_interval": time["output_interval"],
        "precrack.half_length": crack["half_length"],
    }
    invalid = [name for name, value in positive.items() if float(value) <= 0.0]
    if invalid:
        raise ValueError("Values must be positive: " + ", ".join(invalid))

    nu = float(material["poisson_ratio"])
    if not (-1.0 < nu < 0.5):
        raise ValueError("material.poisson_ratio must satisfy -1 < nu < 0.5.")

    alpha = float(material["biot_coefficient"])
    if not (0.0 <= alpha <= 1.0):
        raise ValueError("material.biot_coefficient must lie in [0, 1].")

    if float(flow["fracture_permeability"]) < float(flow["matrix_permeability"]):
        raise ValueError("fracture_permeability must be >= matrix_permeability.")

    half_x = 0.5 * float(domain["length_x"])
    half_y = 0.5 * float(domain["length_y"])
    cx = float(crack["center_x"])
    cy = float(crack["center_y"])
    length = float(crack["half_length"])
    if abs(cx) + length >= half_x or abs(cy) >= half_y:
        raise ValueError("The horizontal pre-crack must lie strictly inside the domain.")
