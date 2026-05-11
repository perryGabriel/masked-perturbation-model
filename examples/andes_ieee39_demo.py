"""Small script demonstrating the packaged ANDES IEEE39 example case.

Run with ANDES installed to build the real ANDES system and linearized LFT-style
container:

    python examples/andes_ieee39_demo.py --build
"""

from __future__ import annotations

import argparse

from masked_perturbation_model.cases import load_ieee39_case


def main() -> None:
    parser = argparse.ArgumentParser(description="Demonstrate the ANDES IEEE39 case API.")
    parser.add_argument("--build", action="store_true", help="Load ANDES and build the linearized model.")
    args = parser.parse_args()

    case = load_ieee39_case()
    print("IEEE39 case metadata:")
    for key, value in case.summary().items():
        print(f"  {key}: {value}")

    if not args.build:
        print("\nPass --build after installing mpmgame[andes] to construct the ANDES system.")
        return

    system = case.build_system()
    lft = case.build_lft(system_model=system)
    print("\nBuilt ANDES system:")
    print(system.summary())
    print("\nLFT-style representation:")
    print({"nstates": lft.nstates, "ninputs": lft.ninputs, "noutputs": lft.noutputs})


if __name__ == "__main__":
    main()
