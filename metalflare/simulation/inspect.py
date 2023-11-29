import argparse
import os
import subprocess
import tempfile

from jinja2 import Environment, FileSystemLoader


def cli_vmd_inspect() -> None:
    r"""VMD inspection of MD setup."""
    parser = argparse.ArgumentParser(
        description="Inspect system for MD simulations using VMD."
    )
    parser.add_argument(
        "topo",
        type=str,
        nargs="?",
        help="Path to topology file.",
    )
    parser.add_argument(
        "coord",
        type=str,
        nargs="?",
        help="Path to coordinate file.",
    )
    parser.add_argument(
        "--topo_type", type=str, nargs="?", help="VMD file type", default="parm7"
    )
    parser.add_argument(
        "--coord_type", type=str, nargs="?", help="VMD file type", default="rst7"
    )
    args = parser.parse_args()
    template_dir = os.path.join(os.path.dirname(__file__), "..", "templates/")
    environment = Environment(loader=FileSystemLoader(template_dir))
    template = environment.get_template("inspect_structure.vmd")
    vmd_state_content = template.render(
        topology_file=args.topo,
        topology_file_type=args.topo_type,
        coord_file=args.coord,
        coord_file_type=args.coord_type,
    )
    with tempfile.NamedTemporaryFile() as tmp_file:
        with open(tmp_file.name, "w", encoding="utf-8") as f:
            f.write(vmd_state_content)
        subprocess.run(["vmd", "-e", tmp_file.name], check=False)
