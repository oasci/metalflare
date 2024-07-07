import json

import matplotlib as mpl
import matplotlib.pyplot as plt


def use_mpl_rc_params(rc_json_path: str, font_dirs: tuple[str] | None = None) -> None:
    with open(rc_json_path, "r", encoding="utf-8") as f:
        rc_params = json.load(f)
    if font_dirs is not None:
        font_paths = mpl.font_manager.findSystemFonts(
            fontpaths=font_dirs, fontext="ttf"
        )
        for font_path in font_paths:
            mpl.font_manager.fontManager.addfont(font_path)
    for key, params in rc_params.items():
        plt.rc(key, **params)
