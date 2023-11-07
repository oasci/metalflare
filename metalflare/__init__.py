__version__ = "0.0.0"

import os
import sys

from loguru import logger

logger.disable("metalflare")

LOG_FORMAT = "<green>{time:HH:mm:ss}</green> | " \
    "<level>{level: <8}</level> | " \
    "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>"

def enable_logging(level: int, file_path: str | None = None) -> None:
    r"""Enable logging.

    Args:
        level: Requested log level: `10` is debug, `20` is info.
        file_path: Also write logs to files here.
    """
    config = {"handlers": [{"sink": sys.stdout, "level": level, "format": LOG_FORMAT}]}
    if isinstance(file_path, str):
        config["handlers"].append(
            {"sink": file_path, "level": level, "serialize": True, "format": LOG_FORMAT}
        )
    # https://loguru.readthedocs.io/en/stable/api/logger.html#loguru._logger.Logger.configure
    logger.configure(**config)

    logger.enable("metalflare")


if os.environ.get("METALFLARE_LOG", 0):
    level = int(os.environ.get("METALFLARE_LOG_LEVEL", 20))
    file_path = os.environ.get("METALFLARE_LOG_FILE_PATH", None)
    enable_logging(level, file_path)
