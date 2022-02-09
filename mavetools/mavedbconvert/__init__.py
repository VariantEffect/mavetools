import os
import sys
import logging.config

__all__ = [
    "base",
    "constants",
    "empiric",
    "enrich",
    "enrich2",
    "exceptions",
    "utilities",
    "filters",
    "validators",
    "LOGGER",
]

HOMEDIR = os.path.normpath(os.path.expanduser("~/.mavedb_convert/"))
if not os.path.isdir(HOMEDIR):
    os.mkdir(HOMEDIR)  # pragma: no cover

LOGGER = "mavedbconvert"

# Initialize the logging via dictionary configuration
logging.config.dictConfig(
    {
        "version": 1,
        "disable_existing_loggers": False,
        "formatters": {
            "verbose": {"format": "[%(levelname)s] %(asctime)s %(module)s %(message)s"},
            "simple": {"format": "%(levelname)s %(message)s"},
        },
        "handlers": {
            "file": {
                "level": "WARNING",
                "class": "logging.FileHandler",
                "filename": os.path.join(HOMEDIR, "info.log"),
                "formatter": "verbose",
            },
            "console": {
                "level": "INFO",
                "class": "logging.StreamHandler",
                "stream": sys.stdout,
                "formatter": "verbose",
            },
        },
        "loggers": {
            LOGGER: {
                "handlers": ["file", "console"],
                "level": "INFO",
                "propagate": True,
            }
        },
    }
)