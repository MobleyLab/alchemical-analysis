r"""
Alchemical Analysis: An open tool implementing some recommended practices
for analyzing alchemical free energy calculations.
"""

import sys
import os
import logging


class MyFormatter(logging.Formatter):

    err_fmt  = "ERROR: %(msg)s"
    warn_fmt = "WARNING: %(msg)s"
    info_fmt = "%(msg)s"
    dbg_fmt  = "DEBUG: %(module)s: %(lineno)d: %(msg)s"

    def __init__(self, fmt="%(levelno)s: %(msg)s"):
        logging.Formatter.__init__(self, fmt)

    def format(self, record):
        format_orig = self._fmt

        if record.levelno == logging.DEBUG:
            self._fmt = MyFormatter.dbg_fmt
        elif record.levelno == logging.WARN:
            self._fmt = MyFormatter.warn_fmt
        elif record.levelno == logging.INFO:
            self._fmt = MyFormatter.info_fmt
        elif record.levelno == logging.ERROR:
            self._fmt = MyFormatter.err_fmt

        result = logging.Formatter.format(self, record)

        self._fmt = format_orig

        return result


# FIXME: replace with log file?
logging.basicConfig(level=logging.INFO, filename=os.devnull)
logger = logging.getLogger(__name__)

hdlr = logging.StreamHandler(sys.stdout)
hdlr.setFormatter(MyFormatter())
logging.root.addHandler(hdlr)
