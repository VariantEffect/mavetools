import os
import shutil
import unittest
import tempfile

import pandas as pd

import logging

logging.disable(logging.CRITICAL)


__all__ = [
    "test_base",
    "test_empiric",
    "test_enrich",
    "test_enrich2",
    "test_utilities",
    "test_filters",
    "test_validators",
    "ProgramTestCase",
]


# TODO: think up a better name for this class
class ProgramTestCase(unittest.TestCase):
    def setUp(self):
        self._data_dir = tempfile.TemporaryDirectory()  # store the object
        self.data_dir = os.path.join(
            self._data_dir.name, "data"
        )  # store the directory path
        shutil.copytree(
            src=os.path.join(os.path.dirname(os.path.abspath(__file__)), "data"),
            dst=self.data_dir,
        )

    def mock_multi_sheet_excel_file(self, path, data):
        writer = pd.ExcelWriter(path, engine="xlsxwriter")
        for i, di in enumerate(data):
            df = pd.DataFrame(di)
            df.to_excel(writer, sheet_name="Sheet{}".format(i), index=False)
        writer.save()

    def tearDown(self):
        self._data_dir.cleanup()