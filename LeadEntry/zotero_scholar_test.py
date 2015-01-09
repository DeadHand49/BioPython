from __future__ import unicode_literals

__author__ = 'memery'

import unittest
import zotero_scholar as zs


class ZoteroEntryTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.zotero = zs.ZoteroEntry(open('ZoteroTest.csv', 'rb'), 'test')

    def test_fetch_pmid(self):
        self.assertEqual('24411336', self.zotero.read_csv())