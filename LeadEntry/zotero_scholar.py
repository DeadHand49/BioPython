from __future__ import unicode_literals
__author__ = 'memery'

from leadentry import Article, Author, lookup_up_title
import unicodecsv as csv


class ZoteroEntry(object):

    def __init__(self, zotero_csv, specific_source):
        self.zotero_csv = zotero_csv
        self.specific_source = specific_source
        self.pmid_list = []

    def fetch_pmid(self, title):
        self.pmid_list.append(lookup_up_title(title))

    def read_csv(self):
        with self.zotero_csv as zotero:
            reader = csv.DictReader(zotero)
            for row in reader:
                self.fetch_pmid(row['Title'])


if __name__ == '__main__':
    tester = ZoteroEntry(open('ZoteroTest.csv', 'rb'), 'Test')
    tester.read_csv()
    print tester.pmid_list