from __future__ import unicode_literals
__author__ = 'memery'

from leadentry import Article, Author, lookup_up_title
from Bio import Entrez
from bs4 import BeautifulSoup
import unicodecsv as csv


class ZoteroEntry(object):

    def __init__(self, zotero_csv, specific_source):
        self.zotero_csv = zotero_csv
        self.specific_source = specific_source
        self.pmid_list = []
        self.read_csv()
        self.records = self.fetch_from_pubmed()
        self.articles = []
        self.parse_pubmed_soup()


    def fetch_pmid(self, title):
        self.pmid_list.append(lookup_up_title(title))

    def read_csv(self):
        with self.zotero_csv as zotero:
            reader = csv.DictReader(zotero)
            for row in reader:
                self.fetch_pmid(row['Title'])


    def fetch_from_pubmed(self):
        """Returns a list of Pubmed records based on a list of PMIDs"""
        Entrez.email = "matthew.emery@stemcell.com"
        handle = Entrez.efetch(db='pubmed', id=self.pmid_list, retmode='xml')
        return BeautifulSoup(handle.read())

    def parse_pubmed_soup(self):
        """Appends to self.list_of_articles Article objects"""
        for article in self.records('pubmedarticle'):
            self.articles.append((Article(article)))

    def write_csv(self):
        """Complete"""
        field_names = ('Publication Link', 'Publication Date', 'Article Name', 'Abstract', 'Search Term',
                       'Product Use/Assay Type', 'Product Line', 'First Name', 'Last Name', 'Company',
                       'Department', 'Email', 'Lead Source', 'Specific Lead Source')
        pass

if __name__ == '__main__':
    tester = ZoteroEntry(open('ZoteroTest.csv', 'rb'), 'Test')
    print tester.articles